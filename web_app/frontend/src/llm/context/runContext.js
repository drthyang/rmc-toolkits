// The seam where dashboard state becomes LLM input. buildRunContext distills
// the run object the dashboard already holds (parsed plot files, the combined
// R-value history, the folded structure, the detected symmetry) into a compact
// plain object whose JSON form fits comfortably in a small local model's
// context window. Everything here is a pure function of plain data — no React,
// no network — which is what makes it the most testable part of the module.

// Rough character budget for the serialized context: ~3k chars leaves a 4k-token
// local model plenty of room for the system prompt, the question, and its answer.
export const CONTEXT_CHAR_BUDGET = 3000;

const HISTORY_POINTS = 48;

const roundSig = (value, digits = 3) => (
    Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : value
);

// Uniform-stride downsample that always keeps the first and last points, with
// values rounded to 3 significant digits — raw float noise just wastes tokens.
export const downsampleSeries = (values, maxPoints = HISTORY_POINTS) => {
    if (!Array.isArray(values) || !values.length) return [];
    if (values.length <= maxPoints) return values.map((value) => roundSig(value));
    const picked = [];
    const step = (values.length - 1) / (maxPoints - 1);
    for (let i = 0; i < maxPoints; i += 1) {
        picked.push(roundSig(values[Math.round(i * step)]));
    }
    return picked;
};

// Least-squares slope (per step) over the tail of the series — the last 20% of
// points, at least 5 — so "still improving vs flat" reflects recent behavior.
export const recentSlope = (values) => {
    if (!Array.isArray(values) || values.length < 2) return 0;
    const count = Math.max(Math.min(values.length, 5), Math.ceil(values.length * 0.2));
    const tail = values.slice(-count);
    const n = tail.length;
    const meanX = (n - 1) / 2;
    const meanY = tail.reduce((sum, value) => sum + value, 0) / n;
    let covariance = 0;
    let variance = 0;
    tail.forEach((value, index) => {
        covariance += (index - meanX) * (value - meanY);
        variance += (index - meanX) ** 2;
    });
    return variance ? covariance / variance : 0;
};

// Summary statistics computed on the FULL series (before downsampling), so
// nothing important is lost to sampling.
export const seriesStats = (values) => {
    if (!Array.isArray(values) || !values.length) return null;
    return {
        nSteps: values.length,
        first: values[0],
        last: values[values.length - 1],
        min: Math.min(...values),
        max: Math.max(...values),
        recentSlopePerStep: recentSlope(values)
    };
};

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));

const angleBetween = (a, b) => {
    const denominator = Math.max(vectorLength(a) * vectorLength(b), 1e-12);
    const cosine = a.reduce((sum, value, index) => sum + value * b[index], 0) / denominator;
    return Math.acos(Math.max(-1, Math.min(1, cosine))) * (180 / Math.PI);
};

// Cell lengths/angles are recomputed from the lattice vectors with the same
// math as ModelSummary.jsx — duplicated deliberately: importing from a host
// component would break this module's import boundary.
const structureContext = (structure, symmetry) => {
    if (!structure?.latticeVectors || !structure?.supercell) return null;
    const cell = {
        cell_A: structure.latticeVectors.map((vector, index) => (
            roundSig(vectorLength(vector) / Math.max(structure.supercell[index], 1), 4)
        )),
        angles_deg: [
            angleBetween(structure.latticeVectors[1], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[1])
        ].map((angle) => roundSig(angle, 4)),
        supercell: structure.supercell,
        total_atoms: structure.totalAtoms
    };
    if (structure.elementCounts && Object.keys(structure.elementCounts).length) {
        cell.elements = structure.elementCounts;
    }
    if (symmetry?.spaceGroup) {
        cell.space_group = `${symmetry.spaceGroup}`
            + `${symmetry.spaceGroupNumber ? ` (No. ${symmetry.spaceGroupNumber})` : ''}`;
    }
    return cell;
};

const datasetContext = (plotFiles) => (
    (plotFiles || [])
        .filter((file) => file.plotKind && file.plotKind !== 'r_value')
        .map((file) => {
            const entry = { kind: file.plotKind };
            const title = file.plotData?.title;
            if (title && title !== file.plotKind) entry.title = title;
            const rwpValue = file.plotData?.metrics?.rwp;
            if (Number.isFinite(rwpValue)) entry.rwp = roundSig(rwpValue);
            return entry;
        })
);

const convergenceContext = (rValueFile, historyPoints) => {
    const history = rValueFile?.plotData?.series?.[0]?.y;
    const stats = seriesStats(history);
    if (!stats) return null;
    const convergence = {
        // The dashboard stores ln(chi^2), not raw chi^2 (browserData.js applies
        // Math.log when parsing the .log files) — say so, or the model will
        // misread the magnitudes.
        quantity: 'ln of chi^2 goodness metric (natural log; lower is better)',
        n_steps: stats.nSteps,
        first: roundSig(stats.first),
        last: roundSig(stats.last),
        min: roundSig(stats.min),
        max: roundSig(stats.max),
        recent_slope_per_step: roundSig(stats.recentSlopePerStep, 2),
        history: downsampleSeries(history, historyPoints)
    };
    const finalChi = rValueFile?.plotData?.metrics?.final_chi_r;
    if (Number.isFinite(finalChi)) convergence.final_chi_squared = roundSig(finalChi);
    return convergence;
};

// Distill dashboard props into the compact context object sent to the model.
// Every input is optional: missing pieces (no structure file, no logs) simply
// omit their section, and the prompt templates handle absence explicitly.
export const buildRunContext = ({
    runName,
    plotFiles = [],
    rValueFile = null,
    structure = null,
    symmetry = null,
    liveData = false
} = {}) => {
    const context = {
        run: runName || 'unnamed run',
        live_mode: Boolean(liveData)
    };
    const structureInfo = structureContext(structure, symmetry);
    if (structureInfo) context.structure = structureInfo;
    const datasets = datasetContext(plotFiles);
    if (datasets.length) context.datasets = datasets;
    const convergence = convergenceContext(rValueFile, HISTORY_POINTS);
    if (convergence) context.convergence = convergence;
    return context;
};

// Serialize with the character budget enforced: first shrink the history, then
// truncate the dataset list (recording how many were dropped) — never silently.
export const contextToJson = (context, budget = CONTEXT_CHAR_BUDGET) => {
    let current = context;
    let json = JSON.stringify(current, null, 1);
    for (const historyPoints of [24, 12]) {
        if (json.length <= budget) return json;
        if (current.convergence?.history?.length > historyPoints) {
            current = {
                ...current,
                convergence: {
                    ...current.convergence,
                    history: downsampleSeries(current.convergence.history, historyPoints)
                }
            };
            json = JSON.stringify(current, null, 1);
        }
    }
    if (json.length > budget && current.datasets?.length > 12) {
        current = {
            ...current,
            datasets: current.datasets.slice(0, 12),
            datasets_omitted: current.datasets.length - 12
        };
        json = JSON.stringify(current, null, 1);
    }
    return json;
};
