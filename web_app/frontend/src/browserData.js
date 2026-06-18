const SUPPORTED_NAMES = new Set(['scale_ft.gr', 'scale_ft.sq', 'scale_ft_rmc.fq', 'stog_input.dat']);

export const isStaticMode = () => {
    if (import.meta.env.VITE_STATIC_MODE === 'true') return true;
    return window.location.hostname.endsWith('github.io');
};

export const detectPlotKind = (name) => {
    if (name.endsWith('_FT_XFQ1.csv')) return 'xpdf';
    if (name.includes('PDF') && name.endsWith('.csv')) {
        return name.includes('PDFpartials') ? 'pdf_partials' : 'npdf';
    }
    if (name.endsWith('_FQ1.csv')) return 'xray_sq';
    if (name.endsWith('_SQ1.csv')) return 'neutron_sq';
    if (name.endsWith('_bragg.csv')) return 'bragg';
    if (name.endsWith('.log')) return 'r_value';
    if (['scale_ft.gr', 'scale_ft.sq', 'scale_ft_rmc.fq'].includes(name)) return 'stog';
    return null;
};

const isSupportedFile = (name) => (
    name.endsWith('.csv')
    || name.endsWith('.log')
    || name.endsWith('.rmc6f')
    || name.startsWith('Frac')
    || SUPPORTED_NAMES.has(name)
);

const basename = (path) => path.split('/').pop() || path;

const parseNumberRows = (lines, startIndex = 0, separator = /\s+/) => {
    const rows = [];
    lines.slice(startIndex).forEach((line) => {
        const values = line.trim().split(separator).filter(Boolean);
        if (!values.length) return;
        const parsed = values.map(Number);
        if (parsed.every(Number.isFinite)) rows.push(parsed);
    });
    return rows;
};

const transpose = (rows) => rows[0].map((_, index) => rows.map((row) => row[index]));

const readRmcCsv = (text, name) => {
    const lines = text.split(/\r?\n/).filter((line) => line.trim());
    if (!lines.length) throw new Error(`${name} is empty`);
    const labels = lines[0].split(',').map((label) => label.trim());
    const rows = lines.slice(1).map((line, index) => {
        const values = line.split(',').map((value) => value.trim()).filter(Boolean);
        if (!values.length) return null;
        if (values.length !== labels.length) {
            throw new Error(`${name} line ${index + 2} has ${values.length} values; expected ${labels.length}`);
        }
        return values.map(Number);
    }).filter(Boolean);
    if (!rows.length) throw new Error(`${name} does not contain numeric rows`);
    return { labels, data: transpose(rows) };
};

const readChi = (text) => {
    const chiR = [];
    text.split(/\r?\n/).slice(2).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (parts.length >= 2) {
            const value = Number(parts[parts.length - 1]);
            if (Number.isFinite(value)) chiR.push(value);
        }
    });
    return chiR;
};

const readStog = (text, name) => {
    const rows = parseNumberRows(text.split(/\r?\n/), 2);
    if (!rows.length) throw new Error(`${name} does not contain STOG numeric rows`);
    return transpose(rows);
};

const rwp = (x, observed, fitted) => {
    void x;
    const denom = observed.reduce((sum, value) => sum + value * value, 0);
    if (!denom) return 0;
    const residual = fitted.reduce((sum, value, index) => sum + (value - observed[index]) ** 2, 0);
    return Math.sqrt(residual / denom);
};

const pdfIndex = (name) => {
    const match = name.match(/PDF(\d+)\.csv$/);
    return match ? Number(match[1]) : 0;
};

const cleanAxisLabel = (label) => {
    const normalized = label.trim();
    if (normalized === 'Q') return 'Q (Å^{-1})';
    if (normalized === 'r' || normalized === 'R') return 'r (Å)';
    return normalized
        .replace('(A^-1)', '(Å^{-1})')
        .replace('(A^{-1})', '(Å^{-1})')
        .replace('(A)', '(Å)');
};

export const plotMetadataFromFile = (file) => {
    const kind = file.plotKind;
    if (kind === 'xpdf') return { kind, title: 'xPDF', metrics: file.plotData?.metrics || {} };
    if (kind === 'npdf') return { kind, title: file.name.replace(/\.[^.]+$/, '').split('_').pop(), metrics: file.plotData?.metrics || {} };
    if (kind === 'pdf_partials') return { kind, title: file.name.replace(/\.[^.]+$/, '').split('_').pop(), metrics: file.plotData?.metrics || {} };
    if (kind === 'xray_sq') return { kind, title: 'S(Q) (x-ray)', metrics: file.plotData?.metrics || {} };
    if (kind === 'neutron_sq') return { kind, title: 'S(Q) (neutron)', metrics: file.plotData?.metrics || {} };
    if (kind === 'bragg') return { kind, title: 'BRAGG', metrics: file.plotData?.metrics || {} };
    if (kind === 'r_value') return { kind, title: 'R-value', metrics: file.plotData?.metrics || {} };
    if (kind === 'stog') return { kind, title: file.name, metrics: file.plotData?.metrics || {} };
    return null;
};

export const plotDataFromText = (file) => {
    const kind = file.plotKind;
    if (!kind) return null;

    if (kind === 'r_value') {
        const yValues = readChi(file.text);
        if (!yValues.length) throw new Error(`${file.name} does not contain chi values`);
        return {
            kind,
            title: 'R-value',
            metrics: { final_chi_r: yValues[yValues.length - 1] },
            xLabel: 'Time steps',
            yLabel: 'log(χ)',
            series: [{
                label: 'R',
                x: yValues.map((_, index) => index),
                y: yValues.map((value) => Math.log(Math.max(value, 1e-12)))
            }]
        };
    }

    if (kind === 'stog') {
        const data = readStog(file.text, file.name);
        return {
            kind,
            title: file.name,
            metrics: {},
            xLabel: file.name.endsWith('.gr') ? 'r (Å)' : 'Q (Å^{-1})',
            yLabel: file.name.endsWith('.gr') ? 'G(r)' : 'S(Q)',
            series: [{ label: file.name, x: data[0], y: data[1] }]
        };
    }

    const csv = readRmcCsv(file.text, file.name);
    const metrics = {};
    if (['xpdf', 'npdf', 'xray_sq', 'neutron_sq', 'bragg'].includes(kind) && csv.data.length >= 3) {
        metrics.rwp = rwp(csv.data[0], csv.data[1], csv.data[2]);
    }
    if (kind === 'npdf') metrics.pdf_index = pdfIndex(file.name);

    let xLabel = csv.labels[0] || 'x';
    if (['xpdf', 'npdf', 'pdf_partials'].includes(kind)) xLabel = 'r (Å)';
    else if (kind === 'bragg') xLabel = 'Q (Å^{-1})';
    else xLabel = cleanAxisLabel(xLabel);

    const title = plotMetadataFromFile({ ...file, plotData: { metrics } })?.title || file.name;
    return {
        kind,
        title,
        metrics,
        xLabel,
        yLabel: 'data',
        series: csv.labels.slice(1).map((label, index) => ({
            label: label.trim() || `Series ${index + 1}`,
            x: csv.data[0],
            y: csv.data[index + 1]
        })).filter((series) => series.y)
    };
};

const readCellVectors = (text) => {
    const lines = text.split(/\r?\n/);
    let latticeVectors = null;
    let supercell = null;
    lines.forEach((line, index) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Supercell') supercell = parts.slice(-3).map(Number);
        if (parts[0] === 'Lattice') {
            latticeVectors = [lines[index + 1], lines[index + 2], lines[index + 3]]
                .map((row) => row.trim().split(/\s+/).map(Number));
        }
    });
    if (!latticeVectors || !supercell) throw new Error('Missing lattice or supercell metadata');
    return { latticeVectors, supercell };
};

export const structureFromRmc6f = (file, maxPoints = 100) => {
    const { latticeVectors, supercell } = readCellVectors(file.text);
    const counts = {};
    const atomIndices = {};
    const atoms = [];
    let inAtoms = false;
    file.text.split(/\r?\n/).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Atoms:') {
            inAtoms = true;
            return;
        }
        if (!inAtoms || parts.length < 10) return;
        const referenceNumber = Number(parts[6]);
        const element = parts[1];
        const coords = parts.slice(3, 6).map(Number);
        const cellIndices = parts.slice(7, 10).map(Number);
        if (![referenceNumber, ...coords, ...cellIndices].every(Number.isFinite)) return;
        counts[element] = (counts[element] || 0) + 1;
        if (!atomIndices[element]) atomIndices[element] = new Set();
        atomIndices[element].add(referenceNumber);
        atoms.push({ element, referenceNumber, coords, cellIndices });
    });

    const stride = Math.max(1, Math.ceil(atoms.length / maxPoints));
    const points = atoms.filter((_, index) => index % stride === 0).slice(0, maxPoints).map((atom) => {
        const reduced = atom.coords.map((value, index) => value - atom.cellIndices[index] / supercell[index]);
        const unitCell = reduced.map((value, index) => ((value * supercell[index]) % 1 + 1) % 1);
        return {
            element: atom.element,
            referenceNumber: atom.referenceNumber,
            boxX: atom.coords[0],
            boxY: atom.coords[1],
            boxZ: atom.coords[2],
            x: unitCell[0],
            y: unitCell[1],
            z: unitCell[2]
        };
    });

    return {
        source: file.path,
        totalAtoms: atoms.length,
        sampledAtoms: points.length,
        sampleStride: stride,
        elements: Object.keys(counts).sort(),
        elementCounts: counts,
        atomIndices: Object.fromEntries(Object.entries(atomIndices).map(([key, value]) => [key, [...value].sort((a, b) => a - b)])),
        supercell,
        latticeVectors,
        points
    };
};

export const buildLocalRun = async (fileList) => {
    const files = await Promise.all([...fileList]
        .filter((file) => isSupportedFile(file.name))
        .map(async (file) => {
            const path = file.webkitRelativePath || file.name;
            return {
                name: basename(path),
                path,
                type: 'file',
                plotKind: detectPlotKind(basename(path)),
                text: await file.text()
            };
        }));
    if (!files.length) {
        throw new Error('No supported RMCprofile files were selected');
    }

    const enrichedFiles = files.map((file) => {
        if (!file.plotKind) return file;
        try {
            return { ...file, plotData: plotDataFromText(file) };
        } catch (error) {
            return { ...file, parseError: error.message };
        }
    });

    const rmc6f = enrichedFiles.find((file) => file.name.endsWith('.rmc6f'));
    const directoryRoot = enrichedFiles
        .map((file) => file.path)
        .find((path) => path.includes('/'))
        ?.split('/')[0];

    return {
        name: directoryRoot || 'Local files',
        files: enrichedFiles,
        structure: null,
        structureFile: rmc6f || null,
        structureError: rmc6f ? 'Open KDE / 3D to load structure data' : 'No model structure detected'
    };
};
