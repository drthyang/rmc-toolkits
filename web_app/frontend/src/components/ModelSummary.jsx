import React, { useMemo, useState } from 'react';
import { describeSymmetry, toleranceLadder, orbitLabel } from '../symmetryModel';
import './ModelSummary.css';

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));

const angleBetween = (a, b) => {
    const denominator = Math.max(vectorLength(a) * vectorLength(b), 1e-12);
    const cosine = a.reduce((sum, value, index) => sum + value * b[index], 0) / denominator;
    const clamped = Math.max(-1, Math.min(1, cosine));
    return Math.acos(clamped) * (180 / Math.PI);
};

const formatNumber = (value, digits = 3) => Number(value).toLocaleString(undefined, {
    maximumFractionDigits: digits
});

// Ladder brick style: more operations → deeper accent fill (theme-aware via
// color-mix over the panel), with a legible label colour for the fill darkness.
const brickStyle = (nSpace, maxOps) => {
    const level = maxOps > 1 ? Math.log(nSpace) / Math.log(maxOps) : 0;
    const pct = Math.round(12 + level * 74);
    return {
        background: `color-mix(in srgb, var(--accent) ${pct}%, var(--panel-raised))`,
        color: pct > 50 ? '#fff' : 'var(--text)'
    };
};

const ModelSummary = ({ structure }) => {
    const [symTol, setSymTol] = useState(0.2);   // Å tolerance for symmetry detection

    // Symmetry finder (FINDSYM-like): space group + Wyckoff orbits at `symTol`,
    // plus the full symmetry-vs-tolerance ladder. Runs client-side, no backend.
    const symmetry = useMemo(() => describeSymmetry(structure, symTol), [structure, symTol]);
    const ladder = useMemo(() => toleranceLadder(structure, 1.0), [structure]);
    const maxOps = ladder.length ? Math.max(...ladder.map((b) => b.nSpace)) : 1;

    // Brick widths are NOT the raw tolerance range — the full-symmetry rung holds
    // over most of the axis and would dominate. Cap the widest rung at ~1/3 and
    // split the rest evenly, so the progression reads clearly.
    const widestBrick = ladder.reduce((best, b, i) => ((b.to - b.from) > (ladder[best].to - ladder[best].from) ? i : best), 0);
    const brickWidth = (i) => (ladder.length <= 1 ? 100 : i === widestBrick ? 34 : 66 / (ladder.length - 1));

    const summary = useMemo(() => {
        if (!structure?.latticeVectors || !structure?.supercell) return null;

        const boxLengths = structure.latticeVectors.map(vectorLength);
        const cellLengths = boxLengths.map((length, index) => length / Math.max(structure.supercell[index], 1));
        const angles = [
            angleBetween(structure.latticeVectors[1], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[1])
        ];
        const elementEntries = Object.entries(structure.elementCounts || {})
            .sort(([a], [b]) => a.localeCompare(b))
            .map(([element, count]) => ({
                element,
                count,
                referenceSites: structure.atomIndices?.[element]?.length || 0
            }));
        const referenceSites = Object.values(structure.atomIndices || {}).reduce((sum, indices) => sum + indices.length, 0);

        return {
            source: structure.source?.split('/').pop() || 'Structure file',
            totalAtoms: structure.totalAtoms,
            referenceSites,
            supercell: structure.supercell,
            cellLengths,
            angles,
            elementEntries
        };
    }, [structure]);

    if (!summary) return null;

    return (
        <section className="model-summary" aria-label="Model information">
            <h2 className="model-summary-title" title={summary.source}>Model information</h2>
            <dl className="model-stats">
                <div className="model-stat">
                    <dt>Cell (Å)</dt>
                    <dd>{summary.cellLengths.map((value) => formatNumber(value)).join(' × ')}</dd>
                </div>
                <div className="model-stat">
                    <dt>Angles</dt>
                    <dd>{summary.angles.map((value) => `${formatNumber(value, 1)}°`).join(' · ')}</dd>
                </div>
                <div className="model-stat">
                    <dt>Supercell</dt>
                    <dd>{summary.supercell.map((value) => formatNumber(value, 0)).join(' × ')}</dd>
                </div>
                {summary.elementEntries.map(({ element, count, referenceSites }) => (
                    <div className="model-stat" key={element}>
                        <dt>{element}</dt>
                        <dd>
                            {formatNumber(count, 0)}
                            <span className="model-stat-sub">{formatNumber(referenceSites, 0)} sites</span>
                        </dd>
                    </div>
                ))}
                <div className="model-stat">
                    <dt>Total atoms</dt>
                    <dd>
                        {formatNumber(summary.totalAtoms, 0)}
                        <span className="model-stat-sub">{formatNumber(summary.referenceSites, 0)} sites</span>
                    </dd>
                </div>
            </dl>

            {symmetry && (
                <div className="model-symmetry">
                    <h3 className="model-summary-title model-symmetry-title">Symmetry Analysis</h3>
                    <div className="model-symmetry-body">
                        <div className="sg-line">
                            <span className="sg-symbol" title={`Point group ${symmetry.pointGroup} · ${symmetry.nSpace} operations · fits to ${symmetry.maxResidual.toFixed(3)} Å`}>
                                {symmetry.spaceGroup}
                            </span>
                            {symmetry.spaceGroupNumber && <span className="sg-number">No. {symmetry.spaceGroupNumber}</span>}
                            <span className="sg-meta">{symmetry.pointGroup} · {symmetry.nSpace} ops</span>
                        </div>

                        {ladder.length > 0 && (
                            <div className="sym-ladder" role="group" aria-label="Space group vs. tolerance — click to select">
                                {ladder.map((b, i) => {
                                    const active = symTol >= b.from && symTol < b.to;
                                    return (
                                        <button
                                            key={i}
                                            type="button"
                                            className={`sym-brick${active ? ' is-active' : ''}`}
                                            style={{ width: `${brickWidth(i)}%`, ...brickStyle(b.nSpace, maxOps) }}
                                            title={`${b.spaceGroup}${b.spaceGroupNumber ? ` (No. ${b.spaceGroupNumber})` : ''} · holds ${b.from.toFixed(2)}–${b.to.toFixed(2)} Å · ${b.nSpace} ops — click to select`}
                                            onClick={() => setSymTol(Math.max(0.02, Math.min(1, (b.from + b.to) / 2)))}
                                        >
                                            <span className="sym-brick-label">{b.spaceGroup}</span>
                                        </button>
                                    );
                                })}
                            </div>
                        )}

                        {symmetry.orbits.length > 0 && (
                            <div className="sym-orbits">
                                {symmetry.orbits.map((o, i) => (
                                    <span className="sym-orbit" key={i} title={`site symmetry ${o.site}`}>
                                        <span className="sym-orbit-el">{o.element}</span>
                                        {orbitLabel(o)}
                                    </span>
                                ))}
                            </div>
                        )}
                    </div>
                </div>
            )}
        </section>
    );
};

export default ModelSummary;
