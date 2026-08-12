// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useContext, useMemo, useState } from 'react';
import { describeSymmetry, toleranceLadder } from '../symmetryModel';
import { SymTolContext } from '../symTolContext';
import InfoBadge from './InfoBadge';
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
    // Tolerance is shared via context (kept across page switches); fall back to
    // local state if no provider is present.
    const sharedSymTol = useContext(SymTolContext);
    const localSymTol = useState(0.2);
    const [symTol, setSymTol] = sharedSymTol ?? localSymTol;

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
            sourcePath: structure.source,
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
        <div className="model-cards">
            <section className="model-summary" aria-label="Model information">
                <h2 className="model-summary-title" title={summary.sourcePath || summary.source}>
                    Model information
                    <span className="model-summary-source">{summary.source}</span>
                </h2>
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
                                {referenceSites > 0 && (
                                    <span className="model-stat-sub">{formatNumber(referenceSites, 0)} sites</span>
                                )}
                            </dd>
                        </div>
                    ))}
                    <div className="model-stat">
                        <dt>Total atoms</dt>
                        <dd>
                            {formatNumber(summary.totalAtoms, 0)}
                            {summary.referenceSites > 0 && (
                                <span className="model-stat-sub">{formatNumber(summary.referenceSites, 0)} sites</span>
                            )}
                        </dd>
                    </div>
                </dl>
            </section>

            {symmetry && (
                <section className="model-summary model-symmetry" aria-label="Detected space group">
                    <h2 className="model-summary-title">
                        Detected SG
                        <InfoBadge label="How space-group detection works">
                            <p>
                                The average (reference) site positions from the <code>.rmc6f</code> model
                                are tested for the symmetry operations {'{R | t}'} that map them onto
                                themselves within a position tolerance (Å). Each matching operation is
                                split into a rotation and a translation part, so screw axes and glide
                                planes are named as such, giving the Hermann–Mauguin symbol and number.
                            </p>
                            <p>
                                It runs entirely in your browser — no fitting or external service.
                                Loosening the tolerance admits more operations, so higher symmetry
                                appears; the ladder shows which space group holds over each tolerance range.
                            </p>
                            <p>
                                Unlike FINDSYM this does not search for a smaller or differently-shaped
                                cell: the unit cell is taken from the <code>.rmc6f</code> supercell as
                                given, and the symbol is reported in that setting (up to an axis
                                permutation), with no origin shift or idealized structure.
                            </p>
                        </InfoBadge>
                    </h2>
                    <dl className="model-stats">
                        <div className="model-stat">
                            <dt>Space group</dt>
                            <dd title={`Point group ${symmetry.pointGroup} · fits to ${symmetry.maxResidual.toFixed(3)} Å`}>
                                {symmetry.spaceGroup}
                                <span className="model-stat-sub">
                                    {symmetry.spaceGroupNumber ? `No. ${symmetry.spaceGroupNumber} · ` : ''}{symmetry.pointGroup}
                                </span>
                            </dd>
                        </div>
                        <div className="model-stat">
                            <dt>Operations</dt>
                            <dd>
                                {symmetry.nSpace}
                                <span className="model-stat-sub">symmetry ops</span>
                            </dd>
                        </div>
                        {ladder.length > 0 && (
                            <div className="model-stat sym-ladder-cell">
                                <dt className="sym-ladder-dt">
                                    <span>Space group vs. tolerance</span>
                                    <span className="sym-tol-arrow" title="Bricks run from tight (left) to loose (right) atomic-position tolerance">atom pos. tol. →</span>
                                </dt>
                                <dd>
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
                                                    onClick={() => setSymTol((b.from + b.to) / 2)}
                                                >
                                                    <span className="sym-brick-label">{b.spaceGroup}</span>
                                                </button>
                                            );
                                        })}
                                    </div>
                                </dd>
                            </div>
                        )}
                    </dl>
                </section>
            )}
        </div>
    );
};

export default ModelSummary;
