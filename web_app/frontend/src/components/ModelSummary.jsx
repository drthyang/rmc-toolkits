// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useContext, useMemo, useState } from 'react';
import { describeSymmetry, toleranceLadder } from '../symmetryModel';
import { moveRatios } from '../moveStats';
import { SymTolContext } from '../symTolContext';
import { Stat, StatRail } from '../ui';
import InfoBadge from '../ui/InfoBadge';
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
        color: pct > 50 ? 'var(--on-accent)' : 'var(--text)'
    };
};

// What the atom-line parser could not use (both runtimes report the same
// counts: browserData.structureFromRmc6f / Flask /api/structure `parseReport`),
// condensed for the card; the full sentence goes in the tooltip. Null when the
// atom section parsed cleanly.
const parseSummary = (structure) => {
    const report = structure?.parseReport;
    const warning = structure?.parseWarning;
    if (!report || !warning) return null;
    const skipped = (report.invalidLines || 0) + (report.nonFiniteLines || 0);
    const accepted = (report.parsedAtoms || 0) + (report.coordsOnlyAtoms || 0);
    const missing = Number.isFinite(report.declaredAtoms) ? report.declaredAtoms - accepted : null;
    return { warning, skipped, missing, declared: report.declaredAtoms };
};

// `showSymmetry={false}` renders the model card alone — pages that want the
// structure facts without the Detected SG card (Bond Geometry) opt out, and
// the symmetry finder is skipped entirely rather than computed and hidden.
const ModelSummary = ({ structure, showSymmetry = true }) => {
    // Tolerance is shared via context (kept across page switches); fall back to
    // local state if no provider is present.
    const sharedSymTol = useContext(SymTolContext);
    const localSymTol = useState(0.2);
    const [symTol, setSymTol] = sharedSymTol ?? localSymTol;

    // Symmetry finder (FINDSYM-like): space group + Wyckoff orbits at `symTol`,
    // plus the full symmetry-vs-tolerance ladder. Runs client-side, no backend.
    const symmetry = useMemo(
        () => (showSymmetry ? describeSymmetry(structure, symTol) : null),
        [structure, symTol, showSymmetry]
    );
    const ladder = useMemo(
        () => (showSymmetry ? toleranceLadder(structure, 1.0) : []),
        [structure, showSymmetry]
    );
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
            elementEntries,
            moves: moveRatios(structure.moves, structure.totalAtoms),
            parse: parseSummary(structure)
        };
    }, [structure]);

    if (!summary) return null;

    return (
        <div className="ui-stack">
            <StatRail
                aria-label="Model information"
                headingProps={{ title: summary.sourcePath || summary.source }}
                heading={(
                    <>
                        Model information
                        <span className="ui-stat-rail__source">{summary.source}</span>
                    </>
                )}
            >
                <Stat label="Cell (Å)">
                    {summary.cellLengths.map((value) => formatNumber(value)).join(' × ')}
                </Stat>
                <Stat label="Angles">
                    {summary.angles.map((value) => `${formatNumber(value, 1)}°`).join(' · ')}
                </Stat>
                <Stat label="Supercell">
                    {summary.supercell.map((value) => formatNumber(value, 0)).join(' × ')}
                </Stat>
                {summary.elementEntries.map(({ element, count, referenceSites }) => (
                    <Stat key={element} label={element}>
                        {formatNumber(count, 0)}
                        {referenceSites > 0 && (
                            <span className="ui-stat__sub">{formatNumber(referenceSites, 0)} sites</span>
                        )}
                    </Stat>
                ))}
                <Stat label="Total atoms">
                    {formatNumber(summary.totalAtoms, 0)}
                    {summary.referenceSites > 0 && (
                        <span className="ui-stat__sub">{formatNumber(summary.referenceSites, 0)} sites</span>
                    )}
                </Stat>
                {/* Atom lines the parser could not use: unparsed layouts, non-finite
                    coordinates, or fewer atoms than the header declares (e.g. a
                    Live Data read of a file still being written). */}
                {summary.parse && (
                    <Stat role="status" label="Parse warning" ddProps={{ title: summary.parse.warning }}>
                        {summary.parse.skipped > 0
                            ? `${formatNumber(summary.parse.skipped, 0)} lines skipped`
                            : `${formatNumber(Math.abs(summary.parse.missing ?? 0), 0)} atoms ${summary.parse.missing > 0 ? 'missing' : 'extra'}`}
                        <span className="ui-stat__sub">
                            {summary.parse.declared != null
                                ? `header declares ${formatNumber(summary.parse.declared, 0)}`
                                : 'hover for details'}
                        </span>
                    </Stat>
                )}
                {/* Move counters per atom — the raw totals mean little without the
                    box size. Absent for configurations whose header omits them. */}
                {summary.moves?.generatedPerAtom !== undefined && (
                    <Stat
                        end
                        label="Generated / atom"
                        ddProps={{ title: `${formatNumber(summary.moves.generated, 0)} moves generated over ${formatNumber(summary.totalAtoms, 0)} atoms` }}
                    >
                        {formatNumber(summary.moves.generatedPerAtom, 1)}
                        <span className="ui-stat__sub">{formatNumber(summary.moves.generated, 0)} moves</span>
                    </Stat>
                )}
                {summary.moves?.acceptedPerAtom !== undefined && (
                    <Stat
                        end
                        label="Accepted / atom"
                        ddProps={{ title: `${formatNumber(summary.moves.accepted, 0)} moves accepted over ${formatNumber(summary.totalAtoms, 0)} atoms` }}
                    >
                        {formatNumber(summary.moves.acceptedPerAtom, 2)}
                        <span className="ui-stat__sub">{formatNumber(summary.moves.accepted, 0)} moves</span>
                    </Stat>
                )}
                {summary.moves?.acceptedPerGenerated !== undefined && (
                    <Stat
                        end
                        label="Accepted / generated"
                        ddProps={{ title: 'Acceptance ratio: accepted moves as a fraction of those generated' }}
                    >
                        {formatNumber(summary.moves.acceptedPerGenerated, 3)}
                        <span className="ui-stat__sub">
                            {formatNumber(summary.moves.acceptedPerGenerated * 100, 1)}% accepted
                        </span>
                    </Stat>
                )}
            </StatRail>

            {symmetry && (
                <StatRail
                    aria-label="Detected space group"
                    heading={(
                        <>
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
                                    The unit cell is the <code>.rmc6f</code> supercell divided by its supercell
                                    dimensions, but the symbol is reported in its standard setting: the finder also
                                    tries the other axis orders and cells built from the symmetry elements (a centred
                                    or primitive cell, or the true cell of a supercell). Where no standard setting is
                                    found the crystal class is shown, without a number; a symbol marked ≥ is a lower
                                    bound. Unlike FINDSYM there is no origin shift or idealized structure.
                                </p>
                            </InfoBadge>
                        </>
                    )}
                >
                    {/* A skipped (not analysed) structure has no fit -- maxResidual is
                        NaN -- so its tooltip gives the reason, never 'fits to NaN Å'. */}
                    <Stat
                        label="Space group"
                        ddProps={{
                            title: symmetry.skipped
                                ? symmetry.reason
                                : `Point group ${symmetry.pointGroup}${Number.isFinite(symmetry.maxResidual)
                                    ? ` · fits to ${symmetry.maxResidual.toFixed(3)} Å`
                                    : ''}`
                        }}
                    >
                        {symmetry.spaceGroup}
                        <span className="ui-stat__sub">
                            {symmetry.spaceGroupNumber ? `No. ${symmetry.spaceGroupNumber} · ` : ''}{symmetry.pointGroup}
                        </span>
                    </Stat>
                    <Stat label="Operations">
                        {symmetry.nSpace}
                        <span className="ui-stat__sub">symmetry ops</span>
                    </Stat>
                    {ladder.length > 0 && (
                        <Stat
                            className="sym-ladder-cell"
                            dtProps={{ className: 'sym-ladder-dt' }}
                            label={(
                                <>
                                    <span>Space group vs. tolerance</span>
                                    <span className="sym-tol-arrow" title="Bricks run from tight (left) to loose (right) atomic-position tolerance">atom pos. tol. →</span>
                                </>
                            )}
                        >
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
                        </Stat>
                    )}
                </StatRail>
            )}
        </div>
    );
};

export default ModelSummary;
