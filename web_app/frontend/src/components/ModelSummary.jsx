import React, { useMemo } from 'react';
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

const ModelSummary = ({ structure }) => {
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
        </section>
    );
};

export default ModelSummary;
