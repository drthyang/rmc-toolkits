// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { saveSvgFigure } from '../figureExport';
import SaveMenu from './SaveMenu';
import './InteractivePlot.css';

const palette = ['#1f6fd6', '#e8590c', '#099268', '#d6336c', '#6741d9', '#66a80f', '#0c8599', '#e67700'];

const CHART_SAVE_OPTIONS = [
    { id: 'png', label: 'PNG image', hint: '.png' },
    { id: 'svg', label: 'SVG vector', hint: '.svg' },
];

// Series whose label mentions "exp" hold measured data (Experiment, F(Q)_Expt, X_ray_exp_renorm, ...).
const isExperimental = (label) => /exp/i.test(label);

const MARKER_RADIUS = 2.8;

const formatNumber = (value) => {
    const abs = Math.abs(value);
    if (abs >= 1e6 || (abs > 0 && abs < 0.01)) return value.toExponential(2);
    if (abs >= 1000) return Math.round(value).toLocaleString();
    return value.toPrecision(4);
};

const SUPERSCRIPTS = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };

// Render `Q (Å^{-1})` style labels as plain text with unicode superscripts.
const labelToText = (label) => label.replace(/\^\{([^}]+)\}/g, (_, exponent) =>
    exponent.split('').map((ch) => SUPERSCRIPTS[ch] || ch).join('')
);

// Decimal places needed so consecutive ticks of `step` stay distinct.
const decimalsForStep = (step) => {
    if (!Number.isFinite(step) || step <= 0) return 2;
    const exponent = Math.floor(Math.log10(step));
    return exponent >= 0 ? 0 : Math.min(6, -exponent);
};

// Compact tick labels: 12.5k / 3.2M instead of 1.25e4, exponents only at extremes.
// The unit (k/M) follows the axis maximum so every tick on an axis shares it.
const formatTick = (value, step, axisMax) => {
    if (Math.abs(value) < 1e-12) return '0';
    const magnitude = Math.max(Math.abs(axisMax ?? value), Math.abs(value));
    if (magnitude >= 1e9 || magnitude < 1e-4) return value.toExponential(1).replace('e+', 'e');
    if (magnitude >= 1e6) return `${(value / 1e6).toFixed(decimalsForStep(step / 1e6))}M`;
    if (magnitude >= 1e4) return `${(value / 1e3).toFixed(decimalsForStep(step / 1e3))}k`;
    return value.toFixed(decimalsForStep(step));
};

const niceDomain = (values) => {
    const finite = values.filter(Number.isFinite);
    if (!finite.length) return [0, 1];
    let min = Math.min(...finite);
    let max = Math.max(...finite);
    if (min === max) {
        min -= 1;
        max += 1;
    }
    const pad = (max - min) * 0.05;
    return [min - pad, max + pad];
};

// Round ticks at 1/2/5 x 10^n steps, clipped to the domain.
const niceTicks = (domain, count = 5) => {
    const [min, max] = domain;
    const span = max - min;
    if (!Number.isFinite(span) || span <= 0) return { ticks: [min], step: 1 };
    const raw = span / Math.max(1, count);
    const magnitude = 10 ** Math.floor(Math.log10(raw));
    const normalized = raw / magnitude;
    const multiplier = normalized < 1.5 ? 1 : normalized < 3 ? 2 : normalized < 7 ? 5 : 10;
    const step = multiplier * magnitude;
    const start = Math.ceil(min / step) * step;
    const ticks = [];
    for (let value = start; value <= max + step * 1e-6; value += step) {
        ticks.push(Math.abs(value) < step * 1e-9 ? 0 : value);
    }
    return { ticks, step };
};

const AxisLabel = ({ label, x, y, textAnchor = 'middle', rotate = false }) => {
    const superscriptMatch = label.match(/^(.*)\^\{([^}]+)\}(.*)$/);
    const transform = rotate ? `rotate(-90 ${x} ${y})` : undefined;
    if (!superscriptMatch) {
        return <text className="axis-label" x={x} y={y} textAnchor={textAnchor} transform={transform}>{label}</text>;
    }
    const [, before, superscript, after] = superscriptMatch;
    return (
        <text className="axis-label" x={x} y={y} textAnchor={textAnchor} transform={transform}>
            <tspan>{before}</tspan>
            <tspan baselineShift="super" fontSize="70%">{superscript}</tspan>
            <tspan>{after}</tspan>
        </text>
    );
};

const InteractivePlot = ({ file, variant, plotData, refreshKey }) => {
    const wide = variant === 'wide';
    const [plot, setPlot] = useState(null);
    const [error, setError] = useState(null);
    const [hidden, setHidden] = useState(() => new Set());
    const [xDomain, setXDomain] = useState(null);
    const [hover, setHover] = useState(null);
    const [drag, setDrag] = useState(null);
    const svgRef = useRef(null);
    const loadedPathRef = useRef(file.path);
    const effectivePlot = plotData || plot;

    useEffect(() => {
        if (plotData) {
            return;
        }

        const fetchData = async () => {
            setError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/plot/data`, {
                    params: { path: file.path }
                });
                if (loadedPathRef.current !== file.path) {
                    setHidden(new Set());
                    setXDomain(null);
                    setHover(null);
                    setDrag(null);
                    loadedPathRef.current = file.path;
                }
                setPlot(response.data);
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to load interactive plot');
            }
        };

        fetchData();
    }, [file.path, plotData, refreshKey]);

    // Measured data first (drawn underneath, first palette color); the
    // calculated curve follows and is drawn on top of the hollow markers.
    const orderedSeries = useMemo(() => {
        const series = effectivePlot?.series || [];
        const experimental = series.filter((entry) => isExperimental(entry.label));
        const calculated = series.filter((entry) => !isExperimental(entry.label));
        const paired = experimental.length > 0 && calculated.length > 0;
        const reordered = paired ? [...experimental, ...calculated] : series;
        return reordered.map((entry, index) => ({
            ...entry,
            marker: paired && isExperimental(entry.label),
            color: palette[index % palette.length]
        }));
    }, [effectivePlot]);

    const visibleSeries = useMemo(() => {
        return orderedSeries.filter((series) => !hidden.has(series.label));
    }, [orderedSeries, hidden]);

    const domains = useMemo(() => {
        const allX = visibleSeries.flatMap((series) => series.x);
        const baseX = niceDomain(allX);
        const currentX = xDomain || baseX;
        const allY = visibleSeries.flatMap((series) =>
            series.y.filter((_, index) => series.x[index] >= currentX[0] && series.x[index] <= currentX[1])
        );
        return { x: currentX, y: niceDomain(allY.length ? allY : visibleSeries.flatMap((series) => series.y)), baseX };
    }, [visibleSeries, xDomain]);

    // 8:5 (golden-ish) for grid cards, a slim strip for the wide variant.
    const view = wide
        ? { width: 1440, height: 320, left: 64, right: 20, top: 18, bottom: 58 }
        : { width: 720, height: 450, left: 60, right: 18, top: 16, bottom: 58 };
    const plotWidth = view.width - view.left - view.right;
    const plotHeight = view.height - view.top - view.bottom;

    const xScale = (x) => view.left + ((x - domains.x[0]) / (domains.x[1] - domains.x[0] || 1)) * plotWidth;
    const yScale = (y) => view.top + plotHeight - ((y - domains.y[0]) / (domains.y[1] - domains.y[0] || 1)) * plotHeight;
    const xInvert = (px) => domains.x[0] + ((px - view.left) / plotWidth) * (domains.x[1] - domains.x[0]);

    const yTicks = niceTicks(domains.y, wide ? 4 : 6);
    const xTicks = niceTicks(domains.x, wide ? 11 : 7);
    const yAxisMax = Math.max(...yTicks.ticks.map(Math.abs), 0);
    const xAxisMax = Math.max(...xTicks.ticks.map(Math.abs), 0);

    // Geometry is memoized so hover re-renders skip rebuilding the (large)
    // marker paths; only domain or series changes recompute them.
    const seriesShapes = useMemo(() => {
        const sx = (x) => view.left + ((x - domains.x[0]) / (domains.x[1] - domains.x[0] || 1)) * plotWidth;
        const sy = (y) => view.top + plotHeight - ((y - domains.y[0]) / (domains.y[1] - domains.y[0] || 1)) * plotHeight;
        return visibleSeries.map((series) => {
            const points = [];
            series.x.forEach((x, index) => {
                if (x < domains.x[0] || x > domains.x[1]) return;
                const y = series.y[index];
                if (!Number.isFinite(x) || !Number.isFinite(y)) return;
                points.push([sx(x), sy(y)]);
            });
            if (series.marker) {
                const r = MARKER_RADIUS;
                const commands = points.map(([cx, cy]) =>
                    `M ${(cx - r).toFixed(2)} ${cy.toFixed(2)} a ${r} ${r} 0 1 0 ${r * 2} 0 a ${r} ${r} 0 1 0 ${-r * 2} 0`
                );
                return { label: series.label, color: series.color, marker: true, d: commands.join(' ') };
            }
            const d = points.map(([px, py], index) => `${index ? 'L' : 'M'} ${px.toFixed(2)} ${py.toFixed(2)}`).join(' ');
            return { label: series.label, color: series.color, marker: false, d };
        });
    }, [visibleSeries, domains, view.left, view.top, plotWidth, plotHeight]);

    const pointerToViewX = (event) => {
        const svg = svgRef.current;
        if (!svg) return view.left;
        const transform = svg.getScreenCTM();
        if (!transform) return view.left;
        const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(transform.inverse());
        return point.x;
    };

    const clampPlotX = (x) => Math.max(view.left, Math.min(view.width - view.right, x));

    const nearestHover = (event) => {
        if (!effectivePlot || !visibleSeries.length) return;
        const x = pointerToViewX(event);
        if (x < view.left || x > view.width - view.right) {
            setHover(null);
            return;
        }
        const dataX = xInvert(x);
        const values = visibleSeries.map((series) => {
            let best = 0;
            let bestDistance = Infinity;
            series.x.forEach((value, pointIndex) => {
                const distance = Math.abs(value - dataX);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = pointIndex;
                }
            });
            return {
                label: series.label,
                color: series.color,
                x: series.x[best],
                y: series.y[best],
                cx: xScale(series.x[best]),
                cy: yScale(series.y[best])
            };
        });
        setHover({ x: values[0]?.x ?? dataX, px: xScale(values[0]?.x ?? dataX), values });
    };

    const startDrag = (event) => {
        const x = pointerToViewX(event);
        if (x < view.left || x > view.width - view.right) return;
        try {
            event.currentTarget.setPointerCapture(event.pointerId);
        } catch {
            // Pointer capture is best-effort; drag still works without it.
        }
        setDrag({ start: x, current: x });
        setHover(null);
    };

    const moveDrag = (event) => {
        if (!drag) {
            nearestHover(event);
            return;
        }
        const x = clampPlotX(pointerToViewX(event));
        setDrag((current) => ({ ...current, current: x }));
    };

    const finishDrag = (event) => {
        if (!drag) return;
        const current = clampPlotX(pointerToViewX(event));
        const start = Math.min(drag.start, current);
        const end = Math.max(drag.start, current);
        setDrag(null);
        if (end - start > 8) {
            setXDomain([xInvert(start), xInvert(end)]);
        }
        try {
            event.currentTarget.releasePointerCapture?.(event.pointerId);
        } catch {
            // No capture was held for this pointer.
        }
    };

    const zoom = (event) => {
        event.preventDefault();
        const px = clampPlotX(pointerToViewX(event));
        const center = xInvert(px);
        const factor = event.deltaY > 0 ? 1.22 : 0.82;
        const span = (domains.x[1] - domains.x[0]) * factor;
        let next = [center - span / 2, center + span / 2];
        next = [Math.max(domains.baseX[0], next[0]), Math.min(domains.baseX[1], next[1])];
        if (next[1] - next[0] > 1e-9) setXDomain(next);
    };

    const saveFigure = async (format) => {
        if (!svgRef.current) return;
        try {
            await saveSvgFigure(svgRef.current, effectivePlot?.title || file.name, format);
        } catch (saveError) {
            setError(saveError.message || 'Could not save the figure');
        }
    };

    if (error && !effectivePlot) return <div className="interactive-plot-error">{error}</div>;
    if (!effectivePlot) return <div className="interactive-plot-loading">Loading plot...</div>;

    // Keep the tooltip on the emptier side of the crosshair.
    const hoverOnLeftHalf = hover && hover.px < view.width / 2;

    return (
        <div className={`interactive-plot${wide ? ' interactive-plot--wide' : ''}`}>
            <div className="plot-toolbar">
                <div className="plot-legend">
                    {orderedSeries.map((series) => (
                        <button
                            key={series.label}
                            type="button"
                            className={hidden.has(series.label) ? 'muted' : ''}
                            onClick={() => {
                                setHidden((current) => {
                                    const next = new Set(current);
                                    if (next.has(series.label)) next.delete(series.label);
                                    else next.add(series.label);
                                    return next;
                                });
                            }}
                        >
                            <span
                                className={series.marker ? 'swatch-hollow' : ''}
                                style={series.marker ? { borderColor: series.color } : { background: series.color }}
                            />
                            {series.label}
                        </button>
                    ))}
                </div>
                <div className="plot-actions">
                    {xDomain && (
                        <button type="button" className="plot-reset" onClick={() => setXDomain(null)}>
                            Reset zoom
                        </button>
                    )}
                    <SaveMenu onSave={saveFigure} options={CHART_SAVE_OPTIONS} label="Save" align="right" />
                </div>
            </div>
            <div className="plot-stage">
                <svg
                    ref={svgRef}
                    viewBox={`0 0 ${view.width} ${view.height}`}
                    role="img"
                    aria-label={effectivePlot.title}
                    onPointerDown={startDrag}
                    onPointerMove={moveDrag}
                    onPointerUp={finishDrag}
                    onPointerCancel={() => setDrag(null)}
                    onPointerLeave={() => setHover(null)}
                    onWheel={zoom}
                    onDoubleClick={() => setXDomain(null)}
                >
                    <rect className="plot-bg" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                    {yTicks.ticks.map((tick) => (
                        <g key={`y-${tick}`}>
                            <line className="plot-grid-line" x1={view.left} x2={view.width - view.right} y1={yScale(tick)} y2={yScale(tick)} />
                            <text className="plot-tick" x={view.left - 10} y={yScale(tick) + 4.5} textAnchor="end">
                                {formatTick(tick, yTicks.step, yAxisMax)}
                            </text>
                        </g>
                    ))}
                    {xTicks.ticks.map((tick) => (
                        <g key={`x-${tick}`}>
                            <line
                                className="plot-tick-mark"
                                x1={xScale(tick)}
                                x2={xScale(tick)}
                                y1={view.top + plotHeight}
                                y2={view.top + plotHeight + 5}
                            />
                            <text className="plot-tick" x={xScale(tick)} y={view.height - 36} textAnchor="middle">
                                {formatTick(tick, xTicks.step, xAxisMax)}
                            </text>
                        </g>
                    ))}
                    {seriesShapes.map((series) => (
                        <path
                            key={series.label}
                            className={series.marker ? 'series-markers' : 'series-path'}
                            d={series.d}
                            stroke={series.color}
                        />
                    ))}
                    <rect className="plot-frame" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                    <AxisLabel label={effectivePlot.xLabel} x={view.left + plotWidth / 2} y={view.height - 10} />
                    <AxisLabel label={effectivePlot.yLabel} x={18} y={view.top + plotHeight / 2} rotate />
                    {hover && (
                        <g>
                            <line className="hover-line" x1={hover.px} x2={hover.px} y1={view.top} y2={view.top + plotHeight} />
                            {hover.values.map((value) => (
                                <circle key={value.label} className="hover-dot" cx={value.cx} cy={value.cy} r="3.6" fill={value.color} />
                            ))}
                        </g>
                    )}
                    {drag && (
                        <rect
                            className="zoom-selection"
                            x={Math.min(drag.start, drag.current)}
                            y={view.top}
                            width={Math.abs(drag.current - drag.start)}
                            height={plotHeight}
                        />
                    )}
                </svg>
                {hover && (
                    <div
                        className="plot-tooltip"
                        style={hoverOnLeftHalf
                            ? { left: `calc(${(hover.px / view.width) * 100}% + 14px)` }
                            : { right: `calc(${100 - (hover.px / view.width) * 100}% + 14px)` }}
                    >
                        <strong>{labelToText(effectivePlot.xLabel)}: {formatNumber(hover.x)}</strong>
                        {hover.values.map((value) => (
                            <span key={value.label}>
                                <i style={{ background: value.color }} />
                                {value.label}
                                <em>{formatNumber(value.y)}</em>
                            </span>
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
};

export default InteractivePlot;
