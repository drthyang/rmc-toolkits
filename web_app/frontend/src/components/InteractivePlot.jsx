// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useId, useMemo, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import axios from 'axios';
import API_BASE_URL from '../api';
import { saveSvgFigure } from '../figureExport';
import { nearestFiniteIndex, niceDomain, plotPayloadError } from '../plotDomain';
import { GUIDE_STROKE, PLOT_PALETTE } from '../plotPalette';
import { Banner, Pill } from '../ui';
import { uiScale } from '../uiScale';
import SaveMenu from '../ui/SaveMenu';
import './InteractivePlot.css';

const palette = PLOT_PALETTE;

const CHART_SAVE_OPTIONS = [
    { id: 'png', label: 'PNG image', hint: '.png' },
    { id: 'svg', label: 'SVG vector', hint: '.svg' },
];

// Series whose label mentions "exp" hold measured data (Experiment, F(Q)_Expt, X_ray_exp_renorm, ...).
const isExperimental = (label) => /exp/i.test(label);

// Theory/reference lines (series with role: 'guide') are drawn dashed and
// muted, and are excluded from the palette rotation and from hover snapping. A
// caller may set an explicit `color` to tie one to the series it belongs to.
//
// Opt-in axis and series options (each does nothing when absent, so a payload
// without them renders exactly as before — pinned by interactivePlotAxes.test):
//   plotData.xDomain    [min, max]: the un-zoomed x range, unpadded; also
//                       bounds the wheel zoom.
//   plotData.xTicks     labelled major ticks while un-zoomed (nice ticks after
//                       a zoom).
//   plotData.xMinorStep unlabelled minor tick marks at this step, un-zoomed.
//   plotData.xGrid      vertical grid lines at the major ticks.
//   plotData.yMin       the un-zoomed y range starts here (e.g. 0).
//   series.curve        'step': a histogram outline, flat across each x ± half
//                       a bin (series.binWidth, else the x spacing).
//   series.fill         a light area under the curve, down to y = 0.
//   series.width        stroke width in px.
//   series.legend       false: left out of the legend (e.g. window guides).
//
// Opt-in toolbar props (absent, the toolbar row renders as before):
//   legend={false}      no series legend — for a card whose title already
//                       names the curves; the actions stay on the right.
//   actionsTarget       an element the caller owns (a card header slot): the
//                       actions (Reset zoom, Save) render there, and with
//                       legend={false} the plot has no toolbar row at all, so
//                       a short card gives that height to the plot. While the
//                       prop is present but null (the slot not mounted yet)
//                       the actions render nowhere, never inline.

const MARKER_RADIUS = 2.8;

const formatNumber = (value) => {
    const abs = Math.abs(value);
    if (abs >= 1e6 || (abs > 0 && abs < 0.01)) return value.toExponential(2);
    if (abs >= 1000) return Math.round(value).toLocaleString();
    return value.toPrecision(4);
};

const SUPERSCRIPTS = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };

// Render `Q (Å^{-1})` style labels as plain text with unicode superscripts.
const labelToText = (label = '') => String(label ?? '').replace(/\^\{([^}]+)\}/g, (_, exponent) =>
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

const finitePair = (pair) => Array.isArray(pair) && pair.length === 2
    && pair.every(Number.isFinite) && pair[1] > pair[0];

// Half the bin width of a step series: its own binWidth, else the smallest
// positive spacing of its x values (bin centres of a uniform histogram).
const stepHalfWidth = (series) => {
    if (Number.isFinite(series.binWidth) && series.binWidth > 0) return series.binWidth / 2;
    let spacing = Infinity;
    for (let index = 1; index < series.x.length; index += 1) {
        const gap = series.x[index] - series.x[index - 1];
        if (gap > 0 && gap < spacing) spacing = gap;
    }
    return Number.isFinite(spacing) ? spacing / 2 : 0.5;
};

const AxisLabel = ({ label: rawLabel, x, y, textAnchor = 'middle', rotate = false }) => {
    const label = String(rawLabel ?? '');
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

const InteractivePlot = ({ file, variant, plotData, refreshKey, legend = true, actionsTarget }) => {
    const wide = variant === 'wide';
    // 'fit' takes its viewBox from the rendered box instead of a fixed aspect,
    // so the drawing fills the card rather than letterboxing inside it. Opt-in:
    // it only suits callers that give the plot a definite height.
    const fit = variant === 'fit';
    const stageRef = useRef(null);
    const [stageSize, setStageSize] = useState(null);
    useEffect(() => {
        if (!fit) return undefined;
        const node = stageRef.current;
        if (!node) return undefined;
        const observer = new ResizeObserver(([entry]) => {
            const { width, height } = entry.contentRect;
            // Round to whole pixels: sub-pixel churn would rebuild the viewBox
            // (and every tick label) on any reflow. The UI scale (2K / 4K root
            // type size) is read with the size: a window resize changes both.
            setStageSize((current) => {
                const next = { width: Math.round(width), height: Math.round(height), scale: uiScale() };
                return current && current.width === next.width && current.height === next.height
                    && current.scale === next.scale
                    ? current
                    : next;
            });
        });
        observer.observe(node);
        return () => observer.disconnect();
    }, [fit]);
    const [plot, setPlot] = useState(null);
    const [error, setError] = useState(null);
    // A failed figure save: shown under the toolbar (the chart stays), cleared
    // by the next save. Separate from `error`, which replaces an unloaded chart.
    const [saveError, setSaveError] = useState(null);
    const [hidden, setHidden] = useState(() => new Set());
    const [xDomain, setXDomain] = useState(null);
    const [yDomain, setYDomain] = useState(null);
    const [hover, setHover] = useState(null);
    const [drag, setDrag] = useState(null);
    // Unique clip-path id so a rectangle-zoomed series is clipped to the plot area
    // (and does not draw over the axes) even with several charts on the page.
    const clipId = `plot-clip-${useId().replace(/[^a-zA-Z0-9_-]/g, '')}`;
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
                // A body that is not valid JSON arrives as a raw string; never
                // hand that to the renderer — say what went wrong instead.
                const payloadError = plotPayloadError(response.data);
                if (payloadError) {
                    setPlot(null);
                    setError(payloadError);
                    return;
                }
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

    // Reset view state (zoom, hover, hidden series) when direct plotData is
    // swapped for a different dataset, honoring per-series defaultHidden flags.
    // Render-time adjustment (not an effect) per the React "adjusting state
    // when a prop changes" pattern.
    const [lastPlotData, setLastPlotData] = useState(null);
    // With an initialYDomain, double-click toggles out to the full data
    // extent (and back); without one it just resets zoom as before.
    const [fullExtent, setFullExtent] = useState(false);
    if (plotData && plotData !== lastPlotData) {
        setLastPlotData(plotData);
        setHidden(new Set((plotData.series || [])
            .filter((entry) => entry.defaultHidden)
            .map((entry) => entry.label)));
        setXDomain(null);
        setYDomain(null);
        setHover(null);
        setDrag(null);
        setFullExtent(false);
    }

    // Measured data first (drawn underneath, first palette color); the
    // calculated curve follows and is drawn on top of the hollow markers.
    // Guide series draw last (on top, dashed) and never consume palette slots.
    const orderedSeries = useMemo(() => {
        const series = effectivePlot?.series || [];
        const guides = series.filter((entry) => entry.role === 'guide');
        const data = series.filter((entry) => entry.role !== 'guide');
        const experimental = data.filter((entry) => isExperimental(entry.label));
        const calculated = data.filter((entry) => !isExperimental(entry.label));
        const paired = experimental.length > 0 && calculated.length > 0;
        const reordered = paired ? [...experimental, ...calculated] : data;
        return [
            ...reordered.map((entry, index) => ({
                ...entry,
                marker: paired && isExperimental(entry.label),
                color: entry.color || palette[index % palette.length]
            })),
            ...guides.map((entry) => ({
                ...entry,
                marker: false,
                guide: true,
                color: entry.color || GUIDE_STROKE
            }))
        ];
    }, [effectivePlot]);

    const visibleSeries = useMemo(() => {
        return orderedSeries.filter((series) => !hidden.has(series.label));
    }, [orderedSeries, hidden]);

    const domains = useMemo(() => {
        const fixedX = finitePair(effectivePlot?.xDomain) ? effectivePlot.xDomain : null;
        const baseX = fixedX || niceDomain(visibleSeries.flatMap((series) => series.x));
        const currentX = xDomain || baseX;
        const allY = visibleSeries.flatMap((series) =>
            series.y.filter((_, index) => series.x[index] >= currentX[0] && series.x[index] <= currentX[1])
        );
        let baseY = niceDomain(allY.length ? allY : visibleSeries.flatMap((series) => series.y));
        const yMin = effectivePlot?.yMin;
        if (Number.isFinite(yMin)) baseY = [yMin, Math.max(baseY[1], yMin + 1e-9)];
        // A caller-supplied initial y window (e.g. the Auto StoG low-r zoom)
        // acts as the un-zoomed default; user zooms override it, and
        // double-click toggles the full data extent.
        const initialY = fullExtent ? null : effectivePlot?.initialYDomain;
        return { x: currentX, y: yDomain || initialY || baseY, baseX, baseY };
    }, [visibleSeries, xDomain, yDomain, effectivePlot, fullExtent]);

    // Under 'fit' the box is whatever the card gives, so tick density follows
    // it — roughly one y tick per 70px and one x tick per 95px, clamped to the
    // range the fixed variants use. A short card otherwise crams in six y ticks.
    // The fitted box is in user units: one unit is uiScale() CSS pixels, so
    // margins, tick labels and tick density grow with the rest of the UI on
    // 2K / 4K monitors (one unit is one pixel at the 15px design size).
    const fitted = fit && stageSize?.width && stageSize?.height
        ? { width: stageSize.width / stageSize.scale, height: stageSize.height / stageSize.scale }
        : null;
    const clampTicks = (value, low, high) => Math.max(low, Math.min(high, Math.round(value)));
    const yTicks = niceTicks(
        domains.y,
        fitted ? clampTicks(fitted.height / 70, 3, 8) : wide ? 4 : 6
    );
    // Caller-fixed major ticks hold only while un-zoomed; a zoom hands the
    // axis back to nice ticks.
    const fixedTicks = !xDomain && Array.isArray(effectivePlot?.xTicks)
        ? effectivePlot.xTicks.filter((tick) => Number.isFinite(tick) && tick >= domains.x[0] && tick <= domains.x[1])
        : null;
    const xTicks = fixedTicks?.length
        ? { ticks: fixedTicks, step: fixedTicks.length > 1 ? fixedTicks[1] - fixedTicks[0] : 1 }
        : niceTicks(
            domains.x,
            fitted ? clampTicks(fitted.width / 95, 4, 12) : wide ? 11 : 7
        );
    const minorStep = effectivePlot?.xMinorStep;
    const xMinorTicks = [];
    if (!xDomain && Number.isFinite(minorStep) && minorStep > 0) {
        const span = domains.x[1] - domains.x[0];
        if (span / minorStep <= 400) {
            const first = Math.ceil(domains.x[0] / minorStep - 1e-9);
            const last = Math.floor(domains.x[1] / minorStep + 1e-9);
            for (let k = first; k <= last; k += 1) {
                const tick = k * minorStep;
                if (!xTicks.ticks.some((major) => Math.abs(major - tick) < minorStep * 1e-6)) xMinorTicks.push(tick);
            }
        }
    }
    const yAxisMax = Math.max(...yTicks.ticks.map(Math.abs), 0);
    const xAxisMax = Math.max(...xTicks.ticks.map(Math.abs), 0);

    // The default left margin fits ~3-character y ticks beside the rotated
    // axis label (the historical look). Wider ticks — e.g. counts like 7000 —
    // push the plot area right by one 14.5px tabular digit per extra
    // character, so they never run into the label.
    const widestYTick = Math.max(
        ...yTicks.ticks.map((tick) => formatTick(tick, yTicks.step, yAxisMax).length),
        1
    );
    // 8.7px per extra 14.5px tabular digit, plus a fixed 6px of breathing
    // room once padding engages (the default margin fits 3 characters with
    // no slack left over).
    const leftPad = widestYTick > 3 ? (widestYTick - 3) * 8.7 + 6 : 0;

    // 8:5 (golden-ish) for grid cards, a slim strip for the wide variant, the
    // measured box for 'fit'. Under 'fit' one user unit is one CSS pixel at
    // the 15px design size (uiScale() pixels above it), so the margins below
    // keep the size they have relative to the UI; before the first
    // measurement it falls back to the 8:5 box.
    const view = fitted
        ? {
            width: fitted.width,
            height: fitted.height,
            left: 60 + leftPad,
            right: 18,
            top: 16,
            bottom: 58
        }
        : wide
            ? { width: 1440, height: 320, left: 64 + leftPad, right: 20, top: 18, bottom: 58 }
            : { width: 720, height: 450, left: 60 + leftPad, right: 18, top: 16, bottom: 58 };
    const plotWidth = view.width - view.left - view.right;
    const plotHeight = view.height - view.top - view.bottom;

    const xScale = (x) => view.left + ((x - domains.x[0]) / (domains.x[1] - domains.x[0] || 1)) * plotWidth;
    const yScale = (y) => view.top + plotHeight - ((y - domains.y[0]) / (domains.y[1] - domains.y[0] || 1)) * plotHeight;
    const xInvert = (px) => domains.x[0] + ((px - view.left) / plotWidth) * (domains.x[1] - domains.x[0]);
    const yInvert = (py) => domains.y[0] + ((view.top + plotHeight - py) / plotHeight) * (domains.y[1] - domains.y[0]);

    // Geometry is memoized so hover re-renders skip rebuilding the (large)
    // marker paths; only domain or series changes recompute them.
    const seriesShapes = useMemo(() => {
        const sx = (x) => view.left + ((x - domains.x[0]) / (domains.x[1] - domains.x[0] || 1)) * plotWidth;
        const sy = (y) => view.top + plotHeight - ((y - domains.y[0]) / (domains.y[1] - domains.y[0] || 1)) * plotHeight;
        const fmt = (value) => value.toFixed(2);
        return visibleSeries.map((series) => {
            const extra = {
                width: Number.isFinite(series.width) ? series.width : undefined,
                step: series.curve === 'step'
            };
            if (series.curve === 'step') {
                // Flat across each bin, vertical at the shared edges; a gap
                // (non-finite y) or a missing bin breaks the outline.
                const half = stepHalfWidth(series);
                const runs = [];
                let run = null;
                series.x.forEach((x, index) => {
                    const y = series.y[index];
                    if (!Number.isFinite(x) || !Number.isFinite(y)
                        || x + half < domains.x[0] || x - half > domains.x[1]) {
                        run = null;
                        return;
                    }
                    const left = x - half;
                    if (!run || Math.abs(left - run.right) > half * 1e-6) {
                        run = { right: x + half, points: [] };
                        runs.push(run);
                    }
                    run.points.push([sx(left), sx(x + half), sy(y)]);
                    run.right = x + half;
                });
                const d = runs.map(({ points }) => points.map(([px0, px1, py], index) => (
                    `${index ? 'L' : 'M'} ${fmt(px0)} ${fmt(py)} L ${fmt(px1)} ${fmt(py)}`
                )).join(' ')).join(' ');
                const base = sy(0);
                const area = series.fill
                    ? runs.map(({ points }) => {
                        const first = points[0];
                        const last = points[points.length - 1];
                        const top = points.map(([px0, px1, py]) => `L ${fmt(px0)} ${fmt(py)} L ${fmt(px1)} ${fmt(py)}`).join(' ');
                        return `M ${fmt(first[0])} ${fmt(base)} ${top} L ${fmt(last[1])} ${fmt(base)} Z`;
                    }).join(' ')
                    : null;
                return { label: series.label, color: series.color, marker: false, guide: Boolean(series.guide), d, area, ...extra };
            }
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
                return { label: series.label, color: series.color, marker: true, guide: false, d: commands.join(' ') };
            }
            const d = points.map(([px, py], index) => `${index ? 'L' : 'M'} ${px.toFixed(2)} ${py.toFixed(2)}`).join(' ');
            const base = sy(0);
            const area = series.fill && points.length
                ? `M ${fmt(points[0][0])} ${fmt(base)} ${points.map(([px, py]) => `L ${fmt(px)} ${fmt(py)}`).join(' ')} L ${fmt(points[points.length - 1][0])} ${fmt(base)} Z`
                : null;
            return { label: series.label, color: series.color, marker: false, guide: Boolean(series.guide), d, area, ...extra };
        });
    }, [visibleSeries, domains, view.left, view.top, plotWidth, plotHeight]);

    const pointerToView = (event) => {
        const svg = svgRef.current;
        if (!svg) return { x: view.left, y: view.top };
        const transform = svg.getScreenCTM();
        if (!transform) return { x: view.left, y: view.top };
        const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(transform.inverse());
        return { x: point.x, y: point.y };
    };
    const pointerToViewX = (event) => pointerToView(event).x;

    const clampPlotX = (x) => Math.max(view.left, Math.min(view.width - view.right, x));
    const clampPlotY = (y) => Math.max(view.top, Math.min(view.height - view.bottom, y));

    const nearestHover = (event) => {
        const hoverSeries = visibleSeries.filter((series) => !series.guide);
        if (!effectivePlot || !hoverSeries.length) return;
        const x = pointerToViewX(event);
        if (x < view.left || x > view.width - view.right) {
            setHover(null);
            return;
        }
        const dataX = xInvert(x);
        // Only drawable (finite x AND y) points can be snapped to: a masked
        // region arrives as null/NaN gaps, which must not win the search.
        const values = hoverSeries.map((series) => {
            const best = nearestFiniteIndex(series.x, series.y, dataX);
            if (best < 0) return null;
            return {
                label: series.label,
                color: series.color,
                x: series.x[best],
                y: series.y[best],
                cx: xScale(series.x[best]),
                cy: yScale(series.y[best])
            };
        }).filter(Boolean);
        if (!values.length) {
            setHover(null);
            return;
        }
        setHover({ x: values[0].x, px: xScale(values[0].x), values });
    };

    const startDrag = (event) => {
        const { x, y } = pointerToView(event);
        if (x < view.left || x > view.width - view.right) return;
        try {
            event.currentTarget.setPointerCapture(event.pointerId);
        } catch {
            // Pointer capture is best-effort; drag still works without it.
        }
        setDrag({ x0: clampPlotX(x), y0: clampPlotY(y), x1: clampPlotX(x), y1: clampPlotY(y) });
        setHover(null);
    };

    const moveDrag = (event) => {
        if (!drag) {
            nearestHover(event);
            return;
        }
        const { x, y } = pointerToView(event);
        setDrag((current) => ({ ...current, x1: clampPlotX(x), y1: clampPlotY(y) }));
    };

    const finishDrag = (event) => {
        if (!drag) return;
        const { x, y } = pointerToView(event);
        const xLo = Math.min(drag.x0, clampPlotX(x));
        const xHi = Math.max(drag.x0, clampPlotX(x));
        const yLo = Math.min(drag.y0, clampPlotY(y));
        const yHi = Math.max(drag.y0, clampPlotY(y));
        setDrag(null);
        // Zoom whichever dimension(s) the drag actually spans (> 8px): a real box
        // zooms both axes to it, a thin horizontal/vertical drag zooms just that
        // axis. Screen y grows downward, so the top pixel is the higher value.
        const zoomX = xHi - xLo > 8;
        const zoomY = yHi - yLo > 8;
        if (zoomX) setXDomain([xInvert(xLo), xInvert(xHi)]);
        if (zoomY) setYDomain([yInvert(yHi), yInvert(yLo)]);
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
        setSaveError(null);
        try {
            await saveSvgFigure(svgRef.current, effectivePlot?.title || file.name, format);
        } catch (failure) {
            setSaveError(failure?.message || 'Could not save the figure');
        }
    };

    if (error && !effectivePlot) return <div className="ui-loading ui-loading--error">{error}</div>;
    if (!effectivePlot) return <div className="ui-loading">Loading plot…</div>;

    // Keep the tooltip on the emptier side of the crosshair.
    const hoverOnLeftHalf = hover && hover.px < view.width / 2;

    const actions = (
        <div className="plot-actions">
            {(xDomain || yDomain || fullExtent) && (
                <Pill
                    tint
                    onClick={() => { setXDomain(null); setYDomain(null); setFullExtent(false); }}
                >
                    Reset zoom
                </Pill>
            )}
            <SaveMenu onSave={saveFigure} options={CHART_SAVE_OPTIONS} label="Save" align="right" />
        </div>
    );
    const actionsElsewhere = actionsTarget !== undefined;

    return (
        <div className={`interactive-plot${wide ? ' interactive-plot--wide' : ''}${fit ? ' interactive-plot--fit' : ''}`}>
            {actionsElsewhere && actionsTarget && createPortal(actions, actionsTarget)}
            {(legend || !actionsElsewhere) && (
                <div className="plot-toolbar">
                    <div className="plot-legend">
                        {legend && orderedSeries.filter((series) => series.legend !== false).map((series) => (
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
                                    className={series.marker ? 'swatch-hollow' : series.guide ? 'swatch-guide' : ''}
                                    style={series.marker
                                        ? { borderColor: series.color }
                                        : series.guide
                                            ? { borderColor: series.color }
                                            : { background: series.color }}
                                />
                                {series.label}
                            </button>
                        ))}
                    </div>
                    {!actionsElsewhere && actions}
                </div>
            )}
            {saveError && (
                <Banner tone="danger" sm role="alert" className="plot-save-error" onDismiss={() => setSaveError(null)}>
                    {saveError}
                </Banner>
            )}
            <div className="plot-stage" ref={stageRef}>
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
                    onDoubleClick={() => {
                        const atDefault = !xDomain && !yDomain;
                        setXDomain(null);
                        setYDomain(null);
                        if (effectivePlot?.initialYDomain) {
                            setFullExtent(atDefault ? !fullExtent : false);
                        }
                    }}
                >
                    <defs>
                        <clipPath id={clipId}>
                            <rect x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                        </clipPath>
                    </defs>
                    <rect className="plot-bg" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                    {effectivePlot.xGrid && xTicks.ticks.map((tick) => (
                        <line key={`xg-${tick}`} className="plot-grid-line" x1={xScale(tick)} x2={xScale(tick)} y1={view.top} y2={view.top + plotHeight} />
                    ))}
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
                    {xMinorTicks.map((tick) => (
                        <line
                            key={`xm-${tick}`}
                            className="plot-tick-mark plot-tick-mark--minor"
                            x1={xScale(tick)}
                            x2={xScale(tick)}
                            y1={view.top + plotHeight}
                            y2={view.top + plotHeight + 5}
                        />
                    ))}
                    <g clipPath={`url(#${clipId})`}>
                        {seriesShapes.filter((series) => series.area).map((series) => (
                            <path key={`area-${series.label}`} className="series-area" d={series.area} fill={series.color} />
                        ))}
                        {seriesShapes.map((series) => (
                            <path
                                key={series.label}
                                className={series.marker
                                    ? 'series-markers'
                                    : `${series.guide ? 'series-path series-path--guide' : 'series-path'}${series.step ? ' series-path--step' : ''}`}
                                d={series.d}
                                stroke={series.color}
                                style={series.width !== undefined ? { strokeWidth: series.width } : undefined}
                            />
                        ))}
                    </g>
                    <rect className="plot-frame" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                    <AxisLabel label={effectivePlot.xLabel} x={view.left + plotWidth / 2} y={view.height - 10} />
                    <AxisLabel label={effectivePlot.yLabel} x={18} y={view.top + plotHeight / 2} rotate />
                    {hover && (
                        <g clipPath={`url(#${clipId})`}>
                            <line className="hover-line" x1={hover.px} x2={hover.px} y1={view.top} y2={view.top + plotHeight} />
                            {hover.values.map((value) => (
                                <circle key={value.label} className="hover-dot" cx={value.cx} cy={value.cy} r="3.6" fill={value.color} />
                            ))}
                        </g>
                    )}
                    {drag && (
                        <rect
                            className="zoom-selection"
                            x={Math.min(drag.x0, drag.x1)}
                            y={Math.min(drag.y0, drag.y1)}
                            width={Math.abs(drag.x1 - drag.x0)}
                            height={Math.abs(drag.y1 - drag.y0)}
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
