import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import './InteractivePlot.css';

const palette = ['#2f74d0', '#de8f35', '#3aa66f', '#c94f70', '#7b61d1', '#87923b', '#40a4b8', '#bc6e31'];

const formatNumber = (value) => {
    const abs = Math.abs(value);
    if (abs >= 1000 || (abs > 0 && abs < 0.01)) return value.toExponential(2);
    return value.toPrecision(4);
};

const formatInteger = (value) => String(Math.round(value));

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

const InteractivePlot = ({ file }) => {
    const [plot, setPlot] = useState(null);
    const [error, setError] = useState(null);
    const [hidden, setHidden] = useState(() => new Set());
    const [xDomain, setXDomain] = useState(null);
    const [hover, setHover] = useState(null);
    const [drag, setDrag] = useState(null);
    const svgRef = useRef(null);

    useEffect(() => {
        const fetchData = async () => {
            setPlot(null);
            setError(null);
            setHidden(new Set());
            setXDomain(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/plot/data`, {
                    params: { path: file.path }
                });
                setPlot(response.data);
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to load interactive plot');
            }
        };

        fetchData();
    }, [file.path]);

    const visibleSeries = useMemo(() => {
        return (plot?.series || []).filter((series) => !hidden.has(series.label));
    }, [plot, hidden]);

    const domains = useMemo(() => {
        const allX = visibleSeries.flatMap((series) => series.x);
        const baseX = niceDomain(allX);
        const currentX = xDomain || baseX;
        const allY = visibleSeries.flatMap((series) =>
            series.y.filter((_, index) => series.x[index] >= currentX[0] && series.x[index] <= currentX[1])
        );
        return { x: currentX, y: niceDomain(allY.length ? allY : visibleSeries.flatMap((series) => series.y)), baseX };
    }, [visibleSeries, xDomain]);

    const view = { width: 720, height: 405, left: 72, right: 14, top: 14, bottom: 62 };
    const plotWidth = view.width - view.left - view.right;
    const plotHeight = view.height - view.top - view.bottom;

    const xScale = (x) => view.left + ((x - domains.x[0]) / (domains.x[1] - domains.x[0] || 1)) * plotWidth;
    const yScale = (y) => view.top + plotHeight - ((y - domains.y[0]) / (domains.y[1] - domains.y[0] || 1)) * plotHeight;
    const xInvert = (px) => domains.x[0] + ((px - view.left) / plotWidth) * (domains.x[1] - domains.x[0]);

    const ticks = (domain, count = 4) => {
        const step = (domain[1] - domain[0]) / count;
        return Array.from({ length: count + 1 }, (_, index) => domain[0] + step * index);
    };

    const integerTicks = (domain, count = 4) => {
        const start = Math.ceil(domain[0]);
        const end = Math.floor(domain[1]);
        if (end <= start) return [Math.round((domain[0] + domain[1]) / 2)];
        const step = Math.max(1, Math.ceil((end - start) / count));
        const values = [];
        for (let value = start; value <= end; value += step) values.push(value);
        if (values[values.length - 1] !== end) values.push(end);
        return values.slice(0, count + 2);
    };

    const seriesPaths = visibleSeries.map((series, seriesIndex) => {
        const commands = [];
        series.x.forEach((x, index) => {
            if (x < domains.x[0] || x > domains.x[1]) return;
            const y = series.y[index];
            if (!Number.isFinite(x) || !Number.isFinite(y)) return;
            commands.push(`${commands.length ? 'L' : 'M'} ${xScale(x).toFixed(2)} ${yScale(y).toFixed(2)}`);
        });
        return { label: series.label, color: palette[seriesIndex % palette.length], d: commands.join(' ') };
    });

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
        if (!plot || !visibleSeries.length) return;
        const x = pointerToViewX(event);
        if (x < view.left || x > view.width - view.right) {
            setHover(null);
            return;
        }
        const dataX = xInvert(x);
        const values = visibleSeries.map((series, index) => {
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
                color: palette[index % palette.length],
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
        event.currentTarget.setPointerCapture(event.pointerId);
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
        if (end - start > 8) {
            setXDomain([xInvert(start), xInvert(end)]);
        }
        event.currentTarget.releasePointerCapture?.(event.pointerId);
        setDrag(null);
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

    if (error) return <div className="interactive-plot-error">{error}</div>;
    if (!plot) return <div className="interactive-plot-loading">Loading plot...</div>;

    return (
        <div className="interactive-plot">
            <div className="plot-actions">
                <button type="button" onClick={() => setXDomain(null)}>Reset</button>
            </div>
            <div className="plot-legend">
                {(plot.series || []).map((series, index) => (
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
                        <span style={{ background: palette[index % palette.length] }} />
                        {series.label}
                    </button>
                ))}
            </div>
            <svg
                ref={svgRef}
                viewBox={`0 0 ${view.width} ${view.height}`}
                role="img"
                aria-label={plot.title}
                onPointerDown={startDrag}
                onPointerMove={moveDrag}
                onPointerUp={finishDrag}
                onPointerCancel={() => setDrag(null)}
                onPointerLeave={() => setHover(null)}
                onWheel={zoom}
                onDoubleClick={() => setXDomain(null)}
            >
                <rect className="plot-bg" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                {ticks(domains.y).map((tick) => (
                    <g key={`y-${tick}`}>
                        <line className="plot-grid-line" x1={view.left} x2={view.width - view.right} y1={yScale(tick)} y2={yScale(tick)} />
                        <text className="plot-tick" x={view.left - 12} y={yScale(tick) + 4} textAnchor="end">{formatNumber(tick)}</text>
                    </g>
                ))}
                {integerTicks(domains.x).map((tick) => (
                    <g key={`x-${tick}`}>
                        <line className="plot-grid-line soft" x1={xScale(tick)} x2={xScale(tick)} y1={view.top} y2={view.top + plotHeight} />
                        <text className="plot-tick" x={xScale(tick)} y={view.height - 30} textAnchor="middle">{formatInteger(tick)}</text>
                    </g>
                ))}
                {seriesPaths.map((series) => (
                    <path key={series.label} className="series-path" d={series.d} stroke={series.color} />
                ))}
                <rect className="plot-frame" x={view.left} y={view.top} width={plotWidth} height={plotHeight} />
                <AxisLabel label={plot.xLabel} x={view.left + plotWidth / 2} y={view.height - 7} />
                <AxisLabel label={plot.yLabel} x={18} y={view.top + plotHeight / 2} rotate />
                {hover && (
                    <g>
                        <line className="hover-line" x1={hover.px} x2={hover.px} y1={view.top} y2={view.top + plotHeight} />
                        {hover.values.map((value) => (
                            <circle key={value.label} cx={value.cx} cy={value.cy} r="3.2" fill={value.color} />
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
                <div className="plot-tooltip">
                    <strong>{plot.xLabel}: {formatNumber(hover.x)}</strong>
                    {hover.values.map((value) => (
                        <span key={value.label}>{value.label}: {formatNumber(value.y)}</span>
                    ))}
                </div>
            )}
        </div>
    );
};

export default InteractivePlot;
