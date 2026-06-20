import { computeDensityGpu, shouldUseGpu } from './gpuKde.js';

const CUBE_CORNERS = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]
];

const CUBE_EDGES = [
    [0, 1], [0, 2], [1, 3], [2, 3],
    [4, 5], [4, 6], [5, 7], [6, 7],
    [0, 4], [1, 5], [2, 6], [3, 7]
];

const dot = (a, b) => a.reduce((sum, value, index) => sum + value * b[index], 0);
const add = (a, b) => a.map((value, index) => value + b[index]);
const subtract = (a, b) => a.map((value, index) => value - b[index]);
const scale = (vector, factor) => vector.map((value) => value * factor);
const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
];

const vectorLength = (vector) => Math.sqrt(dot(vector, vector));
const normalize = (vector, fallback = [0, 0, 1]) => {
    const length = vectorLength(vector);
    if (length <= 1e-9) return fallback;
    return vector.map((value) => value / length);
};

const makeFreePlaneBasis = (normal) => {
    const reference = Math.abs(normal[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
    const u = normalize(subtract(reference, scale(normal, dot(reference, normal))), [0, 1, 0]);
    const v = normalize(cross(normal, u), [0, 0, 1]);
    return { u, v };
};

const planeSectionVertices = (normal, offset) => {
    const vertices = [];
    CUBE_EDGES.forEach(([start, end]) => {
        const p0 = CUBE_CORNERS[start];
        const p1 = CUBE_CORNERS[end];
        const d0 = dot(p0, normal) - offset;
        const d1 = dot(p1, normal) - offset;
        if (Math.abs(d0) <= 1e-9) vertices.push(p0);
        if (Math.abs(d1) <= 1e-9) vertices.push(p1);
        if (d0 * d1 < 0) {
            const t = d0 / (d0 - d1);
            vertices.push(add(p0, scale(subtract(p1, p0), t)));
        }
    });

    const unique = [];
    vertices.forEach((vertex) => {
        if (!unique.some((existing) => vectorLength(subtract(vertex, existing)) <= 1e-8)) {
            unique.push(vertex);
        }
    });
    if (unique.length < 3) return [];

    const { u, v } = makeFreePlaneBasis(normal);
    const center = scale(unique.reduce((acc, vertex) => add(acc, vertex), [0, 0, 0]), 1 / unique.length);
    return unique.sort((a, b) => {
        const aDelta = subtract(a, center);
        const bDelta = subtract(b, center);
        return Math.atan2(dot(aDelta, v), dot(aDelta, u)) - Math.atan2(dot(bDelta, v), dot(bDelta, u));
    });
};

const makeSlab = ({ points, normal, uVector, vVector, range, zCenter, thickness }) => {
    const depthSpan = range[1] - range[0] || 1;
    const slab = [];
    points.forEach((point) => {
        const fraction = [point.x, point.y, point.z];
        const normalizedDepth = (dot(fraction, normal) - range[0]) / depthSpan;
        if (Math.abs(normalizedDepth - zCenter) <= thickness / 2) {
            slab.push([dot(fraction, uVector), dot(fraction, vVector)]);
        }
    });
    return slab;
};

const randomUnit = (seed) => {
    let value = seed >>> 0;
    return () => {
        value += 0x6D2B79F5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
};

const sampleWithoutReplacement = (items, limit, seed = 0) => {
    if (items.length <= limit) return items;
    const random = randomUnit(seed);
    const indices = Array.from({ length: items.length }, (_, index) => index);
    for (let index = 0; index < limit; index += 1) {
        const swapIndex = index + Math.floor(random() * (indices.length - index));
        [indices[index], indices[swapIndex]] = [indices[swapIndex], indices[index]];
    }
    return indices.slice(0, limit).map((index) => items[index]);
};

const covariance = (samples) => {
    const n = samples.length;
    const mean = samples.reduce((acc, sample) => [acc[0] + sample[0], acc[1] + sample[1]], [0, 0])
        .map((value) => value / Math.max(n, 1));
    let c00 = 0;
    let c01 = 0;
    let c11 = 0;
    samples.forEach(([x, y]) => {
        const dx = x - mean[0];
        const dy = y - mean[1];
        c00 += dx * dx;
        c01 += dx * dy;
        c11 += dy * dy;
    });
    const denom = Math.max(n - 1, 1);
    return { c00: c00 / denom, c01: c01 / denom, c11: c11 / denom };
};

const makeKernel = (samples, bandwidth) => {
    const cov = covariance(samples);
    const factor = Math.max(Number(bandwidth) || 0.03, 1e-4);
    const scaleFactor = factor * factor;
    const regularizer = 1e-8;
    let k00 = cov.c00 * scaleFactor + regularizer;
    let k01 = cov.c01 * scaleFactor;
    let k11 = cov.c11 * scaleFactor + regularizer;
    let det = k00 * k11 - k01 * k01;

    if (det <= 1e-12) {
        const spread = Math.max(cov.c00, cov.c11, 1e-4) * scaleFactor + 1e-6;
        k00 += spread;
        k11 += spread;
        k01 = 0;
        det = k00 * k11;
    }

    const inv00 = k11 / det;
    const inv01 = -k01 / det;
    const inv11 = k00 / det;
    const normalizer = 1 / (2 * Math.PI * Math.sqrt(det) * samples.length);
    return { inv00, inv01, inv11, normalizer };
};

const extractContours = ({ density, grid, xMin, xMax, yMin, yMax, vmin, vmax, levels = 8 }) => {
    if (levels <= 0 || !(vmax > vmin)) return [];
    const xStep = (xMax - xMin) / Math.max(grid - 1, 1);
    const yStep = (yMax - yMin) / Math.max(grid - 1, 1);
    const contours = [];

    const interpolate = (level, a, b) => {
        const denom = b.value - a.value;
        const t = Math.abs(denom) <= 1e-12 ? 0.5 : (level - a.value) / denom;
        return [a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t];
    };

    for (let levelIndex = 1; levelIndex <= levels; levelIndex += 1) {
        const level = vmin + (levelIndex / (levels + 1)) * (vmax - vmin);
        const lines = [];
        for (let y = 0; y < grid - 1; y += 1) {
            for (let x = 0; x < grid - 1; x += 1) {
                const corners = [
                    { x: xMin + x * xStep, y: yMin + y * yStep, value: density[y][x] },
                    { x: xMin + (x + 1) * xStep, y: yMin + y * yStep, value: density[y][x + 1] },
                    { x: xMin + (x + 1) * xStep, y: yMin + (y + 1) * yStep, value: density[y + 1][x + 1] },
                    { x: xMin + x * xStep, y: yMin + (y + 1) * yStep, value: density[y + 1][x] }
                ];
                const edgePoints = [];
                [[0, 1], [1, 2], [2, 3], [3, 0]].forEach(([start, end]) => {
                    const a = corners[start];
                    const b = corners[end];
                    if ((a.value < level && b.value >= level) || (b.value < level && a.value >= level)) {
                        edgePoints.push(interpolate(level, a, b));
                    }
                });
                if (edgePoints.length === 2) {
                    lines.push(edgePoints);
                } else if (edgePoints.length === 4) {
                    lines.push([edgePoints[0], edgePoints[1]], [edgePoints[2], edgePoints[3]]);
                }
            }
        }
        if (lines.length) contours.push({ level, lines });
    }
    return contours;
};

// The density map is the worker's hot loop: for every grid cell, sum an
// anisotropic 2D Gaussian over every sample (O(grid^2 * samples)). This is the
// CPU implementation, used directly on devices without WebGPU and as the
// fallback whenever the GPU path is unavailable or errors.
const computeDensityCpu = ({ samples, kernel, grid, xMin, yMin, xStep, yStep }) => {
    const density = Array.from({ length: grid }, () => new Array(grid).fill(0));
    for (let y = 0; y < grid; y += 1) {
        const gy = yMin + y * yStep;
        for (let x = 0; x < grid; x += 1) {
            const gx = xMin + x * xStep;
            let sum = 0;
            for (let index = 0; index < samples.length; index += 1) {
                const dx = gx - samples[index][0];
                const dy = gy - samples[index][1];
                const exponent = -0.5 * (kernel.inv00 * dx * dx + 2 * kernel.inv01 * dx * dy + kernel.inv11 * dy * dy);
                if (exponent > -60) sum += Math.exp(exponent);
            }
            density[y][x] = sum * kernel.normalizer;
        }
    }
    return density;
};

const computeKde = async (payload) => {
    const {
        points,
        normal,
        uVector,
        vVector,
        range,
        zCenter,
        thickness,
        bandwidth,
        gridSize,
        logScale
    } = payload;
    const projectedCorners = CUBE_CORNERS.map((corner) => [dot(corner, uVector), dot(corner, vVector)]);
    const xValues = projectedCorners.map(([x]) => x);
    const yValues = projectedCorners.map(([, y]) => y);
    const xMin = Math.min(...xValues);
    const xMax = Math.max(...xValues);
    const yMin = Math.min(...yValues);
    const yMax = Math.max(...yValues);
    const slab = makeSlab({ points, normal, uVector, vVector, range, zCenter, thickness });
    const grid = Math.max(16, Math.min(Number(gridSize) || 120, 260));
    let density = null;
    let fitCount = 0;
    let backend = 'cpu';

    if (slab.length >= 5) {
        const fitLimit = 6000;
        const samples = sampleWithoutReplacement(slab, fitLimit, 0);
        fitCount = samples.length;

        if (samples.length >= 5) {
            const kernel = makeKernel(samples, bandwidth);
            const xStep = (xMax - xMin) / Math.max(grid - 1, 1);
            const yStep = (yMax - yMin) / Math.max(grid - 1, 1);
            const args = { samples, kernel, grid, xMin, yMin, xStep, yStep };

            // Run the density map on the GPU when the workload is large enough to
            // amortize the setup cost. Any failure or unavailability falls back to
            // the CPU loop so the output is identical on every device.
            let mapped = null;
            if (shouldUseGpu(grid, samples.length)) {
                try {
                    mapped = await computeDensityGpu(args);
                } catch {
                    mapped = null;
                }
            }
            density = mapped ?? computeDensityCpu(args);
            if (mapped) backend = 'gpu';
        }
    }

    // Degenerate slabs never reach the density solver; emit a flat grid so the
    // scaling and contour passes below always have a well-formed array.
    if (!density) {
        density = Array.from({ length: grid }, () => new Array(grid).fill(0));
    }

    let vmin = Infinity;
    let vmax = -Infinity;
    for (let y = 0; y < grid; y += 1) {
        for (let x = 0; x < grid; x += 1) {
            if (logScale) density[y][x] = Math.log10(density[y][x] + 1e-12);
            vmin = Math.min(vmin, density[y][x]);
            vmax = Math.max(vmax, density[y][x]);
        }
    }

    const depthSpan = range[1] - range[0] || 1;
    const centerDepth = range[0] + zCenter * depthSpan;
    const planeVertices = planeSectionVertices(normal, centerDepth);
    return {
        density,
        extent: [xMin, xMax, yMin, yMax],
        grid,
        z: zCenter,
        dz: thickness,
        bw: bandwidth,
        log: logScale,
        slabCount: slab.length,
        fitCount,
        vmin: Number.isFinite(vmin) ? vmin : 0,
        vmax: Number.isFinite(vmax) ? vmax : 0,
        contours: extractContours({ density, grid, xMin, xMax, yMin, yMax, vmin, vmax }),
        center: zCenter,
        thickness,
        normal,
        uVector,
        vVector,
        planeVertices,
        planePolygon: planeVertices.map((vertex) => [dot(vertex, uVector), dot(vertex, vVector)]),
        browserKde: true,
        backend
    };
};

self.onmessage = async (event) => {
    try {
        const result = await computeKde(event.data);
        self.postMessage({ id: event.data.id, result });
    } catch (error) {
        self.postMessage({ id: event.data.id, error: error.message || 'Browser KDE computation failed' });
    }
};
