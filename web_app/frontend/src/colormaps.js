// Minimal colormap lookup tables for rendering KDE density grids on a canvas.
// Each map is a list of [r, g, b] anchors interpolated linearly over [0, 1].

const ANCHORS = {
    viridis: [
        [68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]
    ],
    magma: [
        [0, 0, 4], [81, 18, 124], [183, 55, 121], [252, 137, 97], [252, 253, 191]
    ],
    seismic: [
        [0, 0, 76], [0, 0, 255], [255, 255, 255], [255, 0, 0], [127, 0, 0]
    ],
    reds: [
        [255, 245, 240], [252, 187, 161], [251, 106, 74], [203, 24, 29], [103, 0, 13]
    ],
    greys: [
        [12, 13, 15], [70, 74, 80], [130, 136, 144], [195, 200, 207], [245, 247, 250]
    ]
};

export const COLORMAP_NAMES = Object.keys(ANCHORS);

const LUT_SIZE = 256;
const lutCache = {};

function buildLut(name) {
    const anchors = ANCHORS[name] || ANCHORS.viridis;
    const lut = new Uint8ClampedArray(LUT_SIZE * 3);
    const segments = anchors.length - 1;
    for (let i = 0; i < LUT_SIZE; i += 1) {
        const t = i / (LUT_SIZE - 1);
        const scaled = t * segments;
        const lower = Math.min(segments, Math.floor(scaled));
        const upper = Math.min(segments, lower + 1);
        const frac = scaled - lower;
        for (let channel = 0; channel < 3; channel += 1) {
            lut[i * 3 + channel] = anchors[lower][channel] + (anchors[upper][channel] - anchors[lower][channel]) * frac;
        }
    }
    return lut;
}

export function getLut(name) {
    if (!lutCache[name]) {
        lutCache[name] = buildLut(name);
    }
    return lutCache[name];
}

// Map a normalized value in [0, 1] to an [r, g, b] triple for the given map.
export function sampleColormap(name, value) {
    const lut = getLut(name);
    const index = Math.max(0, Math.min(LUT_SIZE - 1, Math.round(value * (LUT_SIZE - 1)))) * 3;
    return [lut[index], lut[index + 1], lut[index + 2]];
}
