// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { buildZip } from './zipArchive';

// Export figures from the dashboard and KDE pages.
//
// Charts are inline <svg>, so they export as PNG (raster) or SVG (true vector).
// The SVG path clones the element and freezes the *computed* value of a small
// set of presentation properties onto each node — this resolves CSS custom
// properties to concrete colors and makes the file self-contained, since the
// charts are otherwise styled through external CSS classes.
//
// The KDE/3D panels render to <canvas> (a pixel heatmap and a WebGL scene), so
// they export as PNG only — at native or higher resolution.

// Presentation properties that actually affect how a chart renders. Keeping the
// list tight avoids bloating the serialized markup with irrelevant declarations.
const STYLE_PROPS = [
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-width',
    'stroke-dasharray',
    'stroke-linejoin',
    'stroke-linecap',
    'opacity',
    'font-family',
    'font-size',
    'font-weight',
    'font-style',
    'letter-spacing',
    'text-anchor',
    'font-variant-numeric',
    'vector-effect',
];

// Walk the source tree and the clone in lockstep, freezing computed styles.
const inlineComputedStyles = (source, target) => {
    const computed = window.getComputedStyle(source);
    let declaration = '';
    for (const property of STYLE_PROPS) {
        const value = computed.getPropertyValue(property);
        if (value) declaration += `${property}:${value};`;
    }
    if (declaration) target.setAttribute('style', declaration);

    const sourceChildren = source.children;
    const targetChildren = target.children;
    for (let index = 0; index < sourceChildren.length; index += 1) {
        inlineComputedStyles(sourceChildren[index], targetChildren[index]);
    }
};

// Clone an inline <svg>, inline its computed styles, and stamp explicit
// dimensions so it stands on its own outside the page.
const standaloneSvg = (svgElement) => {
    const box = svgElement.viewBox?.baseVal;
    const width = box?.width || svgElement.clientWidth || 720;
    const height = box?.height || svgElement.clientHeight || 450;
    const clone = svgElement.cloneNode(true);
    inlineComputedStyles(svgElement, clone);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('width', width);
    clone.setAttribute('height', height);
    return { clone, width, height };
};

const loadImage = (src) => new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error('Could not rasterize the figure'));
    image.src = src;
});

// Turn a filename-ish string into something safe across file systems.
export const sanitizeFilename = (name) => {
    const cleaned = (name || 'figure')
        .replace(/\^\{([^}]+)\}/g, '$1')
        .replace(/[^\w.-]+/g, '_')
        .replace(/^_+|_+$/g, '');
    return cleaned || 'figure';
};

// Prompt the browser to download a Blob under `filename`.
export const downloadBlob = (blob, filename) => {
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    window.setTimeout(() => URL.revokeObjectURL(url), 0);
};

// --- SVG charts ------------------------------------------------------------

// Rasterize an inline <svg> element into a PNG Blob at `scale`x its viewBox.
export const svgToPngBlob = async (svgElement, { scale = 2, background = '#ffffff' } = {}) => {
    const { clone, width, height } = standaloneSvg(svgElement);
    const serialized = new XMLSerializer().serializeToString(clone);
    const svgUrl = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(serialized)}`;

    const image = await loadImage(svgUrl);
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(width * scale);
    canvas.height = Math.round(height * scale);
    const context = canvas.getContext('2d');
    if (background) {
        context.fillStyle = background;
        context.fillRect(0, 0, canvas.width, canvas.height);
    }
    context.drawImage(image, 0, 0, canvas.width, canvas.height);

    return new Promise((resolve, reject) => {
        canvas.toBlob(
            (blob) => (blob ? resolve(blob) : reject(new Error('Could not encode the figure'))),
            'image/png'
        );
    });
};

// Serialize an inline <svg> into a standalone, self-styled SVG document string.
export const svgToSvgString = (svgElement) => {
    const { clone } = standaloneSvg(svgElement);
    const serialized = new XMLSerializer().serializeToString(clone);
    return `<?xml version="1.0" encoding="UTF-8"?>\n${serialized}`;
};

// Convenience: rasterize an inline <svg> and save it as `<name>.png`.
export const saveSvgAsPng = async (svgElement, name, options) => {
    const blob = await svgToPngBlob(svgElement, options);
    downloadBlob(blob, `${sanitizeFilename(name)}.png`);
    return blob;
};

// Convenience: serialize an inline <svg> and save it as a vector `<name>.svg`.
export const saveSvgAsSvg = (svgElement, name) => {
    const blob = new Blob([svgToSvgString(svgElement)], { type: 'image/svg+xml;charset=utf-8' });
    downloadBlob(blob, `${sanitizeFilename(name)}.svg`);
    return blob;
};

// Save an inline <svg> in the requested format ('png' | 'svg').
export const saveSvgFigure = (svgElement, name, format = 'png') => (
    format === 'svg'
        ? saveSvgAsSvg(svgElement, name)
        : saveSvgAsPng(svgElement, name)
);

// Bundle several charts into one .zip download (one entry per figure) so a
// "save all" is a single file instead of N separately-blocked downloads.
// figures: [{ svgElement, name }]; format is 'png' | 'svg'.
export const saveSvgFiguresAsZip = async (figures, format, zipName) => {
    const used = new Map();
    const entries = [];
    for (const { svgElement, name } of figures) {
        // Disambiguate repeated titles (e.g. two "EXAFS Q-space" charts).
        let base = sanitizeFilename(name);
        const seen = used.get(base) || 0;
        used.set(base, seen + 1);
        if (seen > 0) base = `${base}_${seen + 1}`;

        if (format === 'svg') {
            entries.push({ name: `${base}.svg`, data: new TextEncoder().encode(svgToSvgString(svgElement)) });
        } else {
            const blob = await svgToPngBlob(svgElement);
            entries.push({ name: `${base}.png`, data: new Uint8Array(await blob.arrayBuffer()) });
        }
    }
    downloadBlob(buildZip(entries), zipName);
};

// --- Canvas figures (KDE slice, slab, three.js model) ----------------------

// Rasterize a <canvas> element to a PNG Blob.
export const canvasToPngBlob = (canvas) => new Promise((resolve, reject) => {
    canvas.toBlob(
        (blob) => (blob ? resolve(blob) : reject(new Error('Could not encode the figure'))),
        'image/png'
    );
});

// Convenience: save a <canvas> (or already-rendered offscreen canvas) as PNG.
export const saveCanvasAsPng = async (canvas, name) => {
    const blob = await canvasToPngBlob(canvas);
    downloadBlob(blob, `${sanitizeFilename(name)}.png`);
    return blob;
};
