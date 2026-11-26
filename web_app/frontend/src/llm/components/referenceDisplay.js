// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Numbers by reference as the dock shows them: the user never sees a raw
// {{r5.…}}. The verification fills an answer's references before it is
// judged and sends back one whose references name no value (copilot/
// verification/repair.js); this is the dock's last line, for every text the
// dock renders through MarkdownMessage:
//   - a reference the turn's results fill is written in, also while the text
//     still streams (the answer arrives filled; the text before it, and the
//     answer while it is written, are filled here from the steps' data);
//   - while streaming, a reference still being written ("{{r5.data.av") is
//     held back until its braces close;
//   - any other reference — one the repair could not mend, or a turn with no
//     repair left — becomes a "number missing" marker: in prose a ref: link
//     MarkdownMessage draws as a marked span naming the reference in its
//     hover, in math \text{…}, in code the words alone. The answer's checks
//     say which numbers are missing (the `reference` finding).
// Pure, no React.

import { fillReferences, isReferenceName } from '../copilot/references';
import { collapseEchoes } from '../copilot/echoedReferences';

/** The link scheme of a missing number's marker (MarkdownMessage draws it, never as a link). */
export const MISSING_SCHEME = 'ref:';
export const MISSING_LABEL = 'number missing';

// Anything in double braces (references.js reads the same).
const BRACED = /\{\{\s*([^{}]{1,120}?)\s*\}\}/g;
// Code (fenced, then inline) and math (display, then inline), in the order
// remark reads them: what a marker inside may look like.
const SEGMENTS = /(```[\s\S]*?(?:```|$)|`[^`\n]*`|\$\$[\s\S]*?\$\$|\$[^$\n]+\$)/g;
// The tail of a streaming text that is a reference still being written.
const OPEN_REFERENCE = /\{\{\s*(?:r\d*|c|ct|ctx)?(?:[.[\]\w<>+-]*)\s*\}?$/;

const marker = (name, kind) => {
    if (kind === 'code') return `[${MISSING_LABEL}]`;
    if (kind === 'math') return `\\text{[${MISSING_LABEL}]}`;
    return `[${MISSING_LABEL}](${MISSING_SCHEME}${encodeURIComponent(name)})`;
};

const kindOf = (segment) => {
    if (segment.startsWith('`')) return 'code';
    if (segment.startsWith('$')) return 'math';
    return 'prose';
};

/**
 * The text as the dock shows it (see the header): `store` the turn's
 * reference store (createReferenceStore, put with each step's data), or
 * null; `streaming` whether the text is still being written.
 */
export const displayReferences = (text, { store = null, streaming = false } = {}) => {
    if (typeof text !== 'string' || !text.includes('{{')) return text;
    let source = text;
    if (streaming) source = source.replace(OPEN_REFERENCE, '');
    const filled = store ? fillReferences(collapseEchoes(source, store), store).text : source;
    if (!filled.includes('{{')) return filled;
    return filled.split(SEGMENTS).map((segment, index) => {
        const kind = index % 2 === 1 ? kindOf(segment) : 'prose';
        return segment.replace(BRACED, (whole, inner) => (isReferenceName(inner) ? marker(inner.trim(), kind) : whole));
    }).join('');
};

/** A plain text (the screen reader's line) with every reference left in it read as the marker's words. */
export const plainReferences = (text) => (typeof text === 'string' && text.includes('{{')
    ? text.replace(BRACED, (whole, inner) => (isReferenceName(inner) ? MISSING_LABEL : whole))
    : text);

/** The reference a missing number's marker names, from its href, or null for any other href. */
export const missingReferenceOf = (href) => {
    if (typeof href !== 'string' || !href.startsWith(MISSING_SCHEME)) return null;
    try {
        return decodeURIComponent(href.slice(MISSING_SCHEME.length));
    } catch {
        return href.slice(MISSING_SCHEME.length);
    }
};
