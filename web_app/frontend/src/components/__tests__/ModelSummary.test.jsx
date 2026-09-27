// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Model information card keeps its explanations out of the page body:
// the Dashboard's "previous read" state is a chip (reason in its title), and
// the parse warning's full sentence sits in the label's ? help rather than a
// "hover for details" sub-line.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import ModelSummary from '../ModelSummary';

const decode = (html) => html
    .replace(/&#x27;/g, "'")
    .replace(/&quot;/g, '"')
    .replace(/&amp;/g, '&');

const WARNING = 'Skipped 2 atom lines that could not be parsed (first at line 17).';

const structure = (parseReport = null) => ({
    source: 'run/model.rmc6f',
    totalAtoms: 8,
    supercell: [1, 1, 1],
    latticeVectors: [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
    elementCounts: { Se: 8 },
    atomIndices: { Se: [1, 2] },
    ...(parseReport ? { parseReport, parseWarning: WARNING } : {}),
});

const render = (props) => decode(renderToStaticMarkup(<ModelSummary showSymmetry={false} {...props} />));

describe('Model information card', () => {
    it('flags a kept previous read with a chip only when stale', () => {
        expect(render({ structure: structure() })).not.toMatch(/previous read/);
        const html = render({ structure: structure(), stale: true });
        expect(html).toMatch(/class="ui-chip ui-chip--warn"[^>]*>previous read<\/span>/);
        expect(html).toMatch(/title="The \.rmc6f is shorter than its header declares \(still being written\?\); showing the previous complete read\."/);
    });

    it('puts the parse warning in the label help, with no "hover for details"', () => {
        const html = render({ structure: structure({ invalidLines: 2, nonFiniteLines: 0, parsedAtoms: 8, coordsOnlyAtoms: 0, declaredAtoms: null }) });
        expect(html).not.toMatch(/hover for details/);
        expect(html).not.toMatch(/header declares/);
        const help = html.match(/<button[^>]*aria-label="Parse warning details"[^>]*>\?<\/button><span[^>]*role="tooltip"[^>]*>([^<]*)<\/span>/);
        expect(help).not.toBeNull();
        expect(help[1]).toBe(WARNING);
        expect(html).toMatch(/2 lines skipped/);
    });

    it('keeps the declared-count sub-line when the header declares a count', () => {
        const html = render({ structure: structure({ invalidLines: 0, nonFiniteLines: 0, parsedAtoms: 6, coordsOnlyAtoms: 0, declaredAtoms: 8 }) })
            .replace(/<!-- -->/g, '');
        expect(html).toMatch(/2 atoms missing/);
        expect(html).toMatch(/header declares 8/);
    });
});
