// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The cards the AI Copilot can point at: one per Dashboard plot, the χ²
// history with logs, Detected SG with a symmetry analysis, and the fixed
// cards of every page — each of which must exist in the page that renders it
// (a data-copilot-card attribute, or a panel's cardId / registered view).

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { PAGE_CARDS, buildCopilotCards } from '../copilotCards';

const src = join(dirname(fileURLToPath(import.meta.url)), '..');
const read = (path) => readFileSync(join(src, path), 'utf8');

describe('buildCopilotCards', () => {
    it('lists the run\'s plots and history before the fixed cards', () => {
        const cards = buildCopilotCards({
            plotFiles: [
                { name: 'GTS_PDFpartials.csv', plotKind: 'pdf_partials', plotData: { title: 'Partial PDFs' } },
                { name: 'GTS_FQ1.csv', plotKind: 'xray_sq' },
                { name: 'GTS-01.log', plotKind: 'r_value' },
                { name: 'notes.txt' }
            ],
            rValueFile: { name: 'GTS-01.log' },
            symmetry: { spaceGroup: 'F-43m' }
        });
        expect(cards.slice(0, 5).map((card) => card.id)).toEqual([
            'model-info', 'detected-sg', 'chi2-history', 'plot:GTS_PDFpartials.csv', 'plot:GTS_FQ1.csv'
        ]);
        expect(cards[3]).toEqual({
            id: 'plot:GTS_PDFpartials.csv',
            page: 'dashboard',
            kind: 'pdf_partials',
            file: 'GTS_PDFpartials.csv',
            marks: true,
            title: 'Partial PDFs (GTS_PDFpartials.csv)',
            about: 'the partial g(r) of every element pair (x: r in Å)'
        });
        expect(cards[4].title).toBe('GTS_FQ1.csv');
        expect(new Set(cards.map((card) => card.id)).size).toBe(cards.length);
    });

    it('leaves out the history without logs and Detected SG without symmetry', () => {
        const ids = buildCopilotCards().map((card) => card.id);
        expect(ids).not.toContain('chi2-history');
        expect(ids).not.toContain('detected-sg');
        expect(ids).toContain('direction-sphere');
    });

    it('names a card every page really renders', () => {
        const sources = [
            'components/ModelSummary.jsx', 'components/Dashboard.jsx', 'components/StructurePage.jsx',
            'components/BondGeometryPage.jsx', 'components/PcaKdePage.jsx', 'components/OrientationPage.jsx',
            'components/OrientationView.jsx'
        ].map(read).join('\n');
        PAGE_CARDS.forEach((card) => {
            const tagged = sources.includes(`data-copilot-card="${card.id}"`) || sources.includes(`cardId="${card.id}"`);
            expect(tagged, card.id).toBe(true);
            if (card.view3d) {
                // Every 3D card joins the view registry under the same id.
                const viewer = sources.includes(`useCopilotViewer('${card.id}'`) || sources.includes(`cardId="${card.id}"`);
                expect(viewer, `${card.id} registers its 3D view`).toBe(true);
            }
        });
        expect(read('components/Dashboard.jsx')).toContain('data-copilot-card={`plot:${file.name}`}');
        expect(read('components/Dashboard.jsx')).toContain('data-copilot-card="chi2-history"');
        // The panels that take a cardId register it.
        ['components/SiteStructurePanel.jsx', 'components/FoldedCellPanel.jsx'].forEach((path) => {
            expect(read(path)).toContain('useCopilotViewer(cardId,');
            expect(read(path)).toContain('data-copilot-card={cardId || undefined}');
        });
    });
});
