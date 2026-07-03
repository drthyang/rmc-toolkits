import { describe, expect, it } from 'vitest';
import { buildReportMarkdown } from '../report/buildReport';

const fullContext = {
    run: '5K_try1',
    structure: {
        cell_A: [10.02, 10.02, 10.02],
        angles_deg: [90, 90, 90],
        supercell: [6, 6, 6],
        total_atoms: 12288,
        elements: { Ga: 864, Se: 6912 },
        space_group: 'F-43m (No. 216)'
    },
    datasets: [
        { kind: 'neutron_sq', title: 'S(Q) (neutron)', rwp: 0.0412 },
        { kind: 'npdf', title: 'PDF1', rwp: 0.0533 }
    ],
    convergence: {
        n_steps: 500,
        first: 2.3,
        last: 0.87,
        min: 0.86,
        max: 2.3,
        recent_slope_per_step: -1.2e-5,
        final_chi_squared: 2.38
    }
};

describe('buildReportMarkdown', () => {
    it('renders every deterministic section with the data', () => {
        const markdown = buildReportMarkdown(fullContext, null, { generatedAt: new Date('2026-07-03T00:00:00Z') });
        expect(markdown).toContain('# RMCProfile run report — 5K_try1');
        expect(markdown).toContain('## Model summary');
        expect(markdown).toContain('| Cell (Å) | 10.02 × 10.02 × 10.02 |');
        expect(markdown).toContain('| Detected space group | F-43m (No. 216) |');
        expect(markdown).toContain('## Fit metrics');
        expect(markdown).toContain('| S(Q) (neutron) | neutron_sq | 0.0412 |');
        expect(markdown).toContain('## Convergence');
        expect(markdown).toContain('| ln(χ²) first → last | 2.3 → 0.87 |');
        expect(markdown).toContain('2026-07-03');
    });

    it('labels the AI narrative with the model and a disclaimer', () => {
        const markdown = buildReportMarkdown(fullContext, 'The run appears converged.', {
            model: 'llama3.2',
            temperature: 0.2,
            generatedAt: new Date('2026-07-03T00:00:00Z')
        });
        expect(markdown).toContain('## AI assessment');
        expect(markdown).toContain('`llama3.2`');
        expect(markdown).toContain('verify every statement');
        expect(markdown).toContain('The run appears converged.');
        expect(markdown).toContain('temperature 0.2');
    });

    it('degrades gracefully with an empty context and no narrative', () => {
        const markdown = buildReportMarkdown({ run: 'bare' }, null, { generatedAt: new Date('2026-07-03T00:00:00Z') });
        expect(markdown).toContain('_No model structure detected._');
        expect(markdown).toContain('_No fitted datasets detected._');
        expect(markdown).toContain('_No R-value / chi-squared history detected._');
        expect(markdown).not.toContain('## AI assessment');
    });
});
