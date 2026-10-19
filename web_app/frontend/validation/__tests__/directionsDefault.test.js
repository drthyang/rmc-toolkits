// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 2.5, the known answer: the Displacement Directions page and
// analyze_displacement_directions at the same settings. On the dev thermal
// case (thermal motion only, its X flattened across its bond) every site's
// readout at the page's defaults is the tool's result, number for number, and
// arm A reads exactly that. Under the page's old default (ν = 10, the
// readouts of the 2× smoothed map) the flattened X's map test stayed below
// the 3σ bar that the engine's tiling clears. Dev case only.

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../cases';
import { loadRun } from '../runFolder';
import { ARM_A_RULES, directionsOf } from '../armA';
import { summarizeOrientation } from '../../src/copilotTools';
import { unitCellVectors } from '../../src/pcaCrystalFrame';
import { DIRECTION_DEFAULTS, directionMapRequest, directionReadoutRequest } from '../../src/directionSettings';

const loadThermal = () => loadRun(buildCase(DEV_CASES.find((item) => item.id === 'dev-thermal-null')), { ladder: false });
const atDefaults = (run, referenceNumber, request = directionReadoutRequest) => run.engine('orientation', {
    ...request({ ...DIRECTION_DEFAULTS, referenceNumber }), geometry: false
});

describe('the page and the tool at the same settings (plan item 2.5): the dev thermal case', () => {
    it('every site: the page\'s default readout is analyze_displacement_directions\' result, and arm A reads it', async () => {
        const run = await loadThermal();
        const { structure } = run.props;
        const unitCell = unitCellVectors(structure.latticeVectors, structure.supercell);
        const tool = run.tools.find((item) => item.name === 'analyze_displacement_directions');
        const { sites } = await run.engine('sites');
        expect(sites.length).toBeGreaterThan(0);
        for (const site of sites) {
            const ref = site.referenceNumber;
            const page = await atDefaults(run, ref);
            const answer = await tool.run({ site_ref: ref });
            // The engine's own tiling.
            expect(page.frequency, `site ${ref}`).toBe(page.recommendedFrequency);
            // Everything the tool reports (the peak, its direction and σ, the
            // map test, the ± asymmetry, the tensor), from the page's readout.
            expect(summarizeOrientation(page, unitCell).data, `site ${ref}`).toEqual(answer.data);
            // And arm A reads the same numbers off the page (directionsOf takes
            // the run as the arms read it: its engine and structure).
            const readings = [];
            await directionsOf({ engine: run.engine, structure }, site, readings);
            expect(readings[0].value, `site ${ref}`).toEqual({
                frequency: answer.data.frequency,
                peak_sigma: answer.data.peak.significance_sigma,
                enhancement: answer.data.peak.enhancement,
                peak_uvw: answer.data.peak.direction_uvw,
                radius_deg: answer.data.peak.radius_deg,
                map_sigma: answer.data.map_significance_sigma
            });
        }
    });

    it('the flattened X is anisotropic on the engine\'s tiling, and was below the bar on the old page default', async () => {
        const run = await loadThermal();
        const { sites } = await run.engine('sites');
        const flattened = sites.filter((site) => site.element === 'O');
        expect(flattened.length).toBeGreaterThan(0);
        for (const site of flattened) {
            const page = await atDefaults(run, site.referenceNumber);
            const old = await run.engine('orientation', { referenceNumber: site.referenceNumber, frequency: 10, smoothing: 2, geometry: false });
            expect(page.mapSignificance, `site ${site.referenceNumber}`).toBeGreaterThanOrEqual(ARM_A_RULES.SIGMA);
            expect(old.mapSignificance, `site ${site.referenceNumber}`).toBeLessThan(ARM_A_RULES.SIGMA);
            // Thermal motion: no preferred direction on either.
            expect(page.peakSignificance).toBeLessThan(ARM_A_RULES.SIGMA);
        }
    });

    it('smoothing moves only what is drawn: the drawn map\'s peak can move, the raw-count tests never do', async () => {
        const run = await loadThermal();
        const { sites } = await run.engine('sites');
        let moved = 0;
        // The map opens unsmoothed (the drawn map is then the readout itself);
        // a user who smooths it, two passes here, redraws only the map.
        const smoothed = (settings) => directionMapRequest({ ...settings, smoothing: 2 });
        expect(directionMapRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 1 })).toEqual(directionReadoutRequest({ ...DIRECTION_DEFAULTS, referenceNumber: 1 }));
        for (const site of sites) {
            const readout = await atDefaults(run, site.referenceNumber);
            const drawn = await atDefaults(run, site.referenceNumber, smoothed);
            expect(drawn.smoothing).toBe(2);
            // The same tiling and the same raw-count tests.
            expect(drawn.frequency).toBe(readout.frequency);
            expect(drawn.counts).toEqual(readout.counts);
            expect(drawn.mapSignificance).toBe(readout.mapSignificance);
            expect(drawn.antipodalAsymmetryZ).toBe(readout.antipodalAsymmetryZ);
            expect(drawn.orientationAnisotropySignificance).toBe(readout.orientationAnisotropySignificance);
            // The peak cell is picked on the smoothed map.
            if (drawn.peakCell !== readout.peakCell) moved += 1;
        }
        // On this case the smoothed map's peak sits in another cell for at
        // least one site: the readouts must not come from the drawn map.
        expect(moved).toBeGreaterThan(0);
    });
});
