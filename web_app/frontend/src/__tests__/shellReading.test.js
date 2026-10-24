// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One shell reading per pair, per run (architecture candidate #6): on the
// Demo and on every dev case of the validation, each reader shows the run's
// shell catalogue (llm/context/shells.js) and nothing of its own — the run
// context's pair_correlations, find_gr_peaks, the Workbench's clarifier, the
// verification's evidence gate and arm A+'s cluster windows. Before, they searched different
// r ranges, rounded differently and integrated the CN differently, and the
// Demo's Ta–Se first peak read 3.9 Se per Ta in one and 3.86 in another; it
// is 3.25 now, to the minimum between the two Ta–Se peaks (approved
// 2026-10-19). The held-out test cases are never built here.

import { beforeAll, describe, expect, it } from 'vitest';
import { buildClarifier } from '../llm/copilot/clarify';
import { buildRunTools } from '../llm/copilot/runTools';
import { createVerification } from '../llm/copilot/verification';
import { buildRunContext } from '../llm/context/runContext';
import { CLUSTER_CN_MAX, CN_MATCH_TOLERANCE, clusterElementsOf, networkCentresOf, POLYHEDRON_COUNTS, shellCatalogue } from '../llm/context/shells';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { loadRun } from '../../validation/runFolder';
import { runArmSweep } from '../../validation/armA';
import { loadDemoAssistantRun } from './modelVisible/helpers';

const sig = (value, digits = 3) => (Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : value);
const RUNS = ['demo', ...DEV_CASES.map((item) => item.id)];

describe('one shell reading per pair, per run', () => {
    const runs = new Map();

    beforeAll(async () => {
        const demo = await loadDemoAssistantRun();
        runs.set('demo', { props: demo.run });
        for (const item of DEV_CASES) {
            const run = await loadRun(buildCase(item), { ladder: false });
            runs.set(item.id, { props: run.props, copilotRun: run });
        }
    }, 120000);

    it.each(RUNS)('%s: the context, find_gr_peaks, the clarifier and the evidence gate read the catalogue', async (name) => {
        const { props } = runs.get(name);
        const catalogue = shellCatalogue(props);
        expect(catalogue).not.toBeNull();
        // The clarifier's catalogue is the run's (memoised), and the gate's
        // finite clusters are its: "no order" asks analyze_clusters for each.
        const clarifier = buildClarifier(props);
        expect(clarifier.catalog).toBe(catalogue);
        const clusters = clusterElementsOf(clarifier.catalog);
        expect(clusters).toEqual(catalogue.elements.filter((element) => {
            const cn = catalogue.shell(element, element)?.cn;
            return Number.isFinite(cn) && cn >= CN_MATCH_TOLERANCE && cn <= CLUSTER_CN_MAX;
        }));
        const sendBack = createVerification({ catalog: clarifier.catalog })
            .judge('There is no local order beyond thermal motion.', { canSendBack: true }).sendBack;
        expect(sendBack.missing.filter((item) => item.startsWith('analyze_clusters for ')))
            .toEqual(clusters.map((element) => `analyze_clusters for ${element} (its first ${element}–${element} shell is a finite cluster's)`));

        const context = Object.fromEntries((buildRunContext(props).pair_correlations || []).map((entry) => [entry.pair, entry]));
        const grPeaks = buildRunTools(props).find((tool) => tool.name === 'find_gr_peaks');
        let compared = 0;
        for (const shell of catalogue.pairs) {
            const { centre: a, partner: b } = shell;
            const entry = context[shell.label];
            expect(entry, `${name} ${shell.label} in the context`).toBeDefined();
            expect(entry.gr_peaks_A.map((peak) => peak.r)).toEqual(shell.peaks.slice(0, 2).map((peak) => sig(peak.r)));
            const { data } = await grPeaks.run({ pair: shell.label });
            if (!shell.window) {
                expect(entry.first_shell).toBeUndefined();
                expect(data.first_shell).toBeUndefined();
                continue;
            }
            compared += 1;
            // The context.
            expect(entry.first_shell.to_A).toBe(sig(shell.window[1]));
            expect(entry.first_shell[`cn_${b}_per_${a}`]).toBe(sig(shell.cn));
            expect(entry.first_shell.first_peak_to_A).toBe(shell.firstPeak ? sig(shell.firstPeak.window[1]) : undefined);
            expect(entry.first_shell.first_peak_cn).toBe(shell.firstPeak ? sig(shell.firstPeak.cn) : undefined);
            expect('gap_g_fraction' in entry.first_shell).toBe(!shell.separated);
            expect(entry.cut_at_min_distance?.at_A).toBe(shell.constraint?.limited ? shell.constraint.limitA : undefined);
            // find_gr_peaks.
            expect(data.first_shell.window_A).toEqual(shell.window);
            expect(data.first_shell.peaks_in_shell_A).toEqual(shell.peaksInShell.map((r) => sig(r, 4)));
            expect(data.first_shell.cn[`${b}_per_${a}`]).toBe(sig(shell.cn));
            if (a !== b) expect(data.first_shell.cn[`${a}_per_${b}`]).toBe(sig(shell.reverseCn));
            expect(data.first_shell.cn_if_cut_0p2A_earlier_later).toEqual([shell.cut.earlier.cn, shell.cut.later.cn].map((cn) => sig(cn)));
            expect(data.first_shell.first_peak_only?.window_A).toEqual(shell.firstPeak?.window);
            expect(data.first_shell.first_peak_only?.cn?.[`${b}_per_${a}`]).toBe(shell.firstPeak ? sig(shell.firstPeak.cn) : undefined);
            // The split by number (plan item 4.2): only for a joined shell.
            expect(Boolean(data.first_shell.split)).toBe(shell.joined);
            if (shell.split) {
                expect(data.first_shell.split.parts.map((part) => part.window_A)).toEqual(shell.split.parts.map((part) => part.window));
                expect(data.first_shell.split.parts.map((part) => part.mean_d_A)).toEqual(shell.split.parts.map((part) => sig(part.meanR, 4)));
                expect(data.first_shell.split.long_minus_short_A).toBe(sig(shell.split.longMinusShort));
            }
            expect(data.min_distance_constraint_A).toBe(shell.constraint?.limitA);
            expect(data.g_at_constraint_fraction).toBe(shell.constraint?.limited ? sig(shell.constraint.gFraction, 2) : undefined);
            expect(data.first_peak_window_A).toBeUndefined();
        }
        expect(compared).toBeGreaterThan(0);
    });

    // Since candidate 8 A+ tries the evidence gate's cluster elements (the
    // element's own first shell holds 0.5–4.5 neighbours: clusterElementsOf,
    // which the gate's "no order" asks analyze_clusters for, above), not every
    // same-element shell peaking below 3.6 Å.
    it.each(DEV_CASES.map((item) => item.id))('%s: arm A+ runs the cluster engine in the catalogue\'s windows', async (name) => {
        const { props, copilotRun } = runs.get(name);
        const catalogue = shellCatalogue(props);
        const { readings } = await runArmSweep(copilotRun);
        const windows = readings.filter((reading) => reading.page === 'analyze_clusters').map((reading) => reading.what);
        const expected = clusterElementsOf(catalogue)
            .map((element) => catalogue.shell(element, element))
            .map((shell) => `${shell.centre} clusters in ${shell.window[0]}–${shell.window[1]} Å (engine verdict)`);
        expect(windows).toEqual(expected);
    }, 60000);

    // The local-order interview (2026-10-24): "no order" needs
    // polyhedral_distortion on each centre of octahedra or tetrahedra no
    // finite cluster holds (networkCentresOf): its nearest first shell of
    // another element holds 4 or 6, its own holds no finite cluster.
    it('the centres of polyhedra no finite cluster holds: the Demo\'s GaSe₄, a perovskite\'s BO₆, both rock-salt and zincblende sites, no cluster element', () => {
        expect(POLYHEDRON_COUNTS).toEqual([4, 6]);
        const centres = (name) => networkCentresOf(shellCatalogue(runs.get(name).props));
        // The Demo: Ga (4 Se, Ga–Ga twelve); not Ta, a cluster element whose TaSe₆ octahedra analyze_clusters reaches; not Se (½ Ga).
        expect(centres('demo')).toEqual(['Ga']);
        expect(clusterElementsOf(shellCatalogue(runs.get('demo').props))).toEqual(['Ta']);
        // Perovskites: the B site (BO₆; B–B six), not the A site (AO₁₂) nor O (two B).
        for (const [name, centre] of [['dev-zone-boundary', 'Ti'], ['dev-thermal-null', 'Ti'], ['dev-offcentre', 'Nb'], ['dev-symmetry-lowering', 'Ta']]) {
            expect(centres(name), name).toEqual([centre]);
        }
        // MgO (rock salt) and CuBr (zincblende): both sites; isolated CuF₆ (each F bonds one Cu): Cu, F's own F–F shell a cluster's.
        expect(centres('dev-fit-unconverged')).toEqual(['Mg', 'O']);
        expect(centres('dev-chirality-trap')).toEqual(['Br', 'Cu']);
        expect(centres('dev-jahn-teller')).toEqual(['Cu']);
        // Cluster elements only (Mo₃, W₆): no centre outside them.
        expect(centres('dev-breathing')).toEqual([]);
        expect(centres('dev-octahedral-clusters')).toEqual([]);
        // FeSi (seven Si): no octahedron or tetrahedron.
        expect(centres('dev-chiral')).toEqual([]);
        expect(networkCentresOf(null)).toEqual([]);
    });

    // Plan item 4.2: a joined shell's split by number. The Demo joins no
    // shell (its Ta–Se and Se–Se peaks overlap without a deep trough), so it
    // has none; dev-jahn-teller's Cu–F shell is the planted 4 + 2.
    it('the split by number: none on the Demo; dev-jahn-teller\'s Cu–F is 4 short + 2 long, long − short within two sampling errors of the planted 0.32 Å', () => {
        const demo = shellCatalogue(runs.get('demo').props);
        const demoShells = demo.elements.flatMap((centre) => demo.elements.map((partner) => demo.shell(centre, partner))).filter(Boolean);
        expect(demoShells.length).toBeGreaterThan(0);
        expect(demoShells.filter((shell) => shell.joined || shell.split)).toEqual([]);
        const { split } = shellCatalogue(runs.get('dev-jahn-teller').props).shell('Cu', 'F');
        expect(split.direction).toBe('forward');
        expect(split.parts.map((part) => sig(part.cn, 2))).toEqual([4, 2]);
        expect(split.parts.map((part) => sig(part.share, 2))).toEqual([0.66, 0.34]);
        expect(Math.abs(split.longMinusShort - 0.32)).toBeLessThan(2 * split.longMinusShortSe);
        // 216 Cu, 0.07 Å spreads: ± 0.0043 Å.
        expect(sig(split.longMinusShortSe, 2)).toBe(0.0043);
    });

    it('the Demo\'s Ta–Se first peak: 3.25 Se per Ta to the inter-peak minimum (was 3.86–3.9 by r ± FWHM)', () => {
        const catalogue = shellCatalogue(runs.get('demo').props);
        const shell = catalogue.shell('Ta', 'Se');
        expect(shell.peaksInShell).toEqual([2.52, 2.74]);
        expect(shell.window).toEqual([2.21, 3.23]);
        expect(sig(shell.cn)).toBe(6);
        expect(shell.firstPeak.window).toEqual([2.21, 2.65]);
        expect(sig(shell.firstPeak.cn)).toBe(3.25);
        // The window r ± FWHM gave (2.33–2.71 Å) holds 3.86: it reaches into the long bonds.
        expect(sig(catalogue.window('Ta', 'Se', [2.33, 2.71]).cn)).toBe(3.86);
        expect(catalogue.window('Ta', 'Se', [2.33, 2.71]).cnMargin.cut).toBe(true);
    });
});
