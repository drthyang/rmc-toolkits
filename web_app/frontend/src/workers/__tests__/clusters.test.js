// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode cluster engine against the Python reference
// (rmc_toolkits/clusters.py) via the golden fixture written by
// tests/generate_clusters_fixture.py, plus the worker's request checks.

import { describe, expect, it } from 'vitest';
import {
    CHANGE_KINDS, CLUSTER_BECAUSE, CLUSTER_RULES, EVIDENCE_AGAINST, SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE,
    clusterDistortion, evidenceClears, isSplitSite, mix32, shufflePermutation, splitSiteEvidence
} from '../clusters';
import { siteDisplacementsFromRmc6f, siteEllipsoids } from '../pcaKde';
import { handlePcaMessage } from '../pcaKdeWorker';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import fixture from '../../__tests__/fixtures/clusters_fixture.json';

// Integers, strings and booleans must match exactly; every other float to
// 1e-9 (summation order differs between numpy and the port).

const expectClose = (actual, expected, path) => {
    if (expected === null || typeof expected !== 'object') {
        if (typeof expected === 'number' && !Number.isInteger(expected)) {
            expect(actual, path).toBeCloseTo(expected, 9);
        } else {
            expect(actual, path).toEqual(expected);
        }
        return;
    }
    if (Array.isArray(expected)) {
        expect(actual, path).toHaveLength(expected.length);
        expected.forEach((value, index) => expectClose(actual[index], value, `${path}[${index}]`));
        return;
    }
    expect(Object.keys(actual).sort(), path).toEqual(Object.keys(expected).sort());
    Object.entries(expected).forEach(([key, value]) => expectClose(actual[key], value, `${path}.${key}`));
};

describe('Python parity (clusters_fixture.json)', () => {
    it('hashes and permutes like the Python engine', () => {
        Object.entries(fixture.hashes.mix32).forEach(([value, hashed]) => expect(mix32(Number(value))).toBe(hashed));
        const { n, shuffle, template, member, seed, order } = fixture.hashes.permutation;
        expect(shufflePermutation(n, shuffle, template, member, seed)).toEqual(order);
    });

    for (const testCase of fixture.cases) {
        const parsed = siteDisplacementsFromRmc6f(testCase.text);
        testCase.queries.forEach(({ element, window }, index) => {
            it(`${testCase.name}: ${element} ${window.join('–')} Å`, () => {
                expectClose(clusterDistortion(parsed, element, window), testCase.expected[index], 'result');
            });
        });
    }
});

describe('a site at fraction 0 (review 2026-10-12 #1)', () => {
    it('pairs each atom with its partner of the same cell when the site mean is a hair below 0', () => {
        // The rigid V2 dimers' first site averages −0.0001 Å: wrapped to
        // 0.9999 it pulled the partner from the next cell, which hid the
        // dimers' correlation and broke edges.
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'origin-dimer').text);
        expect(parsed.sites[0].siteFractional[0]).toBeGreaterThan(0.99);
        expect(parsed.sites[0].siteCellFractional[0]).toBeCloseTo(-2e-5, 9);
        const result = clusterDistortion(parsed, 'V', [2.0, 2.5]);
        expect(result.templates[0].offsets).toEqual([[0, 0, 0], [0, 0, 0]]);
        expect(result.outliers.brokenEdgeFraction).toBe(0);
        expect(result.z.breathingSd).toBeLessThan(-10);
    });
});

describe('templates with their own edges (review 2026-10-12 #2)', () => {
    it('measures each chain with its own template\'s edges', () => {
        // Chain B (bonded 0-1, 0-2) measured with chain A's edges (0-1, 1-2)
        // counted its 4 Å end-to-end span as an edge: ladderWidth 1.03 Å.
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'chains').text);
        const result = clusterDistortion(parsed, 'Mo', [1.8, 2.2]);
        expect(result.templates.map((template) => template.edges)).toEqual([[[0, 1], [1, 2]], [[0, 1], [0, 2]]]);
        expect(result.observed.ladderWidth).toBeLessThan(0.2);
        expect(result.outliers.brokenEdgeFraction).toBe(0);
    });
});

describe('missing atoms (review 2026-10-12 #16)', () => {
    it('counts the bonds of a cluster that lost an atom as its own, not as contacts between clusters', () => {
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'vacancies').text);
        const result = clusterDistortion(parsed, 'Nb', [2.5, 3.5]);
        expect(result.templates[0].incomplete).toBeGreaterThan(0);
        expect(result.outliers.extraContacts).toBe(0);
    });
});

// Every Evidence record of a verdict, with the finding it belongs to.
const recordsOf = (verdict) => Object.entries(verdict).flatMap(([name, value]) => {
    if (name === 'rules') return [];
    const findings = name === 'window' ? Object.entries(value).map(([key, item]) => [`window.${key}`, item]) : [[name, value]];
    return findings.flatMap(([key, item]) => item.evidence.map((record) => [key, record]));
});
// Every finding as [answer, because], the window's as window.<name>.
const answersOf = (verdict) => Object.fromEntries(Object.entries(verdict).flatMap(([name, value]) => {
    if (name === 'rules') return [];
    if (name === 'window') return Object.entries(value).map(([key, item]) => [`window.${key}`, [item.answer, item.because]]);
    return [[name, [value.answer, value.because]]];
}));
// The Evidence rule as the docs state it, written out apart from evidenceClears.
const clearsAsDocumented = (record) => {
    const entries = Object.entries(record.bar);
    return entries.length > 0 && entries.every(([field, bar]) => {
        const value = record[field];
        if (value == null) return false;
        const pass = (a, b) => (record.strict ? a > b : a >= b);
        if (record.side === 'above') return pass(value, bar);
        if (record.side === 'below') return pass(-value, bar);
        return pass(Math.abs(value), bar);
    });
};

describe('the verdict (the engine decides; clusters.py is the reference)', () => {
    it('pins the closed vocabularies and the rules to the Python engine\'s (the fixture\'s enums)', () => {
        expect(CLUSTER_BECAUSE).toEqual(fixture.enums.because);
        expect([...EVIDENCE_AGAINST]).toEqual(fixture.enums.against);
        expect([...CHANGE_KINDS]).toEqual(fixture.enums.changeKinds);
        expect({ ...CLUSTER_RULES }).toEqual(fixture.enums.rules);
    });

    it('every fixture verdict: because from its closed set, clears recomputable from the record, clearing evidence first', () => {
        let records = 0;
        fixture.cases.forEach((testCase) => testCase.expected.forEach(({ verdict }) => {
            Object.entries(answersOf(verdict)).forEach(([name, [answer, because]]) => {
                expect([true, false, null]).toContain(answer);
                expect(CLUSTER_BECAUSE[name], name).toContain(because);
            });
            recordsOf(verdict).forEach(([name, record]) => {
                records += 1;
                expect(record.clears, `${testCase.name} ${name} ${record.statistic}`).toBe(clearsAsDocumented(record));
                expect(record.clears).toBe(evidenceClears(record));
                expect(EVIDENCE_AGAINST).toContain(record.against);
            });
            Object.entries(verdict).filter(([name]) => !['rules', 'window'].includes(name)).forEach(([name, item]) => {
                const clears = item.evidence.map((record) => record.clears);
                expect(clears, name).toEqual([...clears].sort((a, b) => Number(b) - Number(a)));
            });
        }));
        expect(records).toBeGreaterThan(100);
    });

    it('the split rule: one definition, the PC1 kurtosis on a resolved axis at or below −0.5 and three sampling errors', () => {
        const row = (kappa, count = 4000, resolved = true) => ({ referenceNumber: 1, count, axisResolved: [resolved, true, true], excessKurtosis: [kappa, 0, 0] });
        expect([SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE]).toEqual([-0.5, 3]);
        expect(isSplitSite(row(SPLIT_SITE_KURTOSIS - 0.01))).toBe(true);
        expect(isSplitSite(row(SPLIT_SITE_KURTOSIS + 0.01))).toBe(false);
        // 216 copies: three sampling errors are 3·√(24/216) = 1.
        expect(isSplitSite(row(-0.99, 216))).toBe(false);
        expect(isSplitSite(row(-1.01, 216))).toBe(true);
        expect(isSplitSite(row(-1.5, 4000, false))).toBe(false);
        expect(splitSiteEvidence(row(-1.5, 4000, false)).z).toBeNull();
        expect(isSplitSite(row(null))).toBe(false);
        expect(splitSiteEvidence(row(-1.01, 216))).toMatchObject({ against: 'siteKurtosis', side: 'below', bar: { z: 3, change: 0.5 } });
        // The rule list_sites applied before it moved here, away from the bounds.
        for (const kappa of [-2, -1.2, -0.8, -0.51, -0.49, -0.2, 0.3]) {
            for (const count of [64, 216, 512, 4000]) {
                const old = kappa <= SPLIT_SITE_KURTOSIS && kappa <= -SPLIT_SITE_SE * Math.sqrt(24 / count);
                expect(isSplitSite(row(kappa, count)), `${kappa} ${count}`).toBe(old);
            }
        }
    });

    it('split F6 octahedra: every F site split (the one rule, on the engine\'s own site rows), LESS distorted, the pattern shows', () => {
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'split-octahedra').text);
        const { verdict } = clusterDistortion(parsed, 'F', [2.5, 3.25]);
        expect(verdict.splitMemberSites).toMatchObject({ answer: true, because: 'half_or_more_split', sites: [2, 3, 4, 5, 6, 7] });
        const rows = siteEllipsoids(parsed.sites);
        verdict.splitMemberSites.evidence.forEach((record) => {
            expect(record).toEqual(splitSiteEvidence(rows.find((row) => row.referenceNumber === record.site)));
        });
        expect(answersOf(verdict)).toMatchObject({
            distorted: [true, 'axis_family'], lessDistorted: [true, 'below_null'], rigid: [false, 'split_sites']
        });
    });
});

describe('the verdict on the validation\'s dev cases (the decisions the tool used to make)', () => {
    const verdictOf = (id, element, window) => {
        const built = buildCase(DEV_CASES.find((item) => item.id === id), { partials: false });
        return clusterDistortion(siteDisplacementsFromRmc6f(built.rmc6f), element, window);
    };

    it('dev-breathing: the Mo3 triangles breathe; against independent atoms the edge spread is far below (z < −3), at their own size it is not', () => {
        const result = verdictOf('dev-breathing', 'Mo', [2.4, 3.21]);
        expect(result.z.ladderWidth).toBeLessThan(-3);
        expect(Math.abs(result.verdict.lessDistorted.evidence[0].z)).toBeLessThan(3);
        expect(answersOf(result.verdict)).toMatchObject({
            breathing: [true, 'above_null'], distorted: [false, 'within_null'], lessDistorted: [false, 'within_null'],
            rigid: [false, 'within_null'], splitMemberSites: [false, 'none_split']
        });
        expect(result.verdict.breathing.evidence[0].change).toBeGreaterThan(0.3);
    });

    it('dev-jahn-teller: the F6 octahedra are distorted (an axis family), their F sites split, so not rigid', () => {
        const { verdict } = verdictOf('dev-jahn-teller', 'F', [2.5, 3.25]);
        expect(answersOf(verdict)).toMatchObject({
            distorted: [true, 'axis_family'], lessDistorted: [true, 'below_null'], splitMemberSites: [true, 'half_or_more_split'],
            rigid: [false, 'split_sites']
        });
        expect(verdict.splitMemberSites.sites.length).toBeGreaterThanOrEqual(3);
        const favoured = verdict.favouredAxes.families.filter((family) => family.answer).map((family) => family.family);
        expect(favoured).toContain('100');
    });

    it('dev-cell-wrap: the V2 dimers at fraction 0 are rigid (their size far stiffer than independent atoms)', () => {
        const { verdict } = verdictOf('dev-cell-wrap', 'V', [2.17, 2.35]);
        expect(answersOf(verdict)).toMatchObject({
            rigid: [true, 'stiffer_bonds'], breathing: [false, 'below_null'], lessDistorted: [null, 'single_edge']
        });
    });
});

describe('the worker\'s clusters request', () => {
    const text = fixture.cases[0].text;
    const ask = (params) => handlePcaMessage({ kind: 'clusters', ...params }, async () => text);

    it('answers with the engine payload', async () => {
        const result = await ask({ element: 'nb', rMin: '2.5', rMax: 3.5 });
        expect(result.element).toBe('Nb');
        expect(result.kind).toMatchObject({ members: 4, edges: 6, clusters: 64 });
    });

    it('rejects missing, inverted and oversized windows with the route\'s messages', async () => {
        await expect(ask({ element: 'Nb', rMin: 2.5 })).rejects.toThrow('rMax required');
        await expect(ask({ element: 'Nb', rMin: 2.5, rMax: 20 })).rejects.toThrow('rMax must be <= 8, got 20');
        await expect(ask({ element: 'Nb', rMin: 3.5, rMax: 2.5 })).rejects.toThrow('window must be 0 <= min < max');
        await expect(ask({ element: 'Xx', rMin: 2.5, rMax: 3.5 })).rejects.toThrow('no site holds Xx');
    });
});
