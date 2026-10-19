// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The shared system prompt's ground rules, split by kind (prompts/system.js,
// decided by the manager on the user's delegation, 2026-10-19): every part of
// every sentence is a definition (what a context field is: the quantity, its
// units, how it is computed, its limits as a measurement), an inference (which
// field to use for a conclusion, or what a pattern is evidence of: the
// method) or conduct (how to answer). With method 'full' the prompt is the one
// it always was (the model-visible snapshots); with 'off' (the validation's
// arm C) it holds the definitions and the conduct only.

import { describe, expect, it } from 'vitest';
import { GROUND_RULE_KINDS, GROUND_RULES, SYSTEM_PROMPT, groundRules, systemPrompt } from '../prompts/system';
import { COPILOT_PROMPT_NO_TOOLS, OUTCOME_RULE, copilotSystemPrompt } from '../prompts/copilot';

const kinds = (parts) => parts.map((part) => part.kind).join(' + ');

describe('the ground rules, part by part', () => {
    it('classifies every part of every sentence, and pins the classification', () => {
        for (const parts of GROUND_RULES) {
            expect(parts.length).toBeGreaterThan(0);
            for (const part of parts) {
                expect(GROUND_RULE_KINDS).toContain(part.kind);
                expect(typeof part.text === 'string' && part.text.length > 0).toBe(true);
                // An inference part sends at most the punctuation that ends what remains.
                if (part.kind === 'inference') expect(part.off).toMatch(/^[.;:, ]*$/);
            }
        }
        // One line per sentence, in order: a change to the split shows here.
        expect(GROUND_RULES.map(kinds)).toEqual([
            'definition + inference', // Rwp: the R-factor | the acceptable range, comparisons relative
            'definition', // the convergence history
            'definition + inference + definition', // watchdog_status | use it | change_over_last_10pct
            'definition', // fitted_data, fitted_in, same_measurement_as
            'definition + inference', // partial_weights | say when a conclusion rests on a weak partial
            'definition + inference', // configuration_optimization | accepted moves as sampling sufficiency
            'definition + inference', // acceptance_ratio | what very low values mean
            'definition + inference', // run_settings | treat questions about the .dat as this block
            'definition + inference', // min_distances_A | a peak at the constraint is constraint-limited
            'definition', // max_move_A, fitted_data
            'definition', // symmetry.sites, mean/max_disp_A
            'definition + inference', // static and thermal mixed | predominantly static at low T
            'inference', // the sites above the others distort most
            'definition + inference + definition', // the ladder | what it approximates and names | noise_floor_A
            'definition', // the average structure folds the cells: no superstructure
            'definition + inference + definition', // pair_correlations | the constraint shapes the peak | gap_g_fraction
            'inference + definition', // compare mean_d_A, not the peak | g(r) divides by r^2
            'inference', // a split or displaced shell is evidence of short-range correlations
            'definition', // pca_displacements
            'definition', // U_iso_A2, rms_axes_A, anisotropy
            'definition + inference + definition + inference', // kurtosis and its sign | physical readings
            'definition + inference', // listed by |non_gaussianity| | call out the large ones
            'conduct', // numbers from the context or a tool, quoted
            'conduct', // no uncertainty no tool gave
            'conduct', // say so when the data are lacking
            'conduct' // concise and cautious
        ]);
    });

    it('method full is the prompt as it always read; off holds every definition and conduct part and no inference part', () => {
        expect(SYSTEM_PROMPT).toBe(systemPrompt('full'));
        expect(SYSTEM_PROMPT.endsWith(`Ground rules: ${groundRules('full')}`)).toBe(true);
        expect(groundRules('full')).toBe(GROUND_RULES.map((parts) => parts.map((part) => part.text).join('')).join(' '));
        const off = groundRules('off');
        expect(systemPrompt('off').endsWith(`Ground rules: ${off}`)).toBe(true);
        expect(systemPrompt('off').slice(0, -off.length)).toBe(SYSTEM_PROMPT.slice(0, -groundRules('full').length));
        for (const part of GROUND_RULES.flat()) {
            if (part.kind === 'inference') expect(off).not.toContain(part.text.trim());
            else expect(off).toContain(part.text);
        }
        // What remains reads as sentences: no doubled space or orphaned punctuation.
        expect(off).not.toMatch(/ {2}|[;:,] ?[.;]|\. \./);
        expect(() => systemPrompt('none')).toThrow(TypeError);
        expect(() => groundRules('none')).toThrow(TypeError);
    });

    it('the tool-free prompt follows the method too: the watchdog\'s and the dock\'s stay whole', () => {
        expect(copilotSystemPrompt({ toolsEnabled: false })).toBe(COPILOT_PROMPT_NO_TOOLS);
        expect(COPILOT_PROMPT_NO_TOOLS).toBe(`${SYSTEM_PROMPT} ${OUTCOME_RULE}`);
        expect(copilotSystemPrompt({ toolsEnabled: false, method: 'off' })).toBe(`${systemPrompt('off')} ${OUTCOME_RULE}`);
        for (const part of GROUND_RULES.flat().filter((item) => item.kind === 'inference')) {
            expect(copilotSystemPrompt({ toolsEnabled: false, method: 'off' })).not.toContain(part.text.trim());
            expect(copilotSystemPrompt({ askUser: true, method: 'off' })).not.toContain(part.text.trim());
        }
    });
});
