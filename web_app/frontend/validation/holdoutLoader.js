// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The one way to the sealed test set (docs/algorithms/copilot-validation.md,
// "The sealed test set"). The held-out cases' definitions live in the sealed
// folder beside this file, which no other module may import or read
// (split.test.js scans src/** and validation/** for it); builders and agents
// never open it. This loader imports it, dynamically, only when asked on
// purpose:
//
//   - on the --final path of run.js (arm-a, live and write with --split test;
//     run.js also requires a clean git tree for arm-a and live), or
//   - when the environment sets RMC_VALIDATION_HOLDOUT=1, as CI's frontend job
//     does, so the tests that need the held-out cases (the test-set hash, the
//     pairing with the dev set, the holdout scan's tokens, the first requests
//     of the case-blindness test, arm D) run there and are skipped locally.
//
// Anything else gets an error that says so. The ids alone are public
// (cases.js TEST_CASE_IDS).

import { casesOf } from './cases';

/** The environment flag that opens the sealed test set in CI. */
export const HOLDOUT_ENV = 'RMC_VALIDATION_HOLDOUT';

const environment = () => globalThis.process?.env ?? {};

/** Whether the flag is set (exactly '1'). */
export const holdoutEnabled = (env = environment()) => env[HOLDOUT_ENV] === '1';

/**
 * The sealed folder, relative to web_app/frontend: the one place its name is
 * written in code (split.test.js checks), so the scans that walk validation/
 * skip it by this name and never read it.
 */
export const HOLDOUT_DIR = 'validation/holdout';

/** Where arms A, A+ and D on the test set are recorded (run.js arm-a --split test --final --out …). */
export const ARM_A_TEST_RECORDS = `${HOLDOUT_DIR}/arm-a-test.json`;

export class HoldoutSealedError extends Error {
    constructor() {
        super('the held-out test set is sealed: it loads only on the --final path of run.js (--split test --final), '
            + `or in CI with ${HOLDOUT_ENV}=1; never open its folder to read it`);
        this.name = 'HoldoutSealedError';
    }
}

/**
 * The held-out test cases (as the sealed module defines them), when `final`
 * is true (run.js's --final) or the flag is set; otherwise a HoldoutSealedError,
 * before anything is imported.
 */
export const loadTestCases = async ({ final = false, env = environment() } = {}) => {
    if (final !== true && !holdoutEnabled(env)) throw new HoldoutSealedError();
    const { TEST_CASES } = await import('./holdout/testCases.js');
    return TEST_CASES;
};

/** The cases of a split: the dev set as cases.js has it, the test set through loadTestCases. */
export const casesFor = async (split, options = {}) => (split === 'test' ? loadTestCases(options) : casesOf(split));
