// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The holdout scan's two lists (split.test.js reads them). They name the
// cases and files they cover, by id and path only, so they live here and not
// in a test file: no test file writes a held-out id (holdoutGate.test.js), and
// this file holds its own id-only allowance below.

// Two lists, and nothing in either names a test-case value. An entry covers
// exactly `count` hits of one test case's kind of token in one file (and, for
// a leak, on one line): a different count fails, so a new mention there is
// still caught and a removed one must leave the list. An id-only allowance
// (`case: '*'`, kind 'id') covers every held-out id in its file together.
//
// HOLDOUT_ALLOWED: hits the model can never read, in any arm. Model-visible
// text is the system prompt, the tool descriptions and results, read_skill's
// output and every skill file (the SKILL.md files and the reference files
// copilotSkills.js bundles), so an allowed file is none of those, is in no
// snapshot of what the model sees, and no app module can import it (checked
// in split.test.js). copilot-validation.md lists these as known risks.
export const HOLDOUT_ALLOWED = [
    // The four 2026-10-19 entries for the test breathing case's old cluster
    // (the skills' reference compound's) left with protocol 2026-10-19.2: the
    // case was recomposed and nothing names its new cluster.
    //
    // FOUND 2026-10-19, reviewed: a real-compound fixture of the symmetry
    // finder (an orthorhombic Pnma perovskite) whose composition is one test
    // case's, unrelated to its planted physics. Imported only by the symmetry
    // tests in src/__tests__, never by the app, the tools or the snapshot
    // tests.
    {
        file: 'web_app/frontend/src/__tests__/fixtures/symmetryStructures.js', case: 'test-symmetry-lowering', kind: 'composition', count: 1,
        reason: 'test-only: a symmetry-finder fixture imported only by src/__tests__/*.test.js; never model-visible'
    },
    // The ids, where a file must name its cases (the holdout gate, 2026-10-19):
    {
        file: 'docs/algorithms/copilot-validation.md', case: '*', kind: 'id', count: 10,
        reason: 'the protocol\'s own document: it names every case by its id (the Cases table, the protocol notes); no other kind of token'
    },
    {
        file: 'web_app/frontend/validation/cases.js', case: '*', kind: 'id', count: 15,
        reason: 'TEST_CASE_IDS, the public manifest of the held-out ids (findings and records name their case by it); one id is also inside another'
    },
    {
        file: 'web_app/frontend/validation/findings/decisions/2026-10-19-test-breathing-recomposed.json', case: '*', kind: 'id', count: 2,
        reason: 'a decision record on a test case names its case (the findings schema\'s `case`)'
    },
    {
        file: 'web_app/frontend/validation/__tests__/holdoutLists.js', case: '*', kind: 'id', count: 2,
        reason: 'these lists: each entry names the case and the file it covers'
    }
];

// HOLDOUT_LEAKS_PENDING_USER: hits the model CAN read. They are leaks, not
// allowed hits: printed on every run, pinned line by line so a change shows,
// and left for the user, whose text the skills are. The scan still fails on
// any hit in neither list. copilot-validation.md lists them as open leaks.
// An entry is { file, line, case, kind, count, status }, the status saying
// how the text reaches the model and that it awaits the user's decision.
// Empty since protocol 2026-10-19.2: the three lines of
// coordination-shells/lacunar-spinel-shells.md named the test breathing
// case's old cluster, and the case was recomposed instead of the skill.
export const HOLDOUT_LEAKS_PENDING_USER = [];
