// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The records the scripted live runs keep (harness.test.js, fitCase.test.js,
// judgment.test.js, src/__tests__/turnEndings.test.jsx), pinned as files in
// __snapshots__/records/: one per run, so a change to what a record holds
// shows as a diff. They were first written before the tool transcript
// (architecture candidate #11), which must leave every value as it was.
//
// A file is the record as runLiveCase returns it, with what changes from one
// run to the next written as a placeholder: the clock (`started_at`,
// `finished_at`, `wall_ms`, each step's `ms`), the protocol when it is the
// current PROTOCOL_VERSION, and each request as sent as its length and hash
// (the model-visible snapshots hold the requests themselves). Key order is
// the record's.
//
// The keys the tool transcript added are left out (ADDED): each step's
// `result_status`, `note`, `choice`, `ref` and `analysis`, the steps a
// stopped turn left running (`status: 'stopped'`), and `safety.refused_calls`.
// So the files show every value the records had before it, unchanged; the
// added keys are pinned where they are tested (toolOutcomes.test.jsx,
// harness.test.js).

import { createHash } from 'node:crypto';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { expect } from 'vitest';
import { PROTOCOL_VERSION } from '../harness';

const HERE = dirname(fileURLToPath(import.meta.url));
export const RECORDS_DIR = join(HERE, '__snapshots__/records');

const digest = (value) => {
    const text = JSON.stringify(value);
    return `${text.length} characters, sha256 ${createHash('sha256').update(text, 'utf8').digest('hex').slice(0, 16)}`;
};

/** The keys the tool transcript added to a record (see the header). */
export const ADDED = Object.freeze({
    stepKeys: Object.freeze(['result_status', 'note', 'choice', 'ref', 'analysis']),
    safetyKeys: Object.freeze(['refused_calls'])
});

/** A record without the keys the tool transcript added: the record as it was before it. */
export const withoutAdded = (record) => {
    const copy = JSON.parse(JSON.stringify(record));
    copy.steps = (copy.steps || []).filter((step) => step.status !== 'stopped')
        .map((step) => Object.fromEntries(Object.entries(step).filter(([key]) => !ADDED.stepKeys.includes(key))));
    if (copy.safety) copy.safety = Object.fromEntries(Object.entries(copy.safety).filter(([key]) => !ADDED.safetyKeys.includes(key)));
    return copy;
};

/** The record as its file holds it (see the header). */
export const recordFile = (record) => {
    const pinned = withoutAdded(record);
    if (pinned.protocol === PROTOCOL_VERSION) pinned.protocol = '<PROTOCOL_VERSION>';
    for (const key of ['started_at', 'finished_at']) if (key in pinned) pinned[key] = '<time>';
    if (Number.isFinite(pinned.wall_ms)) pinned.wall_ms = '<ms>';
    pinned.requests = (pinned.requests || []).map((request) => Object.fromEntries(Object.entries(request)
        .map(([key, value]) => [key, key === 'sent' ? digest(value) : value])));
    pinned.steps = pinned.steps.map((step) => (Number.isFinite(step.ms) ? { ...step, ms: '<ms>' } : step));
    return `${JSON.stringify(pinned, null, 1)}\n`;
};

/** Compare `record` with its file, `<name>.json` (written when missing, outside CI). */
export const pinRecord = (name, record) => expect(recordFile(record)).toMatchFileSnapshot(join(RECORDS_DIR, `${name}.json`));
