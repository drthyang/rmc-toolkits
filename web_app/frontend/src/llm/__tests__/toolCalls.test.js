// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    accumulateToolCall,
    canonicalJson,
    finalizeToolCalls,
    parseToolArguments,
    validateToolArguments
} from '../copilot/toolCalls';

const assemble = (fragments) => finalizeToolCalls(fragments.reduce(accumulateToolCall, []), 't');

const fragment = (fields) => ({ index: null, id: null, name: null, arguments: null, ...fields });

describe('accumulateToolCall', () => {
    it('joins OpenAI-style fragments by index, arguments split mid-token', () => {
        expect(assemble([
            fragment({ index: 0, id: 'call_a', name: 'list_sites', arguments: '' }),
            fragment({ index: 0, arguments: '{"sec' }),
            fragment({ index: 0, arguments: 'tion":"dat' }),
            fragment({ index: 0, arguments: 'asets"}' })
        ])).toEqual([{ id: 'call_a', name: 'list_sites', arguments: '{"section":"datasets"}', extraContent: null }]);
    });

    it('keeps parallel calls apart by index, in arrival order', () => {
        const calls = assemble([
            fragment({ index: 0, id: 'a', name: 'one', arguments: '{"x":' }),
            fragment({ index: 1, id: 'b', name: 'two', arguments: '{"y":' }),
            fragment({ index: 0, arguments: '1}' }),
            fragment({ index: 1, arguments: '2}' })
        ]);
        expect(calls).toEqual([
            { id: 'a', name: 'one', arguments: '{"x":1}', extraContent: null },
            { id: 'b', name: 'two', arguments: '{"y":2}', extraContent: null }
        ]);
    });

    it('reads whole calls without an index (Gemini compat) as separate calls', () => {
        const calls = assemble([
            fragment({ id: 'g1', name: 'one', arguments: '{"a":1}', extraContent: null }),
            fragment({ id: 'g2', name: 'one', arguments: '{"a":2}', extraContent: null })
        ]);
        expect(calls.map((call) => call.arguments)).toEqual(['{"a":1}', '{"a":2}']);
    });

    it('splits a reused index when a new id names a function', () => {
        const calls = assemble([
            fragment({ index: 0, id: 'x1', name: 'one', arguments: '{}' }),
            fragment({ index: 0, id: 'x2', name: 'two', arguments: '{}' })
        ]);
        expect(calls.map((call) => call.name)).toEqual(['one', 'two']);
    });

    it('splits id-less calls that name different functions, and continues bare fragments', () => {
        const calls = assemble([
            fragment({ name: 'one', arguments: '{"a"' }),
            fragment({ arguments: ':1}' }),
            fragment({ name: 'two', arguments: '{}' })
        ]);
        expect(calls).toEqual([
            { id: 'call_t0', name: 'one', arguments: '{"a":1}', extraContent: null },
            { id: 'call_t1', name: 'two', arguments: '{}', extraContent: null }
        ]);
    });

    it('keeps a call\'s provider extra_content (Gemini thought_signature)', () => {
        const signature = { google: { thought_signature: 'sig==' } };
        const calls = assemble([
            fragment({ id: 'g1', name: 'list_sites', arguments: '{}', extraContent: signature })
        ]);
        expect(calls).toEqual([{ id: 'g1', name: 'list_sites', arguments: '{}', extraContent: signature }]);
    });

    it('drops calls that never got a function name', () => {
        expect(assemble([fragment({ index: 0, arguments: '{}' })])).toEqual([]);
    });
});

describe('parseToolArguments', () => {
    it('treats empty arguments as an empty object', () => {
        expect(parseToolArguments('')).toEqual({ ok: true, value: {} });
        expect(parseToolArguments('  ')).toEqual({ ok: true, value: {} });
        expect(parseToolArguments(null)).toEqual({ ok: true, value: {} });
    });

    it('accepts a fenced JSON object', () => {
        expect(parseToolArguments('```json\n{"a": 1}\n```')).toEqual({ ok: true, value: { a: 1 } });
    });

    it('rejects invalid JSON and non-objects with a readable error', () => {
        expect(parseToolArguments('{a:1}').ok).toBe(false);
        expect(parseToolArguments('{a:1}').error).toMatch(/not valid JSON/);
        expect(parseToolArguments('[1,2]')).toEqual({ ok: false, error: 'arguments must be a JSON object' });
        expect(parseToolArguments('3')).toEqual({ ok: false, error: 'arguments must be a JSON object' });
    });
});

describe('validateToolArguments', () => {
    const schema = {
        type: 'object',
        properties: {
            pair: { type: 'string', enum: ['Ta-Se', 'Ga-Se'] },
            r_max: { type: 'number', minimum: 1, maximum: 15 },
            n: { type: 'integer', minimum: 1 },
            flag: { type: 'boolean' }
        },
        required: ['pair']
    };

    it('passes valid arguments and coerces numeric strings', () => {
        expect(validateToolArguments(schema, { pair: 'Ta-Se', r_max: '2.5', n: 3, flag: 'true' })).toEqual({
            ok: true, value: { pair: 'Ta-Se', r_max: 2.5, n: 3, flag: true }, errors: [], notes: []
        });
    });

    it('reports every violation the model can fix', () => {
        const result = validateToolArguments(schema, { pair: 'Nb-Se', r_max: 40, n: 1.5, flag: 'yes' });
        expect(result.ok).toBe(false);
        expect(result.errors).toEqual([
            'pair must be one of "Ta-Se", "Ga-Se", got "Nb-Se"',
            'r_max must be ≤ 15, got 40',
            'n must be an integer, got 1.5',
            'flag must be true or false, got "yes"'
        ]);
    });

    it('requires required keys, treats null as absent, and drops unknown keys with a note', () => {
        const result = validateToolArguments(schema, { r_max: null, extra: 1 });
        expect(result.errors).toEqual(['pair is required']);
        expect(result.notes).toEqual(['ignored unknown argument "extra"']);
        expect(result.value).toEqual({});
    });

    it('rejects non-finite numbers', () => {
        expect(validateToolArguments(schema, { pair: 'Ta-Se', r_max: 'abc' }).errors[0]).toMatch(/finite number/);
    });
});

describe('canonicalJson', () => {
    it('is independent of key order', () => {
        expect(canonicalJson({ b: 1, a: { d: 2, c: [1, { f: 1, e: 2 }] } }))
            .toBe(canonicalJson({ a: { c: [1, { e: 2, f: 1 }], d: 2 }, b: 1 }));
    });
});
