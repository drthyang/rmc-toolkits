// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The sealed test set's gate (docs/algorithms/copilot-validation.md, "The
// sealed test set"):
//
//   a. no module but validation/holdoutLoader.js names the sealed folder in
//      code: no static or dynamic import, require, glob or file path into it,
//      anywhere in src/** or validation/** (comments may describe it);
//   b. the default suite never collects a held-out case by name: no test file
//      and no snapshot of what the model sees writes a held-out id, and no
//      test title is built from the held-out list (with the flag, the sealed
//      ids are checked to be exactly the public ones, so the scan is whole);
//   c. the loader refuses without --final or the flag, and the CLI without
//      --final whatever the environment.
//
// This file reads every module's code but never the sealed folder: it skips
// that folder by the loader's own name for it (HOLDOUT_DIR).

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { casesOf, DEV_CASES, TEST_CASE_IDS } from '../cases';
import { casesFor, HOLDOUT_DIR, HOLDOUT_ENV, holdoutEnabled, HoldoutSealedError, loadTestCases } from '../holdoutLoader';
import { selectCases } from '../run';

const FRONTEND = resolve(__dirname, '../..');
const SEALED = join(FRONTEND, HOLDOUT_DIR);
const LOADER = join(FRONTEND, 'validation/holdoutLoader.js');
const HOLDOUT = holdoutEnabled();

const filesUnder = (dir, keep) => readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    if (name === 'node_modules' || name === 'dist' || path === SEALED) return [];
    return statSync(path).isDirectory() ? filesUnder(path, keep) : (keep(path) ? [path] : []);
});
const CODE = /\.(js|jsx|mjs|cjs|ts|tsx)$/;
const modules = () => [
    ...filesUnder(join(FRONTEND, 'src'), (path) => CODE.test(path)),
    ...filesUnder(join(FRONTEND, 'validation'), (path) => CODE.test(path))
];

// --- a small lexer: a module's string and template literals, comments skipped ---

// A `/` starts a regular expression after one of these (or at the start), a
// division otherwise.
const REGEX_AFTER = /[(,=:[!&|?{};+\-*%<>~^]$|(?:^|[^\w$])(?:return|typeof|case|in|of|new|delete|void|throw|else|do|yield|await)$/;

/**
 * `{ literals: [{ text, line }], ok }`: every quoted string and template
 * chunk, and whether the lexer ended cleanly. In JSX (`jsx: true`) a quote
 * whose string does not close on its line is text (an apostrophe in JSX
 * text), not a string.
 */
export const stringLiterals = (code, { jsx = false } = {}) => {
    const literals = [];
    let i = 0;
    let line = 1;
    let last = '';
    const templateDepth = [];
    const take = (quote) => {
        const start = line;
        let text = '';
        i += 1;
        while (i < code.length) {
            const char = code[i];
            if (char === '\\') {
                text += code.slice(i, i + 2);
                if (code[i + 1] === '\n') line += 1;
                i += 2;
                continue;
            }
            if (char === '\n') line += 1;
            if (quote === '`' && char === '$' && code[i + 1] === '{') {
                literals.push({ text, line: start });
                templateDepth.push(0);
                i += 2;
                return 'open';
            }
            if (char === quote) {
                literals.push({ text, line: start });
                i += 1;
                return 'closed';
            }
            if (quote !== '`' && char === '\n') return 'broken';
            text += char;
            i += 1;
        }
        return 'eof';
    };
    while (i < code.length) {
        const char = code[i];
        const next = code[i + 1];
        if (char === '\n') { line += 1; i += 1; continue; }
        if (/\s/.test(char)) { i += 1; continue; }
        if (char === '/' && next === '/') { while (i < code.length && code[i] !== '\n') i += 1; continue; }
        if (char === '/' && next === '*') {
            const end = code.indexOf('*/', i + 2);
            const stop = end < 0 ? code.length : end + 2;
            line += (code.slice(i, stop).match(/\n/g) || []).length;
            i = stop;
            continue;
        }
        if (char === '\'' || char === '"' || char === '`') {
            const [from, fromLine] = [i, line];
            const state = take(char);
            if (jsx && char !== '`' && state === 'broken') {
                [i, line] = [from + 1, fromLine];
                last = char;
                continue;
            }
            if (state === 'broken' || state === 'eof') return { literals, ok: false, line };
            last = state === 'open' ? '{' : 'x';
            continue;
        }
        if (char === '}' && templateDepth.length && templateDepth[templateDepth.length - 1] === 0) {
            templateDepth.pop();
            const state = take('`');
            if (state === 'broken' || state === 'eof') return { literals, ok: false, line };
            last = state === 'open' ? '{' : 'x';
            continue;
        }
        if (char === '{' && templateDepth.length) templateDepth[templateDepth.length - 1] += 1;
        if (char === '}' && templateDepth.length) templateDepth[templateDepth.length - 1] -= 1;
        if (char === '/' && REGEX_AFTER.test(last || '(') && !(jsx && last === '<')) {
            // A regular expression literal: to its closing slash, past classes and escapes.
            const from = i;
            let inClass = false;
            i += 1;
            while (i < code.length && code[i] !== '\n') {
                if (code[i] === '\\') { i += 2; continue; }
                if (code[i] === '[') inClass = true;
                else if (code[i] === ']') inClass = false;
                else if (code[i] === '/' && !inClass) break;
                i += 1;
            }
            if (code[i] !== '/') {
                // In JSX a slash that opens no expression on its line is text (or a tag's).
                if (!jsx) return { literals, ok: false, line };
                i = from + 1;
                last = '/';
                continue;
            }
            i += 1;
            while (/[a-z]/i.test(code[i] || '')) i += 1;
            last = 'x';
            continue;
        }
        // Any other token: remember its tail, for the regex rule.
        const word = /^[\w$]+/.exec(code.slice(i));
        if (word) {
            last = word[0];
            i += word[0].length;
        } else {
            last = char;
            i += 1;
        }
    }
    return { literals, ok: templateDepth.length === 0, line };
};

// A literal that names the sealed folder: "holdout" as a path segment.
const SEALED_NAME = HOLDOUT_DIR.split('/').pop();
const namesSealed = (text) => text === SEALED_NAME || new RegExp(`(^|[\\\\/])${SEALED_NAME}([\\\\/]|$)`).test(text);

describe('a. only the loader reaches the sealed folder', () => {
    const files = modules();
    const lexed = files.map((path) => ({ file: relative(FRONTEND, path), path, ...stringLiterals(readFileSync(path, 'utf8'), { jsx: /x$/.test(path) }) }));

    it('reads every module of src/** and validation/** to its end (the lexer stays in step), never the sealed folder', () => {
        expect(files.length).toBeGreaterThan(300);
        expect(files.some((path) => path.startsWith(`${SEALED}/`))).toBe(false);
        expect(lexed.filter((entry) => !entry.ok).map((entry) => `${entry.file}:${entry.line}`)).toEqual([]);
    });

    it('no module but the loader names it: no import, dynamic import, require, glob or path', () => {
        const named = lexed.filter((entry) => entry.path !== LOADER)
            .flatMap((entry) => entry.literals.filter((literal) => namesSealed(literal.text)).map((literal) => `${entry.file}:${literal.line}`));
        expect(named).toEqual([]);
        // The loader does: its name for the folder and its dynamic import (the check sees them).
        const loader = lexed.find((entry) => entry.path === LOADER);
        expect(loader.literals.filter((literal) => namesSealed(literal.text)).map((literal) => literal.text)).toEqual([HOLDOUT_DIR, `./${SEALED_NAME}/testCases.js`]);
        expect(readFileSync(LOADER, 'utf8')).toMatch(new RegExp(`await import\\('\\./${SEALED_NAME}/testCases\\.js'\\)`));
    });

    it('no glob can reach it: no import.meta.glob in validation/**, and none in src/** into validation/', () => {
        const globs = lexed.flatMap((entry) => {
            const code = readFileSync(entry.path, 'utf8');
            const calls = [...code.matchAll(/import\.meta\.glob(?:Eager)?\s*\(\s*(\[[^\]]*\]|['"`][^'"`]*['"`])/g)].map((match) => match[1]);
            return calls.filter((pattern) => entry.file.startsWith('validation/') || /validation/.test(pattern)).map((pattern) => `${entry.file}: ${pattern}`);
        });
        expect(globs).toEqual([]);
    });

    it('the lexer: strings and templates, not comments or regular expressions', () => {
        const code = [
            '// a comment with \'quotes\' and a path ./x/y',
            'const a = \'one\'; const b = "two"; /* "not" */',
            'const r = /[\'"`]/g; const half = 6 / 2 / 3;',
            'const t = `three ${a + `four`} five`;',
            'import(\'./six.js\'); require("seven");'
        ].join('\n');
        const { literals, ok } = stringLiterals(code);
        expect(ok).toBe(true);
        expect(literals.map((literal) => literal.text)).toEqual(['one', 'two', 'three ', 'four', ' five', './six.js', 'seven']);
        expect(stringLiterals('const broken = \'no end\n').ok).toBe(false);
        // JSX text: an apostrophe is text; an attribute is a string.
        const jsx = stringLiterals('const p = <p title="t">the case\'s data</p>;\nconst q = \'r\';', { jsx: true });
        expect(jsx).toMatchObject({ ok: true });
        expect(jsx.literals.map((literal) => literal.text)).toEqual(['t', 'r']);
        expect(namesSealed(`../${SEALED_NAME}/x.js`)).toBe(true);
        expect(namesSealed(`${SEALED_NAME}`)).toBe(true);
        expect(namesSealed(`${SEALED_NAME}-copy/x`)).toBe(false);
        expect(namesSealed(`the ${SEALED_NAME} set`)).toBe(false);
    });
});

describe('b. the default suite collects no held-out case by name', () => {
    const SNAPSHOTS = join(FRONTEND, 'src/__tests__/modelVisible/__snapshots__');
    const testFiles = [...filesUnder(join(FRONTEND, 'src'), (path) => /\.test\.(js|jsx)$/.test(path)), ...filesUnder(join(FRONTEND, 'validation'), (path) => /\.test\.(js|jsx)$/.test(path))];

    it('no test file and no snapshot writes a held-out id', () => {
        expect(testFiles.length).toBeGreaterThan(100);
        const snapshots = filesUnder(SNAPSHOTS, () => true);
        expect(snapshots.length).toBeGreaterThan(10);
        const found = [...testFiles, ...snapshots].flatMap((path) => {
            const text = readFileSync(path, 'utf8');
            return TEST_CASE_IDS.filter((id) => text.includes(id)).map(() => relative(FRONTEND, path));
        });
        expect(found, 'a test file or snapshot names a held-out case (use TEST_CASE_IDS or the loader)').toEqual([]);
    });

    it('no test title is built from the held-out list', () => {
        const built = testFiles.filter((path) => {
            const code = readFileSync(path, 'utf8');
            return /\.(each|for)\s*\(\s*[^)]*\b(TEST_CASE_IDS|TEST_CASES)\b/.test(code)
                || /\b(it|test|describe)(\.\w+)*(\([^)]*\))?\s*\(\s*`[^`]*\$\{[^}]*\b(TEST_CASE_IDS|TEST_CASES)\b/.test(code);
        }).map((path) => relative(FRONTEND, path));
        expect(built).toEqual([]);
    });

    it.skipIf(!HOLDOUT)('the sealed ids are exactly the public ones, so the scan above is whole', async () => {
        expect((await loadTestCases()).map((item) => item.id)).toEqual([...TEST_CASE_IDS]);
    });
});

describe('c. the loader refuses without --final or the flag', () => {
    it('loadTestCases: only final === true or the flag set to exactly 1', async () => {
        await expect(loadTestCases({ env: {} })).rejects.toBeInstanceOf(HoldoutSealedError);
        await expect(loadTestCases({ final: false, env: { [HOLDOUT_ENV]: '0' } })).rejects.toThrow(/sealed/);
        await expect(loadTestCases({ final: 'true', env: {} })).rejects.toThrow(/--final/);
        await expect(loadTestCases({ env: { [HOLDOUT_ENV]: 'yes' } })).rejects.toThrow(HOLDOUT_ENV);
        expect(holdoutEnabled({})).toBe(false);
        expect(holdoutEnabled({ [HOLDOUT_ENV]: 'true' })).toBe(false);
        expect(holdoutEnabled({ [HOLDOUT_ENV]: '1' })).toBe(true);
    });

    it('the split helpers: dev as before, test only through the loader', async () => {
        expect(casesOf('dev')).toBe(DEV_CASES);
        expect(() => casesOf('test')).toThrow(/sealed/);
        expect(() => casesOf('nope')).toThrow(/unknown split/);
        expect(await casesFor('dev', { env: {} })).toBe(DEV_CASES);
        await expect(casesFor('test', { env: {} })).rejects.toBeInstanceOf(HoldoutSealedError);
    });

    it('the CLI: the test set only with --final, whatever the environment', async () => {
        await expect(selectCases('test', null, { final: false })).rejects.toBeInstanceOf(HoldoutSealedError);
        await expect(selectCases('test', null, {})).rejects.toBeInstanceOf(HoldoutSealedError);
        expect(await selectCases('dev', 'dev-offcentre', {})).toEqual([DEV_CASES[0]]);
    });

    it.skipIf(!HOLDOUT)('opens with the flag, and with final === true (CI)', async () => {
        expect(await loadTestCases()).toHaveLength(TEST_CASE_IDS.length);
        expect(await loadTestCases({ final: true, env: {} })).toHaveLength(TEST_CASE_IDS.length);
        expect(await selectCases('test', null, { final: true })).toHaveLength(TEST_CASE_IDS.length);
    });
});
