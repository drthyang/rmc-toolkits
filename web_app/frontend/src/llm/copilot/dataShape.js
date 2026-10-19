// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The shape of a tool's `data`: the tool spec's `dataShape` (tools.js). With
// numbers by reference the model quotes a value by its path in a result
// ({{r3.sites.0.shift_A}}, references.js), so a result's paths are an
// interface the model writes against, and this is its declaration. It is
// internal: toRequestTools sends a tool's name, description and parameters
// only, so the model never receives it. The tests run every tool on the Demo
// run and on the validation's dev cases and check each result against its
// shape (src/__tests__/toolDataShapes.test.js): a key renamed, dropped or
// added fails there, and a PR that changes what a tool returns changes its
// dataShape in the same diff.
//
// The notation is plain data that reads like the JSON it describes:
//
//   'number'  'string'  'boolean'  'null'   one JSON type ('number' is finite)
//   'number|null'                          any of these types
//   { key: shape }                         an object with exactly the keys declared:
//   { 'key?': shape }                        a key ending in ? may be absent;
//   { 'cn_*_per_*': shape }                  a key holding * stands for any number of
//                                            keys whose names fit it (* = any characters)
//   [shape]                                an array of any length, every item `shape`
//   [shape, shape, …]                      a tuple: exactly these items, in this order
//   oneOf(shape, shape, …)                 any one of these shapes
//
// Values are read as JSON carries them to the model: an object key whose
// value is undefined is absent, and a number that is not finite (or an
// undefined array item) is null. Pure, no dependency.

const TYPES = ['number', 'string', 'boolean', 'null'];
const ONE_OF = Symbol.for('rmc-toolkits.dataShape.oneOf');

/** Any one of `shapes` (two or more). */
export const oneOf = (...shapes) => {
    if (shapes.length < 2) throw new Error('oneOf needs two or more shapes');
    return Object.freeze({ [ONE_OF]: Object.freeze(shapes) });
};

const isOneOf = (shape) => Boolean(shape) && typeof shape === 'object' && ONE_OF in shape;
const isObjectShape = (shape) => Boolean(shape) && typeof shape === 'object' && !Array.isArray(shape) && !isOneOf(shape);

/** A value's type as JSON writes it: 'number' only when finite. */
const jsonType = (value) => {
    if (value === null || value === undefined) return 'null';
    if (typeof value === 'number') return Number.isFinite(value) ? 'number' : 'null';
    if (Array.isArray(value)) return 'array';
    return typeof value;
};

const typeList = (shape) => {
    const types = shape.split('|');
    const unknown = types.filter((type) => !TYPES.includes(type));
    if (unknown.length) throw new Error(`invalid data shape: unknown type ${unknown.map((type) => JSON.stringify(type)).join(', ')} in '${shape}'`);
    return types;
};

const escapeRegExp = (text) => text.replace(/[.+?^${}()|[\]\\]/g, '\\$&');

// An object shape's keys: the named ones (required or optional) and the
// patterns. Compiled once per shape object.
const compiled = new WeakMap();
const objectKeys = (shape) => {
    if (compiled.has(shape)) return compiled.get(shape);
    const named = [];
    const patterns = [];
    Object.keys(shape).forEach((key) => {
        if (key.includes('*')) {
            if (key.endsWith('?')) throw new Error(`invalid data shape: '${key}' — a key holding * is optional already`);
            patterns.push({ key, test: new RegExp(`^${key.split('*').map(escapeRegExp).join('.*')}$`), shape: shape[key] });
        } else {
            const optional = key.endsWith('?');
            const name = optional ? key.slice(0, -1) : key;
            if (!name || named.some((entry) => entry.name === name)) throw new Error(`invalid data shape: key '${name}' is declared twice`);
            named.push({ name, optional, shape: shape[key] });
        }
    });
    const keys = { named, names: new Set(named.map((entry) => entry.name)), patterns };
    compiled.set(shape, keys);
    return keys;
};

/** A shape as one line of text, for messages. */
export const describeShape = (shape) => {
    if (typeof shape === 'string') return shape;
    if (isOneOf(shape)) return `oneOf(${shape[ONE_OF].map(describeShape).join(', ')})`;
    if (Array.isArray(shape)) return `[${shape.map(describeShape).join(', ')}]`;
    if (isObjectShape(shape)) return `{ ${Object.keys(shape).join(', ')} }`;
    return String(shape);
};

const describeValue = (value) => {
    const type = jsonType(value);
    if (type === 'string') return `${JSON.stringify(value.length > 40 ? `${value.slice(0, 40)}…` : value)} (string)`;
    if (type === 'number' || type === 'boolean') return `${value} (${type})`;
    if (type === 'null') return typeof value === 'number' ? `${value} (null in JSON)` : 'null';
    if (type === 'array') return `an array of ${value.length}`;
    return `an ${type}`;
};

// A value's path as a reference names it (sites.0.shift_A); the result itself is '(data)'.
const at = (path) => (path.length ? path.join('.') : '(data)');
const below = (schemaPath, step) => (schemaPath && !step.startsWith('[') ? `${schemaPath}.${step}` : `${schemaPath}${step}`);

const visit = (value, shape, path, schemaPath, problems, seen) => {
    const expected = () => problems.push(`${at(path)}: expected ${describeShape(shape)}, found ${describeValue(value)}`);
    if (typeof shape === 'string') {
        if (!typeList(shape).includes(jsonType(value))) expected();
        return;
    }
    if (isOneOf(shape)) {
        const tries = shape[ONE_OF].map((branch) => {
            const found = [];
            const branchSeen = seen ? new Set() : null;
            visit(value, branch, path, schemaPath, found, branchSeen);
            return { found, branchSeen };
        });
        const fit = tries.find((attempt) => !attempt.found.length);
        if (fit) {
            fit.branchSeen?.forEach((entry) => seen.add(entry));
            return;
        }
        const closest = tries.reduce((best, attempt) => (attempt.found.length < best.found.length ? attempt : best));
        // Wrong at this very value in every form: say so once; else the closest form's problems.
        if (tries.every((attempt) => attempt.found.length === 1 && attempt.found[0].startsWith(`${at(path)}: expected`))) expected();
        else problems.push(...closest.found.map((text) => `${text} (in the closest of ${tries.length} forms)`));
        return;
    }
    if (Array.isArray(shape)) {
        if (!shape.length) throw new Error('invalid data shape: [] declares no item shape');
        if (!Array.isArray(value)) {
            expected();
            return;
        }
        if (shape.length === 1) {
            value.forEach((item, index) => visit(item, shape[0], [...path, index], below(schemaPath, '[]'), problems, seen));
            return;
        }
        if (value.length !== shape.length) {
            problems.push(`${at(path)}: expected a tuple of ${shape.length} ${describeShape(shape)}, found ${describeValue(value)}`);
        }
        shape.forEach((item, index) => {
            if (index < value.length) visit(value[index], item, [...path, index], below(schemaPath, `[${index}]`), problems, seen);
        });
        return;
    }
    if (!isObjectShape(shape)) throw new Error(`invalid data shape: ${describeShape(shape)}`);
    if (jsonType(value) !== 'object') {
        expected();
        return;
    }
    const keys = objectKeys(shape);
    keys.named.forEach(({ name, optional, shape: inner }) => {
        if (!Object.hasOwn(value, name) || value[name] === undefined) {
            if (!optional) problems.push(`${at([...path, name])}: missing (a required key)`);
            return;
        }
        seen?.add(below(schemaPath, name));
        visit(value[name], inner, [...path, name], below(schemaPath, name), problems, seen);
    });
    Object.keys(value).forEach((key) => {
        if (value[key] === undefined || keys.names.has(key)) return;
        const pattern = keys.patterns.find((entry) => entry.test.test(key));
        if (!pattern) {
            problems.push(`${at([...path, key])}: not declared in the shape`);
            return;
        }
        seen?.add(below(schemaPath, pattern.key));
        visit(value[key], pattern.shape, [...path, key], below(schemaPath, pattern.key), problems, seen);
    });
};

/**
 * Every way `value` departs from `shape`, one line each with the value's path
 * as a reference names it ("sites.0.shift_A: expected number|null, found
 * "0.2" (string)"); [] when it fits. With `seen` (a Set), adds the declared
 * key of every value it checked, as declaredKeys() names them. Throws on a
 * shape that is not valid notation.
 */
export const checkShape = (value, shape, { seen = null } = {}) => {
    const problems = [];
    visit(value, shape, [], '', problems, seen);
    return problems;
};

/**
 * Every key a shape declares, with its path in the shape ("sites[].shift_A",
 * "peak.direction_uvw", "first_shell.cn.*_per_*"; a list's items are [], a
 * tuple's [0], [1], …), each once, the keys of every oneOf form included.
 */
export const declaredKeys = (shape, schemaPath = '') => {
    if (typeof shape === 'string') {
        typeList(shape);
        return [];
    }
    if (isOneOf(shape)) return [...new Set(shape[ONE_OF].flatMap((branch) => declaredKeys(branch, schemaPath)))];
    if (Array.isArray(shape)) {
        if (!shape.length) throw new Error('invalid data shape: [] declares no item shape');
        if (shape.length === 1) return declaredKeys(shape[0], below(schemaPath, '[]'));
        return shape.flatMap((item, index) => declaredKeys(item, below(schemaPath, `[${index}]`)));
    }
    if (!isObjectShape(shape)) throw new Error(`invalid data shape: ${describeShape(shape)}`);
    const keys = objectKeys(shape);
    return [
        ...keys.named.flatMap(({ name, shape: inner }) => [below(schemaPath, name), ...declaredKeys(inner, below(schemaPath, name))]),
        ...keys.patterns.flatMap(({ key, shape: inner }) => [below(schemaPath, key), ...declaredKeys(inner, below(schemaPath, key))])
    ];
};
