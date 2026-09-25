// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The RMCProfile `.rmc6f` atom-line grammar — one grammar, shared verbatim with
// `rmc_toolkits/parsers.py` (classify_rmc6f_atom_line / iter_rmc6f_atoms; keep the
// two in sync). An atom line is `id element [label] <data>` where
//
//   id       a non-negative integer (the atom number);
//   element  a token starting with a letter (normalized like Python's
//            str.capitalize(): `SE`/`se` → `Se`);
//   label    optional: a bracket group (`[1]`, or split as `[ 1]`) or one
//            non-numeric token;
//   data     exactly 7 tokens  `x y z ref cx cy cz`   (full layout), or
//            exactly 3 tokens  `x y z`                (legacy coords-only; the
//            site/cell membership is reconstructed downstream by folding and
//            clustering, so referenceNumber / cellIndices come back null).
//
// Numbers accept Fortran `D` exponents. `ref` must be a positive integer and each
// cell index an integer in [0, N_i) (N from the `Supercell` header line, when
// known). A line whose layout is valid but whose coordinates are non-finite
// (NaN, Inf, or a Fortran `****` overflow) is skipped and counted; every other
// line after the Atoms marker that fits no layout is counted as unparsed. Nothing
// is inferred by indexing from the END of the line any more: one extra trailing
// field used to shift every column silently (y,z became x,y; a cell index became
// the reference number) while every value stayed finite.

const NUMBER_RE = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$/;
const NON_FINITE_RE = /^(?:[+-]?(?:nan|inf|infinity)|\*+)$/i;
const INTEGER_RE = /^[+-]?\d+$/;
const ATOMS_MARKER_RE = /^\s*atoms\b/i;
const DECLARED_ATOMS_RE = /^\s*Number of atoms\s*:\s*(\d+)/i;

// Every line-break convention RMCProfile files arrive with (LF, CRLF, bare CR).
export const LINE_BREAK = /\r\n|\r|\n/;

/**
 * A numeric token as a number, Fortran-aware: plain / E / D exponents are
 * numbers, an explicit non-finite token (NaN, Inf, Infinity, `****`) is NaN, and
 * anything else is null. Mirrors parse_fortran_number() in parsers.py.
 */
export const parseFortranNumber = (token) => {
    if (NUMBER_RE.test(token)) return Number(token.replace(/[dD]/, 'e'));
    if (NON_FINITE_RE.test(token)) return NaN;
    return null;
};

/** True for the line that opens the atom list (`Atoms:`, `Atoms :`, `atoms:`, …). */
export const isAtomsMarker = (line) => ATOMS_MARKER_RE.test(line);

const parseInteger = (token) => (INTEGER_RE.test(token) ? Number(token) : null);

// Python str.capitalize(): first character upper-case, the rest lower-case.
const normalizeElement = (token) => token.charAt(0).toUpperCase() + token.slice(1).toLowerCase();

/**
 * Classify one whitespace-split atom line. Returns `{ kind, atom }` with kind
 * 'atom' (full layout), 'coords' (legacy coords-only), 'nonFinite' (valid layout,
 * NaN/Inf coordinate; atom null) or 'invalid'. `supercell` (optional) bounds the
 * cell indices. Mirrors classify_rmc6f_atom_line() in parsers.py.
 */
export const classifyAtomLine = (parts, supercell = null) => {
    const invalid = { kind: 'invalid', atom: null };
    if (parts.length < 5) return invalid;
    const atomNumber = parseInteger(parts[0]);
    if (atomNumber === null || atomNumber < 0 || !/^[A-Za-z]/.test(parts[1])) return invalid;
    const element = normalizeElement(parts[1]);

    let index = 2;
    if (parts[index].startsWith('[')) {
        let closed = false;
        while (index < parts.length) {
            const token = parts[index];
            index += 1;
            if (token.endsWith(']')) { closed = true; break; }
        }
        if (!closed) return invalid;
    } else if (parseFortranNumber(parts[index]) === null) {
        index += 1;
    }
    const data = parts.slice(index);
    if (data.length !== 3 && data.length !== 7) return invalid;

    const coords = data.slice(0, 3).map(parseFortranNumber);
    if (coords.some((value) => value === null)) return invalid;
    let referenceNumber = null;
    let cellIndices = null;
    if (data.length === 7) {
        referenceNumber = parseInteger(data[3]);
        cellIndices = data.slice(4, 7).map(parseInteger);
        if (referenceNumber === null || referenceNumber < 1) return invalid;
        for (let axis = 0; axis < 3; axis += 1) {
            const cell = cellIndices[axis];
            const limit = supercell ? supercell[axis] : null;
            if (cell === null || cell < 0) return invalid;
            if (Number.isFinite(limit) && limit >= 1 && cell >= limit) return invalid;
        }
    }
    if (!coords.every(Number.isFinite)) return { kind: 'nonFinite', atom: null };
    return {
        kind: data.length === 7 ? 'atom' : 'coords',
        atom: { element, coords, referenceNumber, cellIndices }
    };
};

/**
 * One atom line → `{ element, coords, referenceNumber, cellIndices }`, or null for
 * a non-atom / malformed / non-finite line. Coords-only lines return null
 * referenceNumber and cellIndices. `parts` is the whitespace-split, non-empty
 * tokens of the line.
 */
export const parseAtomLine = (parts, supercell = null) => classifyAtomLine(parts, supercell).atom;

/**
 * Lattice vectors (rows, Å) and supercell multiplicities from the header. The
 * error names the file. Mirrors read_cell_vectors() in parsers.py.
 */
export const readRmc6fCellVectors = (text, name = 'structure file') => {
    const lines = text.split(LINE_BREAK);
    let latticeVectors = null;
    let supercell = null;
    lines.forEach((line, index) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Supercell') supercell = parts.slice(-3).map(Number);
        if (parts[0] === 'Lattice') {
            latticeVectors = [lines[index + 1], lines[index + 2], lines[index + 3]]
                .map((row) => (row ?? '').trim().split(/\s+/).map(Number));
        }
    });
    if (!latticeVectors || !supercell) throw new Error(`${name} is missing lattice or supercell metadata`);
    return { latticeVectors, supercell };
};

/**
 * Every atom of an `.rmc6f` text plus a line-by-line report:
 * `{ atoms, report }` where atoms holds both full-layout and coords-only atoms
 * and report is `{ hasAtomsSection, declaredAtoms, atomLines, parsedAtoms,
 * coordsOnlyAtoms, nonFiniteLines, invalidLines, firstInvalidLine,
 * firstNonFiniteLine }` — the same counts as Python's Rmc6fParseReport.
 */
export const parseRmc6fAtoms = (text) => {
    const report = {
        hasAtomsSection: false,
        declaredAtoms: null,
        atomLines: 0,
        parsedAtoms: 0,
        coordsOnlyAtoms: 0,
        nonFiniteLines: 0,
        invalidLines: 0,
        firstInvalidLine: null,
        firstNonFiniteLine: null
    };
    const atoms = [];
    let supercell = null;
    let inAtoms = false;
    for (const line of text.split(LINE_BREAK)) {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) continue;
        if (!inAtoms) {
            if (isAtomsMarker(line)) {
                inAtoms = true;
                report.hasAtomsSection = true;
                continue;
            }
            const declared = DECLARED_ATOMS_RE.exec(line);
            if (declared) {
                report.declaredAtoms = Number(declared[1]);
            } else if (parts[0] === 'Supercell' && parts.length >= 3) {
                const values = parts.slice(-3).map(parseFortranNumber);
                if (values.every((value) => value !== null)) supercell = values;
            }
            continue;
        }
        report.atomLines += 1;
        const { kind, atom } = classifyAtomLine(parts, supercell);
        if (kind === 'atom') {
            report.parsedAtoms += 1;
            atoms.push(atom);
        } else if (kind === 'coords') {
            report.coordsOnlyAtoms += 1;
            atoms.push(atom);
        } else if (kind === 'nonFinite') {
            report.nonFiniteLines += 1;
            if (report.firstNonFiniteLine === null) report.firstNonFiniteLine = line.trim();
        } else {
            report.invalidLines += 1;
            if (report.firstInvalidLine === null) report.firstInvalidLine = line.trim();
        }
    }
    return { atoms, report };
};

/**
 * The human-readable problem list for a parse report, or null when the atom
 * section is clean. Same wording as Rmc6fParseReport.warning() in parsers.py.
 */
export const rmc6fParseWarning = (report) => {
    const problems = [];
    const accepted = report.parsedAtoms + report.coordsOnlyAtoms;
    if (report.declaredAtoms !== null && accepted !== report.declaredAtoms) {
        problems.push(`parsed ${accepted} of ${report.declaredAtoms} atoms declared in the header`);
    }
    if (report.invalidLines) {
        problems.push(`${report.invalidLines} of ${report.atomLines} atom lines unparsed (first: '${report.firstInvalidLine}')`);
    }
    if (report.nonFiniteLines) {
        problems.push(`${report.nonFiniteLines} atom lines skipped for non-finite coordinates (first: '${report.firstNonFiniteLine}')`);
    }
    return problems.length ? problems.join('; ') : null;
};
