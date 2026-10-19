// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Beside the case score (plan item 1.6, part i): what the report adds so that
// a comparison with A+ says something, with the case score unchanged.
//
// The engines alone (arm A+) score 0.99, because most of a case's components
// are facts the engines get right every time: four symmetry fields and the
// spurious-claim test beside one planted feature. So, per split:
// - the case score's parts: each column's score restricted to the planted
//   features (the feature-only score), to the spurious-claim test and to the
//   symmetry fields (scoring.js scoreComponents, restrictedScore), and each
//   live column minus A+ on the cases both have;
// - every component, each live column against A+: its value, A+'s and the
//   difference, with the other model-free columns (arm A) beside them;
// - every disagreement with the planted truth: per case, each feature or
//   fact a run got wrong (missed, read as another kind or element; a
//   direction, magnitude or q off; a spurious claim; a symmetry field wrong
//   or unanswered; no verdict), with the columns and runs that made it.
// It reads run records only: nothing here changes a record or a score.

import { restrictedScore, scoreComponents, scoreVerdict, TOLERANCES } from './scoring';

/** A record's report column; a run given the PCA table (run.js live --pca-table, a dev check) is a column of its own. */
export const columnOf = (record) => `${record.model?.name ?? 'no model'} · ${record.arm}${record.view?.pcaSites ? ' + PCA table' : ''}`;

/** The arm every live column is set against: the engines alone, no model. */
export const REFERENCE_ARM = 'A+';

/** The groups of the case score's components, as the parts table names them. */
export const COMPONENT_GROUPS = Object.freeze([['feature', 'Feature-only score'], ['spurious', 'No spurious claim'], ['symmetry', 'Symmetry fields']]);

const mean = (values) => (values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null);
const fmt = (value) => (value == null ? '–' : value.toFixed(2));
const difference = (value, reference) => (value == null || reference == null ? null : value - reference);
const isZero = (value) => Math.abs(value) < 0.005;

/** A difference as the tables print it: +0.50, −1.00, 0.00, or – when either side is missing. */
export const fmtDifference = (value) => {
    if (value == null) return '–';
    if (isZero(value)) return '0.00';
    return `${value > 0 ? '+' : '−'}${Math.abs(value).toFixed(2)}`;
};

const num = (value) => String(Number(value.toPrecision(3)));
const uvwText = (uvw) => `[${uvw.map((x) => Number(Number(x).toFixed(2))).join(' ')}]`;
const qText = (q) => `(${q.map((x) => Number(Number(x).toFixed(3))).join(', ')})`;
const sameElement = (a, b) => typeof a === 'string' && typeof b === 'string' && a.trim().toLowerCase() === b.trim().toLowerCase();
const nameOf = (item) => [item.kind, item.element].filter(Boolean).join(' ');
const magnitudeTolerance = (planted) => Math.max(TOLERANCES.MAGNITUDE_REL * planted, TOLERANCES.MAGNITUDE_ABS);

// A planted feature's quantities, as the truth column prints them.
const plantedText = (feature) => [
    feature.q ? `q ${qText(feature.q)}` : null,
    feature.direction_uvw ? `${uvwText(feature.direction_uvw)}${feature.direction_mode === 'cubic' ? ' or any cubic image' : ' or its negative'}` : null,
    Number.isFinite(feature.magnitude_A) ? `${num(feature.magnitude_A)} Å` : null
].filter(Boolean).join(', ');

// A claim's quantities, as the answer column prints them.
const claimText = (claim) => [
    claim.q ? `q ${qText(claim.q)}` : null,
    claim.direction_uvw ? uvwText(claim.direction_uvw) : null,
    Number.isFinite(claim.magnitude_A) ? `${num(claim.magnitude_A)} Å` : null
].filter(Boolean).join(', ');

const SYMMETRY_LABELS = { space_group: 'space group', chiral: 'chiral', centrosymmetric: 'centrosymmetric', polar: 'polar' };
const SYMMETRY_ORDER = Object.keys(SYMMETRY_LABELS);
// The order disagreements are listed in, within a case.
const RANK = { verdict: 0, missed: 1, direction: 2, magnitude: 3, q: 4, spurious: 5, symmetry: 6 };

/**
 * Every disagreement of one run with its planted truth: `{ key, rank, where,
 * answer, truth }`, one per component the run lost, so a component below 1
 * has at least one, and each one costs a component: a missed feature (with
 * the spurious claim that took its place, of another kind on its element or
 * of its kind on another element), each quantity of a found feature that is
 * off (direction, magnitude, q), each spurious claim, each symmetry field
 * wrong or unanswered. A run with no parseable verdict has one: that.
 * Tolerated claims (true partial descriptions) are not disagreements.
 */
export const disagreementsOf = (record) => {
    const truth = record.truth ?? null;
    const score = record.score ?? (truth ? scoreVerdict(truth, null) : null);
    if (!score) return [];
    if (score.format_ok === false) {
        return [{
            key: 'verdict', rank: [RANK.verdict, 0], where: 'no parseable verdict',
            answer: record.verdict_error || (record.error ? `${record.error.name}: no answer` : 'none'), truth: 'a VERDICT line'
        }];
    }
    const out = [];
    const spurious = score.spurious || [];
    (score.features || []).forEach((feature, index) => {
        const planted = truth?.features?.[index] ?? { kind: feature.planted ?? feature.kind, element: feature.element };
        const name = nameOf(planted);
        const rank = (type) => [RANK[type], index];
        if (!feature.found) {
            const accept = planted.accept ?? [planted.kind];
            const instead = spurious.flatMap((claim) => {
                const element = sameElement(claim.element, planted.element);
                const kind = accept.includes(claim.kind);
                if (element && !kind) return [`${nameOf(claim)} (wrong kind)`];
                if (kind && !element) return [`${nameOf(claim)} (wrong element)`];
                return [];
            });
            out.push({
                key: `missed|${name}`, rank: rank('missed'), where: `missed: ${name}`,
                answer: instead.length ? `read as ${instead.join(', ')}` : 'not claimed', truth: `planted${plantedText(planted) ? `: ${plantedText(planted)}` : ''}`
            });
            return;
        }
        const claim = record.verdict?.features?.[feature.claim_index] ?? {};
        if (feature.direction_ok === false) {
            out.push({
                key: `direction|${name}`, rank: rank('direction'), where: `direction of ${name}`,
                answer: feature.direction_error_deg == null ? 'not given' : `${claim.direction_uvw ? `${uvwText(claim.direction_uvw)}, ` : ''}${Math.round(feature.direction_error_deg)}° off`,
                truth: planted.direction_uvw
                    ? `${uvwText(planted.direction_uvw)}${planted.direction_mode === 'cubic' ? ' or any cubic image' : ' or its negative'}, within ${TOLERANCES.DIRECTION_DEG}°`
                    : `within ${TOLERANCES.DIRECTION_DEG}°`
            });
        }
        if (feature.magnitude_ok === false) {
            out.push({
                key: `magnitude|${name}`, rank: rank('magnitude'), where: `magnitude of ${name}`,
                answer: Number.isFinite(feature.magnitude_claimed_A) ? `${num(feature.magnitude_claimed_A)} Å` : 'not given',
                truth: Number.isFinite(planted.magnitude_A) ? `${num(planted.magnitude_A)} ± ${num(magnitudeTolerance(planted.magnitude_A))} Å` : 'planted magnitude'
            });
        }
        if (feature.q_ok === false) {
            out.push({
                key: `q|${name}`, rank: rank('q'), where: `q of ${name}`,
                answer: feature.q_claimed ? qText(feature.q_claimed) : 'not given',
                truth: planted.q ? `${qText(planted.q)} (≡ −q, q + G), within ${TOLERANCES.Q_ABS}` : 'planted q'
            });
        }
    });
    const planted = (truth?.features || []).map(nameOf);
    spurious.forEach((claim, index) => out.push({
        key: `spurious|${nameOf(claim)}`, rank: [RANK.spurious, index], where: `spurious claim: ${nameOf(claim)}`,
        answer: claimText(claim) || 'claimed', truth: planted.length ? `not planted (planted: ${planted.join(', ')})` : 'no order planted'
    }));
    Object.entries(score.symmetry || {}).forEach(([field, result]) => {
        if (result === 'right') return;
        const claimed = record.verdict?.[field];
        out.push({
            key: `symmetry|${field}`, rank: [RANK.symmetry, SYMMETRY_ORDER.indexOf(field)], where: SYMMETRY_LABELS[field] ?? field,
            answer: result === 'unanswered' || claimed == null ? 'not answered' : String(claimed),
            truth: truth?.symmetry?.[field] == null ? '–' : String(truth.symmetry[field])
        });
    });
    return out;
};

// One record as the comparison reads it. A record without a score counts as
// no verdict (score 0), as in the report's own tables; one with neither a
// score nor a truth cannot be taken apart and is left out.
const entryOf = (record) => {
    const score = record.score ?? (record.truth ? scoreVerdict(record.truth, null) : null);
    if (!score) return null;
    return {
        record, score, column: columnOf(record), arm: record.arm, live: Boolean(record.model), case: record.case,
        components: scoreComponents(score, record.truth)
    };
};

// The runs of one column on one case: each component's mean over the runs,
// each group's restricted score and the case score, means over the runs.
const cellOf = (entries) => {
    const components = new Map();
    entries.forEach((entry) => entry.components.forEach((component) => {
        if (!components.has(component.id)) components.set(component.id, { id: component.id, group: component.group, values: [] });
        components.get(component.id).values.push(component.value);
    }));
    const groups = Object.fromEntries(COMPONENT_GROUPS.map(([group]) => {
        const values = entries.map((entry) => restrictedScore(entry.components, group)).filter((value) => value != null);
        return [group, mean(values)];
    }));
    return {
        runs: entries.length,
        score: mean(entries.map((entry) => entry.score.score ?? 0)),
        groups,
        components: new Map([...components].map(([id, item]) => [id, { group: item.group, value: mean(item.values) }]))
    };
};

// A column's parts over a set of its cases, each case weighing once.
const partsOver = (cells, cases) => {
    const mine = cases.map((id) => cells.get(id)).filter(Boolean);
    return {
        cases: mine.length,
        score: mean(mine.map((cell) => cell.score)),
        groups: Object.fromEntries(COMPONENT_GROUPS.map(([group]) => {
            const values = mine.map((cell) => cell.groups[group]).filter((value) => value != null);
            return [group, { value: mean(values), cases: values.length }];
        }))
    };
};

/**
 * One split's records taken apart (any split; the report calls it per split):
 * - `columns` (sorted), `reference` (the A+ column, or null), `live` (the
 *   columns with a model), `shown` (the component table's columns: the other
 *   model-free columns but the truth's, then A+, then the live ones);
 * - `parts`: per column, its case score and each group's restricted score
 *   (the feature-only score among them), each case weighing once;
 * - `differences`: per live column, its parts minus A+'s on the cases both have;
 * - `components`: per case, every component with each shown column's value and
 *   each live column's difference from A+ (null where either is missing);
 * - `disagreements`: per case, every disagreement any run made, with the
 *   columns that made it (`made`: runs that did, of the column's runs on the
 *   case) and the answers they gave; `clean`: the cases no run disagrees on.
 */
export const compareWithReference = (records) => {
    const entries = records.map(entryOf).filter(Boolean);
    const columns = [...new Set(entries.map((entry) => entry.column))].sort();
    const grouped = new Map(columns.map((column) => [column, new Map()]));
    entries.forEach((entry) => {
        const byCase = grouped.get(entry.column);
        if (!byCase.has(entry.case)) byCase.set(entry.case, []);
        byCase.get(entry.case).push(entry);
    });
    const cells = new Map(columns.map((column) => [column, new Map([...grouped.get(column)].map(([id, list]) => [id, cellOf(list)]))]));
    const armOf = (column) => entries.find((entry) => entry.column === column);
    const reference = columns.find((column) => !armOf(column).live && armOf(column).arm === REFERENCE_ARM) ?? null;
    const live = columns.filter((column) => armOf(column).live);
    const context = columns.filter((column) => !armOf(column).live && column !== reference && armOf(column).arm !== 'D');
    const shown = [...context, ...(reference ? [reference] : []), ...live];
    const cases = [...new Set(entries.map((entry) => entry.case))].sort();

    const parts = columns.map((column) => ({ column, live: live.includes(column), ...partsOver(cells.get(column), [...cells.get(column).keys()]) }));
    const differences = reference ? live.map((column) => {
        const paired = [...cells.get(column).keys()].filter((id) => cells.get(reference).has(id));
        const mine = partsOver(cells.get(column), paired);
        const theirs = partsOver(cells.get(reference), paired);
        return {
            column, cases: paired.length, score: difference(mine.score, theirs.score),
            groups: Object.fromEntries(COMPONENT_GROUPS.map(([group]) => [group, difference(mine.groups[group].value, theirs.groups[group].value)]))
        };
    }) : [];

    const components = cases.flatMap((id) => {
        const order = [];
        [reference, ...shown].filter(Boolean).forEach((column) => {
            const cell = cells.get(column).get(id);
            cell?.components.forEach((_, componentId) => { if (!order.includes(componentId)) order.push(componentId); });
        });
        return order.map((componentId) => {
            const valueOf = (column) => cells.get(column).get(id)?.components.get(componentId)?.value ?? null;
            const groupOf = shown.map((column) => cells.get(column).get(id)?.components.get(componentId)?.group).find(Boolean);
            return {
                case: id, id: componentId, group: groupOf,
                values: Object.fromEntries(shown.map((column) => [column, valueOf(column)])),
                differences: reference ? Object.fromEntries(live.map((column) => [column, difference(valueOf(column), valueOf(reference))])) : {}
            };
        });
    });

    const found = new Map();
    entries.forEach((entry) => disagreementsOf(entry.record).forEach((item) => {
        const key = `${entry.case}|${item.key}`;
        if (!found.has(key)) found.set(key, { case: entry.case, key: item.key, rank: item.rank, where: item.where, truth: item.truth, answers: [], made: new Map() });
        const group = found.get(key);
        if (!group.answers.includes(item.answer)) group.answers.push(item.answer);
        group.made.set(entry.column, (group.made.get(entry.column) || 0) + 1);
    }));
    const disagreements = [...found.values()]
        .sort((a, b) => a.case.localeCompare(b.case) || a.rank[0] - b.rank[0] || a.rank[1] - b.rank[1] || a.key.localeCompare(b.key))
        .map(({ made, ...item }) => ({
            ...item,
            made: [...made].sort(([a], [b]) => a.localeCompare(b)).map(([column, runs]) => ({ column, runs, of: grouped.get(column).get(item.case).length }))
        }));
    const clean = cases.filter((id) => !disagreements.some((item) => item.case === id));
    return { columns, reference, live, shown, cases, parts, differences, components, disagreements, clean };
};

const MAX_ANSWERS = 4;

/**
 * The comparison as Markdown lines, for one split's records (report.js puts
 * it under each split's tables): the case score's parts, every component
 * against A+, and every disagreement with the planted truth. None without
 * records; a split with no live run or no A+ record says so.
 */
export const renderComparison = (records) => {
    const result = compareWithReference(records);
    if (!result.columns.length) return [];
    const { reference, live, shown } = result;
    const lines = ['### The case score\'s parts', '',
        'Beside the case score, unchanged: each column\'s score restricted to one group of its components (scoring.js `scoreComponents`). '
        + 'The feature-only score keeps the planted features alone, so cases with none are left out of it; the spurious-claim test and '
        + 'the symmetry fields, which every case has, are not in it. Each case weighs once; a live column minus A+ is taken on the cases both have.', '',
        `| Column | Cases | Case score | ${COMPONENT_GROUPS.map(([, label]) => label).join(' | ')} |`,
        `| --- | --- | --- | ${COMPONENT_GROUPS.map(() => '---').join(' | ')} |`];
    result.parts.forEach((part) => {
        const groups = COMPONENT_GROUPS.map(([group]) => {
            const { value, cases } = part.groups[group];
            return `${fmt(value)}${value != null && cases !== part.cases ? ` (${cases} case${cases === 1 ? '' : 's'})` : ''}`;
        });
        lines.push(`| ${part.column} | ${part.cases} | ${fmt(part.score)} | ${groups.join(' | ')} |`);
        const diff = result.differences.find((item) => item.column === part.column);
        if (diff) {
            lines.push(`| ${part.column} − ${REFERENCE_ARM} | ${diff.cases} paired | ${fmtDifference(diff.score)} | `
                + `${COMPONENT_GROUPS.map(([group]) => fmtDifference(diff.groups[group])).join(' | ')} |`);
        }
    });
    lines.push('');
    if (!reference) {
        lines.push(`No ${REFERENCE_ARM} record in this split: the comparison with the engines needs \`run.js arm-a\` records beside the runs.`, '');
    } else if (!live.length) {
        lines.push(`No live run (arm B or C) in this split yet: each live column's difference from ${REFERENCE_ARM} appears here and in the table below once there are.`, '');
    }

    if (reference) {
        lines.push(`### Component by component, against ${REFERENCE_ARM}`, '');
        lines.push(`Every component of the case score: each column's mean over its runs (1 right, 0.5 a quantity off, 0 wrong or missing), and each live column minus ${REFERENCE_ARM}.`, '');
        const header = shown.flatMap((column) => (live.includes(column) ? [column, `${column} − ${REFERENCE_ARM}`] : [column]));
        lines.push(`| Case | Component | ${header.join(' | ')} |`);
        lines.push(`| --- | --- | ${header.map(() => '---').join(' | ')} |`);
        result.components.forEach((row) => {
            const cells = shown.flatMap((column) => {
                const value = fmt(row.values[column]);
                if (!live.includes(column)) return [value];
                const diff = row.differences[column];
                return [value, diff != null && !isZero(diff) ? `**${fmtDifference(diff)}**` : fmtDifference(diff)];
            });
            lines.push(`| ${row.case} | ${row.id} | ${cells.join(' | ')} |`);
        });
        lines.push('');
        live.forEach((column) => {
            const diffs = result.components.map((row) => row.differences[column]).filter((value) => value != null);
            const higher = diffs.filter((value) => !isZero(value) && value > 0).length;
            const lower = diffs.filter((value) => !isZero(value) && value < 0).length;
            lines.push(`- ${column} against ${REFERENCE_ARM}: ${higher} component${higher === 1 ? '' : 's'} higher, ${lower} lower, ${diffs.length - higher - lower} the same.`);
        });
        if (live.length) lines.push('');
    }

    lines.push('### Every disagreement with the planted truth', '');
    lines.push('Each feature or fact a run\'s verdict got wrong, with the columns that made it (runs that did, of the column\'s runs on the case). '
        + 'Each costs a component of the case score; tolerated claims (true partial descriptions) are not disagreements.', '');
    if (result.disagreements.length) {
        lines.push('| Case | Disagreement | The answer | The truth | Made by (runs) |');
        lines.push('| --- | --- | --- | --- | --- |');
        result.disagreements.forEach((item) => {
            const answers = item.answers.length > MAX_ANSWERS
                ? `${item.answers.slice(0, MAX_ANSWERS).join('; ')}; … (${item.answers.length - MAX_ANSWERS} more)`
                : item.answers.join('; ');
            lines.push(`| ${item.case} | ${item.where} | ${answers} | ${item.truth} | ${item.made.map((m) => `${m.column} (${m.runs}/${m.of})`).join(', ')} |`);
        });
        lines.push('');
    } else {
        lines.push('No run disagrees with the planted truth.', '');
    }
    if (result.clean.length && result.disagreements.length) lines.push(`No disagreement in any run: ${result.clean.join(', ')}.`, '');
    return lines;
};
