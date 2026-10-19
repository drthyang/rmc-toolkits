// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Shared by the model-visible snapshot tests (this folder): what the
// Copilot's model is sent, rendered as small, stable, readable files. Every
// change to these files is a change to what the model sees.
//
//   - the system-prompt variants (copilotSystemPrompt);
//   - the Demo run as the Dashboard holds it and the run module
//     (copilotRun.js) hands it to the Copilot (the dock's run);
//   - a fake model port that records each request and answers at once;
//   - the renderers: a prompt with a line break after each sentence, a
//     request as text, a tool message as JSON.

import { createHash } from 'node:crypto';
import { readdirSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
    DASHBOARD_MAX_POINTS, isDashboardPlotFile, loadDemoRun, parseDashboardPlotFile, parseRunSettings
} from '../../browserData';
import { copilotRun, DEFAULT_VIEW } from '../../copilotRun';
import { copilotSkills } from '../../copilotSkills';
import { parseSkillFiles } from '../../llm/copilot/skills';
import { copilotSystemPrompt } from '../../llm/prompts/copilot';
import { handlePcaMessage } from '../../workers/pcaKdeWorker';
import { handleStructureMessage } from '../../workers/localStructureWorker';

// --- the system prompt -----------------------------------------------------------

/** The skills the app bundles (copilotSkills.js), parsed as CopilotChat parses them. */
export const bundledSkills = () => parseSkillFiles(copilotSkills.files, { available: copilotSkills.references });

// Every system prompt the code can build: tools off (one prompt, whatever
// else is asked), or tools on with each combination of the view tools,
// ask_user and the skills' index.
export const PROMPT_VARIANTS = [
    'no-tools',
    'tools', 'tools-ask', 'tools-skills', 'tools-ask-skills',
    'tools-view', 'tools-view-ask', 'tools-view-skills', 'tools-view-ask-skills'
];

export const variantFlags = (name) => ({
    toolsEnabled: name !== 'no-tools',
    viewTools: name.includes('-view'),
    askUser: name.includes('-ask'),
    skills: name.includes('-skills')
});

/** name → the exact system prompt, with `index` as the skills' index. */
export const promptVariants = (index) => new Map(PROMPT_VARIANTS.map((name) => {
    const flags = variantFlags(name);
    return [name, copilotSystemPrompt({ ...flags, skills: flags.skills ? index : '' })];
}));

const sha256 = (text) => createHash('sha256').update(text, 'utf8').digest('hex').slice(0, 16);

// A line break after each sentence's full stop ("e.g." and "i.e." are not
// sentences): the text is otherwise exact, and its length and hash say what
// the model receives.
export const readable = (text) => text.replace(/(?<!\b(?:e\.g|i\.e|etc|vs))\. (?=\S)/g, '.\n');

export const exactLine = (text) => `${text.length} characters, sha256 ${sha256(text)}`;

export const renderPrompt = (name, text) => {
    const flags = variantFlags(name);
    const call = flags.toolsEnabled
        ? `copilotSystemPrompt({ toolsEnabled: true, viewTools: ${flags.viewTools}, askUser: ${flags.askUser}, `
            + `skills: ${flags.skills ? 'skillIndex(<the bundled .claude/skills>)' : '\'\''} })`
        : 'copilotSystemPrompt({ toolsEnabled: false }) — COPILOT_PROMPT_NO_TOOLS';
    return [
        `# The Copilot's system prompt, variant "${name}"`,
        `# ${call}`,
        `# The model receives one string: ${exactLine(text)}.`,
        '# Here a line break follows each sentence\'s full stop, for reading.',
        '',
        readable(text),
        ''
    ].join('\n');
};

// --- the Demo run, as the Dashboard hands it to the Copilot --------------------

// A path built with node:path: Vite rewrites `new URL(…, import.meta.url)`
// into a served asset URL under jsdom, which node:fs cannot read.
const HERE = dirname(fileURLToPath(import.meta.url));
const DEMO = join(HERE, '../../../public/demo');
/** An app source file's text (src/<path>), to check what this file mirrors. */
export const appSource = (path) => readFileSync(join(HERE, '../..', path), 'utf8');
// The Demo button stamps each file with Date.now(); pinned, so nothing about
// the moment of the test can reach a snapshot.
export const DEMO_TIME = Date.UTC(2026, 0, 1);

/** The Demo folder as a picked folder's entries ({ path: 'Demo/<name>', file }), from public/demo. */
export const demoEntries = () => readdirSync(DEMO).sort().map((name) => ({
    path: `Demo/${name}`,
    file: new File([readFileSync(join(DEMO, name))], name, { lastModified: DEMO_TIME })
}));

/** An engine over .rmc6f text: the Copilot worker's handler in-thread, with the worker's message. */
export const inThreadEngine = (rmc6fText) => (kind, params = {}) => handlePcaMessage({ ...params, kind }, async () => rmc6fText);

/**
 * A browser-loaded run as the Dashboard holds it once everything is parsed,
 * as a LoadedRun (copilotRun.js): its files with each plot file parsed
 * (parseDashboardPlotFile), the structure worker's parse (maxPoints 100), the
 * structure file's path, the parsed .dat as useCopilotRun reads it, and an
 * in-thread engine.
 */
export const dashboardLoadedRun = async (localRun) => {
    const files = await Promise.all(localRun.files.map((file) => (isDashboardPlotFile(file) ? parseDashboardPlotFile(file) : file)));
    const structure = await handleStructureMessage({ file: localRun.structureFile, maxPoints: DASHBOARD_MAX_POINTS });
    const settings = localRun.settingsFile ? parseRunSettings(await localRun.settingsFile.sourceFile.text()) : null;
    const rmc6fText = await localRun.structureFile.sourceFile.text();
    return {
        loaded: {
            name: localRun.name, files, structure, structurePath: localRun.structureFile.path, settings, engine: inThreadEngine(rmc6fText)
        },
        rmc6fText
    };
};

/**
 * The Demo as App.jsx passes it to CopilotChat once the Copilot is open: the
 * Demo button's files (browserData.loadDemoRun, served from public/demo) as
 * the Dashboard holds them, through the run module's core (copilotRun.js) at
 * the default view — the Dashboard's plot files in its order, its χ² history,
 * the symmetry at the card's tolerance with its ladder, the parsed .dat.
 * `run` is the eight CopilotChat props, `copilot` the whole CopilotRun,
 * `loaded` its LoadedRun; `localRun` is returned too.
 */
export const loadDemoAssistantRun = async () => {
    const fetchBefore = globalThis.fetch;
    const nowBefore = Date.now;
    globalThis.fetch = async (url) => {
        const name = decodeURIComponent(String(url).split('/').pop());
        const bytes = readFileSync(join(DEMO, name));
        return { ok: true, status: 200, blob: async () => new Blob([bytes]) };
    };
    Date.now = () => DEMO_TIME;
    let localRun;
    try {
        localRun = await loadDemoRun();
    } finally {
        globalThis.fetch = fetchBefore;
        Date.now = nowBefore;
    }
    const { loaded, rmc6fText } = await dashboardLoadedRun(localRun);
    const copilot = copilotRun(loaded, DEFAULT_VIEW);
    return { localRun, rmc6fText, loaded, copilot, run: copilot.props };
};

// --- a fake model port ------------------------------------------------------------

const clone = (value) => (value === undefined ? undefined : JSON.parse(JSON.stringify(value)));

/**
 * A request as the model port received it, frozen at the call (the loop
 * keeps appending to its `messages`): every field, `signal` and `session`
 * by kind only.
 */
export const freezeRequest = (request) => Object.fromEntries(Object.entries(request).map(([key, value]) => {
    if (key === 'signal') return [key, value ? 'AbortSignal' : value];
    if (key === 'session') return [key, value && typeof value === 'object' ? 'session object (contextGuard.js context sizes)' : value];
    return [key, clone(value)];
}));

/**
 * A model that records every request and replies from `replies` in order
 * (the last one repeats): `{ calls: [[name, args], …] }` or `{ text }`.
 * `callModel(request)` records a provider call as it is made (the dock's
 * streamCopilotChat); `port` is a model port (provider/model.js) that records
 * each intent it is sent with its signal, leaving out a null `retry` (no
 * retry: the port's instruction, never sent to the model).
 */
export const recordingModel = (replies = [{ text: 'Recorded.\n\nOutcome: a recorded request.' }]) => {
    const requests = [];
    const callModel = (request) => {
        requests.push(freezeRequest(request));
        const reply = replies[Math.min(requests.length - 1, replies.length - 1)];
        return (async function* stream() {
            if (reply.text) yield { content: reply.text };
            for (const [index, [name, args]] of (reply.calls || []).entries()) {
                yield { toolCall: { index, id: `call_${requests.length}_${index}`, name, arguments: JSON.stringify(args) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
    const port = {
        stream: ({ retry = null, ...intent }, { signal } = {}) => callModel({ ...intent, ...(retry ? { retry } : {}), signal })
    };
    return { callModel, port, requests };
};

// --- a request as text --------------------------------------------------------------

const sameJson = (a, b) => JSON.stringify(a) === JSON.stringify(b);

/**
 * One request as readable text: its fields, the tools by name (their full
 * definitions are in `toolsFile`), and every message verbatim. A system
 * message equal to a prompt variant names its file; with `base`, a message
 * equal to the base request's message at the same index says so.
 */
export const renderRequest = (request, { title, port, variants, toolsFile, base = null, baseName = '' }) => {
    const lines = [`# ${title}`, `# Captured at ${port}.`, '', '## Request fields'];
    for (const key of Object.keys(request).sort()) {
        if (key === 'messages' || key === 'tools') continue;
        lines.push(`${key}: ${JSON.stringify(request[key])}`);
    }
    if (!Object.hasOwn(request, 'temperature')) lines.push('temperature: (not in the request)');
    const tools = request.tools || [];
    if (tools.length) {
        lines.push(`tools: ${tools.length}, in the order sent (each definition in full: ${toolsFile})`);
        tools.forEach((tool, index) => lines.push(`  ${index + 1}. ${tool.function?.name ?? JSON.stringify(tool)}`));
    } else {
        lines.push(`tools: ${JSON.stringify(request.tools)}`);
    }
    (request.messages || []).forEach((message, index) => {
        const extra = Object.keys(message).filter((key) => key !== 'role' && key !== 'content');
        lines.push('', `## messages[${index}] · ${message.role}${extra.length ? ` · ${extra.join(', ')}` : ''}`);
        extra.forEach((key) => lines.push(`${key}: ${JSON.stringify(message[key])}`));
        const content = message.content;
        if (base && sameJson(base.messages?.[index], message)) {
            lines.push(`(identical to ${baseName} messages[${index}])`);
            return;
        }
        const variant = message.role === 'system' && typeof content === 'string'
            ? [...variants.entries()].find(([, text]) => text === content)?.[0]
            : null;
        if (variant) {
            lines.push(`(the system prompt, variant "${variant}": system-prompt/${variant}.txt, ${exactLine(content)})`);
        } else if (message.role === 'system' && typeof content === 'string') {
            lines.push(`(matches no system-prompt variant; ${exactLine(content)}; a line break after each sentence)`, readable(content));
        } else {
            lines.push(typeof content === 'string' ? content : JSON.stringify(content));
        }
    });
    lines.push('');
    return lines.join('\n');
};

/**
 * JSON as the snapshot files hold it: two-space indent, each array of plain
 * values on one line, a final newline. (A bracket followed by a line break
 * never sits inside a JSON string, which escapes its line breaks, so only
 * real arrays are joined; runContext.js compactJson relies on the same.)
 */
export const jsonFile = (value) => `${JSON.stringify(value, null, 2)
    .replace(/\[\n\s*([^[\]{}]*?)\n\s*\]/g, (_, items) => `[${items.split(/,\n\s*/).join(', ')}]`)}\n`;
