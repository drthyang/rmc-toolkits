// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run module's browser adapter (copilotRun.js is its pure core): the run
// as App hands it to the AI Copilot. It takes the run the Dashboard publishes
// (`dashboard`: { files, structure, structurePath, hiddenPlots }), reads the
// run's <stem>.dat, builds the engine, and memoizes each field of the core on
// exactly its own inputs, so a field keeps its identity until they change:
//   symmetry: structure + toleranceA;  ladder: structure;
//   rValueFile: files + hiddenPlots + structurePath;  tools: engine + structure.
// A new PCA table (view.pcaSites) replaces `props` but not `tools` or `engine`.
//
// Ordering: until the Copilot is opened (`enabled`) and the Dashboard has
// published, the run has its name and no data (the Dashboard publishes only
// once the Copilot is open, so its parses and Live Data polls do not
// re-render App before that); then it fills in field by field. The view's
// own props (pcaSites, liveData) pass through from the start, as they always
// have. Its reads never throw: an unreadable .dat gives null, and an engine
// failure rejects that tool call only (the core's RangeError on a
// non-positive tolerance cannot arise from a ladder pick).

import { useEffect, useMemo, useRef, useState } from 'react';
import { isStaticMode, parseRunSettings } from './browserData';
import { requestAnalysis } from './useSiteCloud';
import {
    assembleCopilotRun, DEFAULT_VIEW, normalizeView, runElements, runEngineTools, runLadder, runPlotFiles,
    runRValueFile, runSymmetry
} from './copilotRun';

const NO_FILES = Object.freeze([]);

/**
 * The engine for the open run. A browser-loaded run (a picked folder or the
 * Demo) is analysed from its `.rmc6f` text, read lazily on the first request
 * and cached per file, in a worker of the Copilot's own: a stopped or
 * timed-out request terminates it (the next request starts a fresh one), and
 * the pages' shared worker never queues behind the Copilot. A typed server
 * folder goes to Flask, the request cancelled on abort. Null in the static
 * build with no run open. Rebuilt only when the run or the folder changes.
 */
const useRunEngine = ({ directory, localRun }) => {
    const textRef = useRef({ file: null, text: null });
    const workerRef = useRef(null);
    const localFile = localRun?.structureFile?.sourceFile || null;
    const canRequest = Boolean(localFile) || (!localRun && !isStaticMode());

    useEffect(() => () => {
        workerRef.current?.terminate();
        workerRef.current = null;
    }, []);

    return useMemo(() => {
        if (!canRequest) return null;
        const readText = async () => {
            if (textRef.current.file !== localFile) {
                textRef.current = { file: localFile, text: await localFile.text() };
            }
            return textRef.current.text;
        };
        const workerRequest = async (kind, params, signal) => {
            const text = await readText();
            if (!workerRef.current) {
                workerRef.current = new Worker(new URL('./workers/pcaKdeWorker.js', import.meta.url), { type: 'module' });
            }
            const worker = workerRef.current;
            const onAbort = () => {
                if (workerRef.current === worker) workerRef.current = null;
                worker.terminate();
            };
            signal?.addEventListener('abort', onAbort, { once: true });
            try {
                return await requestAnalysis(kind, params, { rmc6fText: text, worker, signal });
            } finally {
                signal?.removeEventListener('abort', onAbort);
            }
        };
        return (kind, params, { signal } = {}) => (localFile
            ? workerRequest(kind, params, signal)
            : requestAnalysis(kind, params, { directory, signal }));
    }, [canRequest, localFile, directory]);
};

/**
 * The parsed run-control settings (<stem>.dat) of a browser-loaded run, read
 * only once the Copilot is in use and again only when the file's path or time
 * changes (a Live Data poll that leaves it alone does not re-read it); null
 * otherwise (always in Flask mode), while it is read, or when unreadable.
 */
const useRunSettings = ({ localRun, enabled }) => {
    const settingsFile = enabled ? localRun?.settingsFile ?? null : null;
    const source = settingsFile?.sourceFile ?? null;
    const signature = source ? `${settingsFile.path}:${settingsFile.modified}` : '';
    // The last read, tagged with the file it came from: shown only while that
    // is still the file, so a read never outlives its file.
    const [read, setRead] = useState({ signature: '', settings: null });
    const sourceRef = useRef(null);
    useEffect(() => {
        sourceRef.current = source;
    }, [source]);
    useEffect(() => {
        const file = sourceRef.current;
        if (!signature || !file) return undefined;
        let cancelled = false;
        file.text()
            .then((text) => { if (!cancelled) setRead({ signature, settings: parseRunSettings(text) }); })
            .catch(() => { if (!cancelled) setRead({ signature, settings: null }); });
        return () => { cancelled = true; };
    }, [signature]);
    return read.signature === signature ? read.settings : null;
};

/**
 * The CopilotRun ({ props, tools, elements, engine, view }) for App: `props`
 * goes whole to CopilotChat. `directory` and `localRun` are App's open run,
 * `dashboard` what the Dashboard publishes (or null), `view` the state on
 * screen ({ toleranceA, pcaSites, liveData }; hiddenPlots come from the
 * Dashboard), `enabled` whether the Copilot has been opened.
 */
export default function useCopilotRun({ directory, localRun, dashboard = null, view = DEFAULT_VIEW, enabled = false }) {
    const engine = useRunEngine({ directory, localRun });
    const settings = useRunSettings({ localRun, enabled });

    const published = enabled ? dashboard : null;
    const files = published?.files ?? NO_FILES;
    const structure = published?.structure ?? null;
    const structurePath = published?.structurePath ?? null;
    const hiddenInput = published?.hiddenPlots ?? view?.hiddenPlots ?? DEFAULT_VIEW.hiddenPlots;
    const {
        toleranceA = DEFAULT_VIEW.toleranceA, pcaSites = null, liveData = false, ladder = true
    } = view || {};

    const hiddenPlots = useMemo(() => normalizeView({ hiddenPlots: hiddenInput }).hiddenPlots, [hiddenInput]);
    const normalized = useMemo(
        () => normalizeView({ toleranceA, hiddenPlots, pcaSites, liveData, ladder }),
        [toleranceA, hiddenPlots, pcaSites, liveData, ladder]
    );
    const plotFiles = useMemo(() => runPlotFiles(files), [files]);
    const rValueFile = useMemo(() => runRValueFile(plotFiles, hiddenPlots, structurePath), [plotFiles, hiddenPlots, structurePath]);
    const ladderList = useMemo(() => runLadder(structure, normalized.ladder), [structure, normalized.ladder]);
    const symmetry = useMemo(
        () => runSymmetry(structure, normalized.toleranceA, ladderList),
        [structure, normalized.toleranceA, ladderList]
    );
    const elements = useMemo(() => runElements(structure), [structure]);
    const tools = useMemo(() => runEngineTools(engine, structure, elements), [engine, structure, elements]);
    const name = localRun ? localRun.name : directory;
    const runSettings = published ? settings : null;

    return useMemo(() => assembleCopilotRun({
        name, plotFiles, rValueFile, structure, symmetry, settings: runSettings, tools, elements, engine, view: normalized
    }), [name, plotFiles, rValueFile, structure, symmetry, runSettings, tools, elements, engine, normalized]);
}
