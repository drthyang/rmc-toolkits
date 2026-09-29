// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useRef, useState } from 'react';
import axios from 'axios';
import Dashboard from './components/Dashboard';
import StructurePage from './components/StructurePage';
import PcaKdePage from './components/PcaKdePage';
import OrientationPage from './components/OrientationPage';
import BondGeometryPage from './components/BondGeometryPage';
import AutoStogPage from './components/AutoStogPage';
import { AssistantPage } from './llm';
import API_BASE_URL from './api';
import {
  buildLocalRun,
  buildLocalRunFromHandle,
  fileSignature,
  isStaticMode,
  loadDemoRun,
  supportsFileSystemAccess,
  WATCH_INTERVAL_MS,
} from './browserData';
import { SymTolContext } from './symTolContext';
import { IconButton, Segmented, SegmentedButton } from './ui';
import InfoBadge from './ui/InfoBadge';
import './App.css';

const REPO_URL = 'https://github.com/drthyang/rmc-toolkits';
const STATUS_TIMEOUT_MS = 15000;
// Auto StoG is still under development; flip to true to expose the tab again.
const SHOW_AUTO_STOG = false;

function App() {
  const [activePage, setActivePage] = useState('dashboard');
  const [visitedPages, setVisitedPages] = useState({ autostog: false, dashboard: true, structure: false, ellipsoids: false, orientation: false, geometry: false, assistant: false });
  // Parsed run context published by the Dashboard, consumed by the AI Assistant page.
  const [assistantRun, setAssistantRun] = useState(null);
  // Per-site PCA ellipsoid table published by the PCA Ellipsoid page (once opened),
  // also fed to the AI Assistant so it can reason about the thermal displacements.
  const [pcaSites, setPcaSites] = useState(null);
  const [currentDirectory, setCurrentDirectory] = useState('data');
  const [draftDirectory, setDraftDirectory] = useState('data');
  const [browseStatus, setBrowseStatus] = useState(null);
  const [localRun, setLocalRun] = useState(null);
  const [localLoading, setLocalLoading] = useState(false);
  // Whether the bundled demo run is the one currently loaded (drives the header toggle).
  const [demoActive, setDemoActive] = useState(false);
  const [watchFiles, setWatchFiles] = useState(false);
  // Flask-mode configuration epoch: bumped when the run's .rmc6f changes on disk
  // (see the effect below); the analysis pages take it as `dataEpoch` and
  // re-read their data in place when it changes.
  const [configEpoch, setConfigEpoch] = useState(0);
  // Bumped by every Load / Select Folder, even of the folder already shown:
  // React skips an unchanged currentDirectory, so this is what makes loading
  // the same folder again re-check its .rmc6f (with Live Data off, the only
  // way short of a browser reload).
  const [loadRequest, setLoadRequest] = useState(0);
  // Shared "Detected SG" tolerance, so the ladder selection persists across pages.
  const symTolState = useState(0.2);
  const directoryInputRef = useRef(null);
  const dirHandleRef = useRef(null);
  const lastSignatureRef = useRef('');
  const configSignatureRef = useRef({ directory: null, signature: null });
  const runIdRef = useRef(0);
  const shellRef = useRef(null);
  const headerRef = useRef(null);
  const staticMode = isStaticMode();
  const fsAccess = staticMode && supportsFileSystemAccess();

  useEffect(() => {
    document.documentElement.dataset.theme = 'light';
    localStorage.setItem('rmc-theme', 'light');
  }, []);

  useEffect(() => {
    setVisitedPages((current) => (
      current[activePage] ? current : { ...current, [activePage]: true }
    ));
  }, [activePage]);

  // Phones and short windows scroll the whole shell and pin the header by its
  // tab row (App.css). Measure how far the header may slide up — to just above
  // the tabs, keeping its bottom padding above them too — and publish it on the
  // shell as --header-pin (the sticky offset) and --header-pinned-h (the strip
  // left on screen, which a fill-height page such as the AI Assistant
  // subtracts). Elsewhere the header is static and the variables are unused.
  useEffect(() => {
    const shell = shellRef.current;
    const header = headerRef.current;
    if (!shell || !header || typeof ResizeObserver === 'undefined') return undefined;
    const update = () => {
      const tabs = header.querySelector('.page-tabs');
      if (!tabs) return;
      const headerBox = header.getBoundingClientRect();
      const tabsTop = tabs.getBoundingClientRect().top - headerBox.top;
      const padding = parseFloat(window.getComputedStyle(header).paddingBottom) || 0;
      const pin = Math.max(0, Math.round(tabsTop - padding));
      shell.style.setProperty('--header-pin', `${-pin}px`);
      shell.style.setProperty('--header-pinned-h', `${Math.round(headerBox.height) - pin}px`);
    };
    update();
    const observer = new ResizeObserver(update);
    observer.observe(header);
    return () => observer.disconnect();
  }, []);

  // Keep the active tab in view when the tab row scrolls sideways (tablets
  // and phones).
  useEffect(() => {
    const tabs = headerRef.current?.querySelector('.page-tabs');
    const active = tabs?.querySelector('button.is-active');
    if (!tabs || !active || tabs.scrollWidth <= tabs.clientWidth) return;
    const tabsBox = tabs.getBoundingClientRect();
    const activeBox = active.getBoundingClientRect();
    const offset = activeBox.left - tabsBox.left - (tabsBox.width - activeBox.width) / 2;
    tabs.scrollBy?.({ left: offset, behavior: 'smooth' });
  }, [activePage]);

  const handlePageChange = (page) => {
    setVisitedPages((current) => (
      current[page] ? current : { ...current, [page]: true }
    ));
    setActivePage(page);
    // Scrolling shell (phones, short windows): open the new page at its top,
    // right under the pinned tabs, not at the previous page's scroll depth.
    const shell = shellRef.current;
    if (shell && shell.scrollTop > 0) {
      const pin = -parseFloat(shell.style.getPropertyValue('--header-pin')) || 0;
      if (shell.scrollTop > pin) shell.scrollTop = pin;
    }
  };

  useEffect(() => {
    if (!browseStatus || browseStatus.kind === 'loading') return undefined;
    const timer = window.setTimeout(() => setBrowseStatus(null), STATUS_TIMEOUT_MS);
    return () => window.clearTimeout(timer);
  }, [browseStatus]);

  // Static-mode Live Data: re-read the picked folder handle on an interval and rebuild
  // localRun only when files change, mirroring the Flask /api/files poll. Chromium only.
  useEffect(() => {
    if (!fsAccess || !watchFiles || !dirHandleRef.current) return undefined;
    let cancelled = false;
    let inFlight = false;
    const poll = async () => {
      if (inFlight || cancelled) return;
      inFlight = true;
      try {
        const nextRun = await buildLocalRunFromHandle(dirHandleRef.current);
        const signature = fileSignature(nextRun.files);
        if (!cancelled && signature !== lastSignatureRef.current) {
          lastSignatureRef.current = signature;
          // Same folder, updated files: keep runId so the dashboard refreshes in place.
          setLocalRun({ ...nextRun, runId: runIdRef.current });
        }
      } catch (error) {
        if (!cancelled) {
          setWatchFiles(false);
          setBrowseStatus({
            kind: 'error',
            text: error.message || 'Lost access to the run folder — re-select it to resume Live Data.',
          });
        }
      } finally {
        inFlight = false;
      }
    };
    const interval = window.setInterval(poll, WATCH_INTERVAL_MS);
    return () => {
      cancelled = true;
      window.clearInterval(interval);
    };
  }, [fsAccess, watchFiles]);

  // Flask-mode Live Data for the analysis pages. The Dashboard polls /api/files
  // itself, but the Atomic Density, Bond Geometry, PCA Ellipsoid and Displacement
  // Directions pages fetch from the backend on demand, and the backend always
  // reads the file currently on disk: after RMCProfile saves a new .rmc6f, a page
  // left alone would keep its old site table / slab points while its next request
  // (a slider move, a site click) came from the new configuration — two
  // configurations mixed in one view. So watch the .rmc6f signature in the same
  // listing (checked on every Load of a folder -- the same one included -- then
  // every poll while Live Data is on) and bump configEpoch when it changes: the
  // pages receive it as `dataEpoch`, which is in the dependencies of their
  // backend fetches, so they re-read everything from the one new file in
  // place — keeping the picks and view settings that still apply, as static
  // mode does when a picked folder's files change (docs/algorithms/notation.md
  // §3c). A browser-loaded run (Demo, picked folder) is a snapshot and is not
  // watched here.
  useEffect(() => {
    if (staticMode || localRun) return undefined;
    let cancelled = false;
    let inFlight = false;
    const check = async () => {
      if (inFlight || cancelled) return;
      inFlight = true;
      try {
        const response = await axios.get(`${API_BASE_URL}/api/files`, {
          params: { dir: currentDirectory || '.' }
        });
        if (cancelled) return;
        const structureFiles = (response.data.files || []).filter(
          (file) => file.type === 'file' && /\.rmc6f$/i.test(file.name)
        );
        const signature = fileSignature(structureFiles);
        const known = configSignatureRef.current;
        if (known.directory === currentDirectory && known.signature !== null && known.signature !== signature) {
          setConfigEpoch((epoch) => epoch + 1);
        }
        configSignatureRef.current = { directory: currentDirectory, signature };
      } catch {
        // Listing errors surface through the Dashboard's own poll.
      } finally {
        inFlight = false;
      }
    };
    check();
    if (!watchFiles) {
      return () => { cancelled = true; };
    }
    const interval = window.setInterval(check, WATCH_INTERVAL_MS);
    return () => {
      cancelled = true;
      window.clearInterval(interval);
    };
  }, [staticMode, localRun, watchFiles, currentDirectory, loadRequest]);

  const handleDirectorySubmit = (event) => {
    event.preventDefault();
    const nextDirectory = draftDirectory.trim() || '.';
    setCurrentDirectory(nextDirectory);
    setLoadRequest((count) => count + 1);
  };

  const handleNativeBrowse = async () => {
    setBrowseStatus({ kind: 'loading', text: 'Opening folder picker…' });
    try {
      const response = await axios.post(`${API_BASE_URL}/api/dialog/folder`, {
        dir: draftDirectory || currentDirectory || '.'
      });
      const nextPath = response.data.path;
      setDraftDirectory(nextPath);
      setCurrentDirectory(nextPath);
      setLoadRequest((count) => count + 1);
      setBrowseStatus(null);
    } catch (err) {
      const message = err.response?.data?.error || 'Could not open the folder picker';
      if (message !== 'Folder selection cancelled') {
        setBrowseStatus({ kind: 'error', text: message });
      } else {
        setBrowseStatus(null);
      }
    }
  };

  const handleLocalFiles = async (event) => {
    const selectedFiles = event.target.files;
    if (!selectedFiles?.length) return;
    setLocalLoading(true);
    setBrowseStatus({ kind: 'loading', text: 'Indexing selected folder…' });
    try {
      const nextRun = await buildLocalRun(selectedFiles);
      runIdRef.current += 1;
      setLocalRun({ ...nextRun, runId: runIdRef.current });
      setDemoActive(false);
      setCurrentDirectory(nextRun.name);
      setDraftDirectory(nextRun.name);
      setBrowseStatus(null);
      setActivePage('dashboard');
    } catch (error) {
      setBrowseStatus({ kind: 'error', text: error.message || 'Could not read the selected files' });
    } finally {
      setLocalLoading(false);
      event.target.value = '';
    }
  };

  const handleSelectFolderFsAccess = async () => {
    setLocalLoading(true);
    setBrowseStatus({ kind: 'loading', text: 'Opening folder picker…' });
    try {
      const handle = await window.showDirectoryPicker();
      dirHandleRef.current = handle;
      const nextRun = await buildLocalRunFromHandle(handle);
      lastSignatureRef.current = fileSignature(nextRun.files);
      runIdRef.current += 1;
      setLocalRun({ ...nextRun, runId: runIdRef.current });
      setDemoActive(false);
      setCurrentDirectory(nextRun.name);
      setDraftDirectory(nextRun.name);
      setBrowseStatus(null);
      setActivePage('dashboard');
    } catch (error) {
      if (error?.name === 'AbortError') {
        setBrowseStatus(null);
      } else {
        setBrowseStatus({ kind: 'error', text: error.message || 'Could not read the selected folder' });
      }
    } finally {
      setLocalLoading(false);
    }
  };

  const handleToggleDemo = async () => {
    // Second click: tear the demo run back down and return to the empty state.
    if (demoActive) {
      dirHandleRef.current = null;
      lastSignatureRef.current = '';
      setLocalRun(null);
      setDemoActive(false);
      setCurrentDirectory('data');
      setDraftDirectory('data');
      setBrowseStatus(null);
      return;
    }
    setLocalLoading(true);
    setBrowseStatus({ kind: 'loading', text: 'Loading demo dataset…' });
    try {
      const nextRun = await loadDemoRun();
      runIdRef.current += 1;
      setLocalRun({ ...nextRun, runId: runIdRef.current });
      setDemoActive(true);
      setCurrentDirectory(nextRun.name);
      setDraftDirectory(nextRun.name);
      setBrowseStatus(null);
      setActivePage('dashboard');
    } catch (error) {
      setBrowseStatus({ kind: 'error', text: error.message || 'Could not load the demo dataset' });
    } finally {
      setLocalLoading(false);
    }
  };

  const handleStaticLiveDataNotice = () => {
    setBrowseStatus({
      kind: 'info',
      text: 'Live Data needs Chrome, Edge, Arc or Opera — or',
      link: {
        href: REPO_URL,
        label: 'the local app'
      }
    });
  };

  const renderBrowseStatus = () => {
    if (!browseStatus) return null;
    return (
      <div className={`ui-status is-${browseStatus.kind}`} role="status">
        <span>
          {browseStatus.text}
          {browseStatus.link && (
            <>
              {' '}
              <a href={browseStatus.link.href} target="_blank" rel="noreferrer">
                {browseStatus.link.label}
              </a>
            </>
          )}
        </span>
        <IconButton
          variant="close"
          onClick={() => setBrowseStatus(null)}
          aria-label="Close notification"
          title="Close"
        >
          &times;
        </IconButton>
      </div>
    );
  };

  return (
    <div className="app-container">
      <main className="main-content" ref={shellRef}>
        <header className="app-header" ref={headerRef}>
          <div className="header-primary">
            <div className="brand-row">
              <div className="brand-mark" aria-hidden="true">
                <svg className="brand-mark-icon" viewBox="0 0 100 100">
                  <path
                    d="M14 62 C22 62 27 24 37 24 C47 24 45 66 54 66 C62 66 61 44 69 44 C76 44 78 57 86 57"
                    transform="translate(0,5)"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth="9"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
              </div>
              <div className="brand-copy">
                <h1>
                  RMCProfile
                  <span>Workbench</span>
                </h1>
              </div>
            </div>
            <Segmented as="nav" variant="nav" className="page-tabs" aria-label="Workspace pages">
              {SHOW_AUTO_STOG && (
                <SegmentedButton
                  active={activePage === 'autostog'}
                  onClick={() => handlePageChange('autostog')}
                >
                  Auto StoG
                </SegmentedButton>
              )}
              <SegmentedButton
                active={activePage === 'dashboard'}
                onClick={() => handlePageChange('dashboard')}
              >
                Dashboard
              </SegmentedButton>
              <SegmentedButton
                active={activePage === 'structure'}
                onClick={() => handlePageChange('structure')}
              >
                Atomic Density
              </SegmentedButton>
              <SegmentedButton
                active={activePage === 'geometry'}
                onClick={() => handlePageChange('geometry')}
              >
                Bond Geometry
              </SegmentedButton>
              <SegmentedButton
                active={activePage === 'ellipsoids'}
                onClick={() => handlePageChange('ellipsoids')}
              >
                PCA Ellipsoid
              </SegmentedButton>
              <SegmentedButton
                active={activePage === 'orientation'}
                onClick={() => handlePageChange('orientation')}
              >
                Displacement Directions
              </SegmentedButton>
              <SegmentedButton
                active={activePage === 'assistant'}
                onClick={() => handlePageChange('assistant')}
              >
                AI Assistant
              </SegmentedButton>
            </Segmented>
          </div>
          <div className="header-actions">
          {fsAccess ? (
            <div className="path-controls">
              <label className="ui-switch-outline">
                <input
                  type="checkbox"
                  checked={watchFiles}
                  onChange={(event) => setWatchFiles(event.target.checked)}
                />
                <span aria-hidden="true" />
                <b>Live Data</b>
              </label>
              <div className="ui-fieldbar ui-fieldbar--readonly path-bar local-file-bar">
                <label>Local run</label>
                <div className="ui-fieldbar__value">{localRun?.name || 'No folder selected'}</div>
                <button
                  type="button"
                  onClick={handleSelectFolderFsAccess}
                  disabled={localLoading}
                  title="Read locally — nothing is uploaded"
                >
                  {localLoading ? 'Reading' : 'Select Folder'}
                </button>
              </div>
              <InfoBadge label="Where your files go" align="end">
                <p>Files are read locally in your browser and never uploaded. The browser's own folder picker may still label its button “Upload”.</p>
              </InfoBadge>
            </div>
          ) : staticMode ? (
            <div className="path-controls">
              <button
                type="button"
                className="ui-switch-outline ui-switch-outline--button"
                onClick={handleStaticLiveDataNotice}
                aria-pressed="false"
                title="Live Data needs a Chromium browser or the local app"
              >
                <span aria-hidden="true" />
                <b>Live Data</b>
              </button>
              <div className="ui-fieldbar ui-fieldbar--readonly path-bar local-file-bar">
                <label htmlFor="local-run-files">Local run</label>
                <input
                  ref={directoryInputRef}
                  id="local-run-files"
                  className="ui-visually-hidden"
                  type="file"
                  multiple
                  webkitdirectory=""
                  onChange={handleLocalFiles}
                />
                <div className="ui-fieldbar__value">{localRun?.name || 'No folder selected'}</div>
                <button
                  type="button"
                  onClick={() => directoryInputRef.current?.click()}
                  disabled={localLoading}
                  title="Read locally — nothing is uploaded"
                >
                  {localLoading ? 'Reading' : 'Select Folder'}
                </button>
              </div>
              <InfoBadge label="Where your files go" align="end">
                <p>Files are read locally in your browser and never uploaded. The browser's own folder picker may still label its button “Upload”.</p>
              </InfoBadge>
            </div>
          ) : (
            <div className="path-controls">
              <label className="ui-switch-outline">
                <input
                  type="checkbox"
                  checked={watchFiles}
                  onChange={(event) => setWatchFiles(event.target.checked)}
                />
                <span aria-hidden="true" />
                <b>Live Data</b>
              </label>
              <form className="ui-fieldbar path-bar" onSubmit={handleDirectorySubmit}>
                <label htmlFor="data-path">Run folder</label>
                <input
                  id="data-path"
                  type="text"
                  value={draftDirectory}
                  onChange={(event) => setDraftDirectory(event.target.value)}
                  spellCheck="false"
                />
                <button
                  type="button"
                  className="ui-fieldbar__ghost"
                  onClick={handleNativeBrowse}
                >
                  Select Folder
                </button>
                <button type="submit">
                  Load
                </button>
              </form>
            </div>
          )}
          <button
            type="button"
            className={`ui-btn-brand demo-button${demoActive ? ' is-active' : ''}`}
            onClick={handleToggleDemo}
            disabled={localLoading}
            aria-pressed={demoActive}
            title={demoActive ? 'Clear the demo dataset' : 'Load a bundled GaTa4Se8 250 K example run'}
          >
            {localLoading ? 'Loading…' : 'Demo'}
          </button>
          </div>
          {renderBrowseStatus()}
        </header>
        <SymTolContext.Provider value={symTolState}>
        <div className="workspace-pages">
          {SHOW_AUTO_STOG && visitedPages.autostog && (
            <div
              className={`workspace-page${activePage === 'autostog' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'autostog'}
            >
              {/* Auto StoG is pre-processing: page-local uploads, independent
                  of the run folder the post-processing pages share. */}
              <AutoStogPage />
            </div>
          )}
          {visitedPages.dashboard && (
            <div
              className={`workspace-page${activePage === 'dashboard' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'dashboard'}
            >
              <Dashboard
                directory={currentDirectory}
                localRun={localRun}
                watchFiles={watchFiles}
                wantAssistantData={visitedPages.assistant}
                onRunContextChange={setAssistantRun}
              />
            </div>
          )}
          {visitedPages.structure && (
            <div
              className={`workspace-page${activePage === 'structure' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'structure'}
            >
              <StructurePage dataEpoch={configEpoch} directory={currentDirectory} localRun={localRun} theme="light" />
            </div>
          )}
          {visitedPages.ellipsoids && (
            <div
              className={`workspace-page${activePage === 'ellipsoids' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'ellipsoids'}
            >
              <PcaKdePage dataEpoch={configEpoch} directory={currentDirectory} localRun={localRun} theme="light" onSitesChange={setPcaSites} />
            </div>
          )}
          {visitedPages.orientation && (
            <div
              className={`workspace-page${activePage === 'orientation' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'orientation'}
            >
              {/* Displacement-direction histogram — independent of the PCA page
                  (shares only the site picker via useSiteCloud). */}
              <OrientationPage dataEpoch={configEpoch} directory={currentDirectory} localRun={localRun} />
            </div>
          )}
          {visitedPages.geometry && (
            <div
              className={`workspace-page${activePage === 'geometry' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'geometry'}
            >
              {/* Bond-angle (triplet) distribution + bond-length/coordination
                  statistics — the RMCProfile `triplets` workflow. */}
              <BondGeometryPage dataEpoch={configEpoch} directory={currentDirectory} localRun={localRun} />
            </div>
          )}
          {visitedPages.assistant && (
            <div
              className={`workspace-page${activePage === 'assistant' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'assistant'}
            >
              <AssistantPage
                runName={assistantRun?.runName ?? (localRun ? localRun.name : currentDirectory)}
                plotFiles={assistantRun?.plotFiles ?? []}
                rValueFile={assistantRun?.rValueFile ?? null}
                structure={assistantRun?.structure ?? null}
                symmetry={assistantRun?.symmetry ?? null}
                runSettings={assistantRun?.runSettings ?? null}
                pcaSites={pcaSites}
                liveData={watchFiles}
              />
            </div>
          )}
        </div>
        </SymTolContext.Provider>
      </main>
    </div>
  );
}

export default App;
