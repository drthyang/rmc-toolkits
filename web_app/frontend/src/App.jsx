import React, { useEffect, useRef, useState } from 'react';
import axios from 'axios';
import Dashboard from './components/Dashboard';
import StructurePage from './components/StructurePage';
import API_BASE_URL from './api';
import {
  buildLocalRun,
  buildLocalRunFromHandle,
  fileSignature,
  isStaticMode,
  supportsFileSystemAccess,
  WATCH_INTERVAL_MS,
} from './browserData';
import './App.css';

const REPO_URL = 'https://github.com/drthyang/rmc-toolkits';
const STATUS_TIMEOUT_MS = 15000;

function App() {
  const [activePage, setActivePage] = useState('dashboard');
  const [visitedPages, setVisitedPages] = useState({ dashboard: true, structure: false });
  const [currentDirectory, setCurrentDirectory] = useState('data');
  const [draftDirectory, setDraftDirectory] = useState('data');
  const [browseStatus, setBrowseStatus] = useState(null);
  const [localRun, setLocalRun] = useState(null);
  const [localLoading, setLocalLoading] = useState(false);
  const [watchFiles, setWatchFiles] = useState(false);
  const directoryInputRef = useRef(null);
  const dirHandleRef = useRef(null);
  const lastSignatureRef = useRef('');
  const runIdRef = useRef(0);
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

  const handlePageChange = (page) => {
    setVisitedPages((current) => (
      current[page] ? current : { ...current, [page]: true }
    ));
    setActivePage(page);
  };

  // Startup reminder (static mode): reassure users that picking a folder keeps
  // files on their device, since the browser's native picker may say "Upload".
  useEffect(() => {
    if (!staticMode) return;
    setBrowseStatus({
      kind: 'info',
      text: 'Files stay on your device — they are read locally and never uploaded. Your browser may still label the picker button “Upload”.',
    });
  }, [staticMode]);

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
            text: error.message || 'Lost access to the run folder. Re-select it to resume Live Data.',
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

  const handleDirectorySubmit = (event) => {
    event.preventDefault();
    const nextDirectory = draftDirectory.trim() || '.';
    setCurrentDirectory(nextDirectory);
  };

  const handleNativeBrowse = async () => {
    setBrowseStatus({ kind: 'loading', text: 'Opening folder picker...' });
    try {
      const response = await axios.post(`${API_BASE_URL}/api/dialog/folder`, {
        dir: draftDirectory || currentDirectory || '.'
      });
      const nextPath = response.data.path;
      setDraftDirectory(nextPath);
      setCurrentDirectory(nextPath);
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
    setBrowseStatus({ kind: 'loading', text: 'Indexing selected folder...' });
    try {
      const nextRun = await buildLocalRun(selectedFiles);
      runIdRef.current += 1;
      setLocalRun({ ...nextRun, runId: runIdRef.current });
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
    setBrowseStatus({ kind: 'loading', text: 'Opening folder picker...' });
    try {
      const handle = await window.showDirectoryPicker();
      dirHandleRef.current = handle;
      const nextRun = await buildLocalRunFromHandle(handle);
      lastSignatureRef.current = fileSignature(nextRun.files);
      runIdRef.current += 1;
      setLocalRun({ ...nextRun, runId: runIdRef.current });
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

  const handleStaticLiveDataNotice = () => {
    setBrowseStatus({
      kind: 'info',
      text: 'Live Data works in Chromium browsers — Chrome, Edge, Arc, or Opera. You can still open a folder here to view results, or',
      link: {
        href: REPO_URL,
        label: 'install the local app'
      }
    });
  };

  const renderBrowseStatus = () => {
    if (!browseStatus) return null;
    return (
      <div className={`browse-status ${browseStatus.kind}`} role="status">
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
        <button
          type="button"
          className="notification-close"
          onClick={() => setBrowseStatus(null)}
          aria-label="Close notification"
          title="Close"
        >
          &times;
        </button>
      </div>
    );
  };

  return (
    <div className="app-container">
      <main className="main-content">
        <header className="app-header">
          <div className="header-primary">
            <div className="brand-row">
              <div className="brand-mark">R</div>
              <div className="brand-copy">
                <h1>
                  RMCprofile
                  <span>Run Monitor</span>
                </h1>
              </div>
            </div>
            <nav className="page-tabs" aria-label="Workspace pages">
              <button
                className={activePage === 'dashboard' ? 'active' : ''}
                onClick={() => handlePageChange('dashboard')}
              >
                Dashboard
              </button>
              <button
                className={activePage === 'structure' ? 'active' : ''}
                onClick={() => handlePageChange('structure')}
              >
                KDE / 3D
              </button>
            </nav>
          </div>
          {fsAccess ? (
            <div className="path-controls">
              <label className="watch-toggle">
                <input
                  type="checkbox"
                  checked={watchFiles}
                  onChange={(event) => setWatchFiles(event.target.checked)}
                />
                <span aria-hidden="true" />
                <b>Live Data</b>
              </label>
              <div className="path-bar local-file-bar">
                <label>Local run</label>
                <div className="selected-run-name">{localRun?.name || 'No folder selected'}</div>
                <button
                  type="button"
                  onClick={handleSelectFolderFsAccess}
                  disabled={localLoading}
                >
                  {localLoading ? 'Reading' : 'Select Folder'}
                </button>
              </div>
            </div>
          ) : staticMode ? (
            <div className="path-controls">
              <button
                type="button"
                className="watch-toggle watch-toggle-button"
                onClick={handleStaticLiveDataNotice}
                aria-pressed="false"
              >
                <span aria-hidden="true" />
                <b>Live Data</b>
              </button>
              <div className="path-bar local-file-bar">
                <label htmlFor="local-run-files">Local run</label>
                <input
                  ref={directoryInputRef}
                  id="local-run-files"
                  className="visually-hidden"
                  type="file"
                  multiple
                  webkitdirectory=""
                  onChange={handleLocalFiles}
                />
                <div className="selected-run-name">{localRun?.name || 'No folder selected'}</div>
                <button
                  type="button"
                  onClick={() => directoryInputRef.current?.click()}
                  disabled={localLoading}
                >
                  {localLoading ? 'Reading' : 'Select Folder'}
                </button>
              </div>
            </div>
          ) : (
            <div className="path-controls">
              <label className="watch-toggle">
                <input
                  type="checkbox"
                  checked={watchFiles}
                  onChange={(event) => setWatchFiles(event.target.checked)}
                />
                <span aria-hidden="true" />
                <b>Live Data</b>
              </label>
              <form className="path-bar" onSubmit={handleDirectorySubmit}>
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
                  className="secondary-button"
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
          {renderBrowseStatus()}
        </header>
        <div className="workspace-pages">
          {visitedPages.dashboard && (
            <div
              className={`workspace-page${activePage === 'dashboard' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'dashboard'}
            >
              <Dashboard directory={currentDirectory} localRun={localRun} watchFiles={watchFiles} />
            </div>
          )}
          {visitedPages.structure && (
            <div
              className={`workspace-page${activePage === 'structure' ? ' is-active' : ' is-hidden'}`}
              aria-hidden={activePage !== 'structure'}
            >
              <StructurePage directory={currentDirectory} localRun={localRun} theme="light" />
            </div>
          )}
        </div>
      </main>
    </div>
  );
}

export default App;
