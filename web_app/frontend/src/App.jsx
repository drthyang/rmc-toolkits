import React, { useEffect, useRef, useState } from 'react';
import axios from 'axios';
import Dashboard from './components/Dashboard';
import StructurePage from './components/StructurePage';
import API_BASE_URL from './api';
import { buildLocalRun, isStaticMode } from './browserData';
import './App.css';

function App() {
  const [activePage, setActivePage] = useState('dashboard');
  const [currentDirectory, setCurrentDirectory] = useState('data');
  const [draftDirectory, setDraftDirectory] = useState('data');
  const [browseStatus, setBrowseStatus] = useState(null);
  const [localRun, setLocalRun] = useState(null);
  const [localLoading, setLocalLoading] = useState(false);
  const [localLoadingLabel, setLocalLoadingLabel] = useState('Reading');
  const directoryInputRef = useRef(null);
  const zipInputRef = useRef(null);
  const zipRequestRef = useRef(0);
  const staticMode = isStaticMode();

  useEffect(() => {
    document.documentElement.dataset.theme = 'light';
    localStorage.setItem('rmc-theme', 'light');
  }, []);

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
    setLocalLoadingLabel('Reading');
    setBrowseStatus({ kind: 'loading', text: 'Reading local files...' });
    try {
      const nextRun = await buildLocalRun(selectedFiles);
      setLocalRun(nextRun);
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

  const handleLocalZip = async (event) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;
    const id = zipRequestRef.current + 1;
    zipRequestRef.current = id;
    setLocalLoading(true);
    setLocalLoadingLabel('Importing');
    setBrowseStatus({ kind: 'loading', text: 'Importing ZIP archive...' });
    try {
      const buffer = await selectedFile.arrayBuffer();
      const nextRun = await new Promise((resolve, reject) => {
        const worker = new Worker(new URL('./workers/localZipWorker.js', import.meta.url), {
          type: 'module'
        });
        worker.onmessage = (messageEvent) => {
          worker.terminate();
          if (messageEvent.data.id !== id) return;
          if (messageEvent.data.error) {
            reject(new Error(messageEvent.data.error));
            return;
          }
          resolve(messageEvent.data.result);
        };
        worker.onerror = () => {
          worker.terminate();
          reject(new Error('Browser ZIP importer failed'));
        };
        worker.postMessage({ id, name: selectedFile.name, buffer }, [buffer]);
      });
      setLocalRun(nextRun);
      setCurrentDirectory(nextRun.name);
      setDraftDirectory(nextRun.name);
      setBrowseStatus(null);
      setActivePage('dashboard');
    } catch (error) {
      setBrowseStatus({ kind: 'error', text: error.message || 'Could not import ZIP file' });
    } finally {
      setLocalLoading(false);
      event.target.value = '';
    }
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
                onClick={() => setActivePage('dashboard')}
              >
                Dashboard
              </button>
              <button
                className={activePage === 'structure' ? 'active' : ''}
                onClick={() => setActivePage('structure')}
              >
                KDE / 3D
              </button>
            </nav>
          </div>
          {staticMode ? (
            <div className="path-bar local-file-bar">
              <label htmlFor="local-run-zip">Local run</label>
              <input
                ref={zipInputRef}
                id="local-run-zip"
                className="visually-hidden"
                type="file"
                accept=".zip,application/zip"
                onChange={handleLocalZip}
              />
              <input
                ref={directoryInputRef}
                id="local-run-folder"
                className="visually-hidden"
                type="file"
                multiple
                webkitdirectory=""
                onChange={handleLocalFiles}
              />
              <div className="selected-run-name">{localRun?.name || 'No run selected'}</div>
              <button
                type="button"
                onClick={() => zipInputRef.current?.click()}
                disabled={localLoading}
              >
                {localLoading ? localLoadingLabel : 'Import ZIP'}
              </button>
              <button
                type="button"
                className="secondary-button"
                onClick={() => directoryInputRef.current?.click()}
                disabled={localLoading}
              >
                Folder
              </button>
            </div>
          ) : (
            <form className="path-bar" onSubmit={handleDirectorySubmit}>
              <label htmlFor="data-path">Data path</label>
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
                Browse
              </button>
              <button type="submit">
                Load
              </button>
            </form>
          )}
          {browseStatus && <div className={`browse-status ${browseStatus.kind}`}>{browseStatus.text}</div>}
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} localRun={localRun} />}
        {activePage === 'structure' && (
          <StructurePage directory={currentDirectory} localRun={localRun} theme="light" />
        )}
      </main>
    </div>
  );
}

export default App;
