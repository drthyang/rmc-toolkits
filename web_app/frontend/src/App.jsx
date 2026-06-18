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
  const directoryInputRef = useRef(null);
  const fileInputRef = useRef(null);
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
              <label htmlFor="local-run-files">Local run</label>
              <input
                ref={directoryInputRef}
                id="local-run-folder"
                className="visually-hidden"
                type="file"
                multiple
                webkitdirectory=""
                onChange={handleLocalFiles}
              />
              <input
                ref={fileInputRef}
                id="local-run-files"
                className="visually-hidden"
                type="file"
                multiple
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
              <button
                type="button"
                className="secondary-button"
                onClick={() => fileInputRef.current?.click()}
                disabled={localLoading}
              >
                Files
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
