import React, { useEffect, useState } from 'react';
import axios from 'axios';
import Dashboard from './components/Dashboard';
import StructurePage from './components/StructurePage';
import API_BASE_URL from './api';
import './App.css';

function App() {
  const [activePage, setActivePage] = useState('dashboard');
  const [currentDirectory, setCurrentDirectory] = useState('data');
  const [draftDirectory, setDraftDirectory] = useState('data');
  const [theme, setTheme] = useState(() => localStorage.getItem('rmc-theme') || 'dark');
  const [browserOpen, setBrowserOpen] = useState(false);
  const [browserPath, setBrowserPath] = useState('data');
  const [browserRoot, setBrowserRoot] = useState(null);
  const [browserFiles, setBrowserFiles] = useState([]);
  const [browserLoading, setBrowserLoading] = useState(false);
  const [browserError, setBrowserError] = useState(null);

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
    localStorage.setItem('rmc-theme', theme);
  }, [theme]);

  const handleDirectorySubmit = (event) => {
    event.preventDefault();
    const nextDirectory = draftDirectory.trim() || '.';
    setCurrentDirectory(nextDirectory);
  };

  const loadBrowserPath = async (path, commit = false) => {
    setBrowserLoading(true);
    setBrowserError(null);
    try {
      const response = await axios.get(`${API_BASE_URL}/api/files`, {
        params: { dir: path || '.' }
      });
      const nextPath = response.data.currentPath;
      setBrowserPath(nextPath);
      setDraftDirectory(nextPath);
      setBrowserRoot(response.data.root);
      setBrowserFiles(response.data.files || []);
      if (commit) {
        setCurrentDirectory(nextPath);
        setBrowserOpen(false);
      }
    } catch (err) {
      setBrowserError(err.response?.data?.error || 'Failed to load folder');
    } finally {
      setBrowserLoading(false);
    }
  };

  const openBrowser = () => {
    setBrowserOpen((value) => !value);
    if (!browserOpen) {
      loadBrowserPath(draftDirectory || currentDirectory);
    }
  };

  const handleBrowserUp = () => {
    const cleanPath = browserPath.replace(/\/$/, '');
    const parent = cleanPath.split('/').slice(0, -1).join('/') || '/';
    loadBrowserPath(parent);
  };

  const directories = browserFiles.filter((file) => file.type === 'directory');
  const recognizedFiles = browserFiles.filter((file) => file.type !== 'directory');

  return (
    <div className="app-container">
      <main className="main-content">
        <header className="app-header">
          <div className="header-primary">
            <div className="brand-row">
              <div className="brand-mark">R</div>
              <div className="brand-copy">
                <h1>RMCProfile Run Monitor</h1>
                <span>{currentDirectory}</span>
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
            <div className="header-actions">
              <button
                className="theme-toggle"
                type="button"
                onClick={() => setTheme((value) => (value === 'dark' ? 'light' : 'dark'))}
                aria-label={`Switch to ${theme === 'dark' ? 'light' : 'dark'} theme`}
                title={`Switch to ${theme === 'dark' ? 'light' : 'dark'} theme`}
              >
                <span>{theme === 'dark' ? 'Light' : 'Dark'}</span>
              </button>
            </div>
          </div>
          <form className="path-bar" onSubmit={handleDirectorySubmit}>
            <label htmlFor="data-path">Data path</label>
            <input
              id="data-path"
              type="text"
              value={draftDirectory}
              onChange={(event) => setDraftDirectory(event.target.value)}
              spellCheck="false"
            />
            <button type="button" className="secondary-button" onClick={openBrowser}>
              Browse
            </button>
            <button type="submit">
              Load
            </button>
          </form>
          {browserOpen && (
            <div className="path-browser" role="region" aria-label="Data folder browser">
              <div className="path-browser-header">
                <div>
                  <span>Folder</span>
                  <strong>{browserPath}</strong>
                </div>
                <div className="path-browser-actions">
                  <button type="button" onClick={handleBrowserUp} disabled={!browserRoot || browserPath === browserRoot}>
                    Up
                  </button>
                  <button type="button" onClick={() => loadBrowserPath(browserPath, true)}>
                    Use this folder
                  </button>
                </div>
              </div>
              {browserError && <div className="path-browser-error">{browserError}</div>}
              <div className="path-browser-grid">
                <div>
                  <h2>Folders</h2>
                  <div className="browser-list">
                    {directories.map((file) => (
                      <button type="button" key={file.path} onClick={() => loadBrowserPath(file.path)}>
                        <span>DIR</span>
                        <strong>{file.name}</strong>
                      </button>
                    ))}
                    {!browserLoading && directories.length === 0 && <p>No folders here.</p>}
                  </div>
                </div>
                <div>
                  <h2>Detected files</h2>
                  <div className="browser-list file-preview-list">
                    {recognizedFiles.slice(0, 12).map((file) => (
                      <div key={file.path}>
                        <span>{file.plotKind || 'FILE'}</span>
                        <strong>{file.name}</strong>
                      </div>
                    ))}
                    {!browserLoading && recognizedFiles.length === 0 && <p>No supported files in this folder.</p>}
                  </div>
                </div>
              </div>
            </div>
          )}
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} />}
        {activePage === 'structure' && <StructurePage directory={currentDirectory} theme={theme} />}
      </main>
    </div>
  );
}

export default App;
