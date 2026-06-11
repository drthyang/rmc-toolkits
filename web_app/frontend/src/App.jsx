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
  const [browseStatus, setBrowseStatus] = useState(null);

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
          {browseStatus && <div className={`browse-status ${browseStatus.kind}`}>{browseStatus.text}</div>}
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} />}
        {activePage === 'structure' && <StructurePage directory={currentDirectory} theme="light" />}
        <footer className="app-footer">
          &copy; 2026 Tsung-Han Yang. All rights reserved.
        </footer>
      </main>
    </div>
  );
}

export default App;
