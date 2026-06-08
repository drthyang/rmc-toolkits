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
  const [browseStatus, setBrowseStatus] = useState(null);

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
    localStorage.setItem('rmc-theme', theme);
  }, [theme]);

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
            <button type="button" className="secondary-button" onClick={handleNativeBrowse}>
              Browse
            </button>
            <button type="submit">
              Load
            </button>
          </form>
          {browseStatus && <div className={`browse-status ${browseStatus.kind}`}>{browseStatus.text}</div>}
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} />}
        {activePage === 'structure' && <StructurePage directory={currentDirectory} theme={theme} />}
      </main>
    </div>
  );
}

export default App;
