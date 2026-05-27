import React, { useEffect, useState } from 'react';
import Dashboard from './components/Dashboard';
import StructurePage from './components/StructurePage';
import './App.css';

function App() {
  const [activePage, setActivePage] = useState('dashboard');
  const [currentDirectory, setCurrentDirectory] = useState('data');
  const [draftDirectory, setDraftDirectory] = useState('data');
  const [theme, setTheme] = useState(() => localStorage.getItem('rmc-theme') || 'dark');

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
    localStorage.setItem('rmc-theme', theme);
  }, [theme]);

  const handleDirectorySubmit = (event) => {
    event.preventDefault();
    const nextDirectory = draftDirectory.trim() || '.';
    setCurrentDirectory(nextDirectory);
  };

  return (
    <div className="app-container">
      <main className="main-content">
        <header className="app-header">
          <div className="header-primary">
            <div className="brand-row">
              <div className="brand-mark">R</div>
              <div className="brand-copy">
                <h1>RMC++ Workspace</h1>
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
            <button type="submit">
              Load
            </button>
          </form>
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} />}
        {activePage === 'structure' && <StructurePage directory={currentDirectory} theme={theme} />}
      </main>
    </div>
  );
}

export default App;
