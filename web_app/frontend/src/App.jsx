import React, { useState } from 'react';
import Dashboard from './components/Dashboard';
import FileExplorer from './components/FileExplorer';
import PlotViewer from './components/PlotViewer';
import StructurePage from './components/StructurePage';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [refreshKey, setRefreshKey] = useState(0);
  const [activePage, setActivePage] = useState('dashboard');
  const [currentDirectory, setCurrentDirectory] = useState('.');
  const handleFileSelect = (file) => {
    setSelectedFile(file);
    setActivePage('file');
  };

  return (
    <div className="app-container">
      <FileExplorer
        onFileSelect={handleFileSelect}
        onDirectoryChange={setCurrentDirectory}
        refreshKey={refreshKey}
      />
      <main className="main-content">
        <header className="app-header">
          <div className="brand-row">
            <h1>RMC++ Workspace</h1>
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
              <button
                className={activePage === 'file' ? 'active' : ''}
                onClick={() => setActivePage('file')}
                disabled={!selectedFile}
              >
                File
              </button>
            </nav>
          </div>
          <span className="selected-file">{selectedFile?.path || currentDirectory}</span>
        </header>
        {activePage === 'dashboard' && <Dashboard directory={currentDirectory} />}
        {activePage === 'structure' && <StructurePage directory={currentDirectory} />}
        {activePage === 'file' && (
          <PlotViewer
            file={selectedFile}
            onConverted={() => setRefreshKey((value) => value + 1)}
          />
        )}
      </main>
    </div>
  );
}

export default App;
