import React, { useState } from 'react';
import FileExplorer from './components/FileExplorer';
import PlotViewer from './components/PlotViewer';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);

  return (
    <div className="app-container">
      <FileExplorer onFileSelect={setSelectedFile} />
      <main className="main-content">
        <header className="app-header">
          <h1>RMC++ Plotter</h1>
          {selectedFile && <span className="selected-file">{selectedFile}</span>}
        </header>
        <PlotViewer filePath={selectedFile} />
      </main>
    </div>
  );
}

export default App;
