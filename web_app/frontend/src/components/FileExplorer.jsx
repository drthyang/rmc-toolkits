import React, { useCallback, useState, useEffect } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import './FileExplorer.css';

const FileExplorer = ({ onFileSelect, onDirectoryChange, refreshKey }) => {
    const [currentPath, setCurrentPath] = useState('.');
    const [draftPath, setDraftPath] = useState('.');
    const [rootPath, setRootPath] = useState(null);
    const [files, setFiles] = useState([]);
    const [error, setError] = useState(null);

    const fetchFiles = useCallback(async (path) => {
        try {
            const response = await axios.get(`${API_BASE_URL}/api/files`, {
                params: { dir: path }
            });
            setFiles(response.data.files);
            setCurrentPath(response.data.currentPath);
            setDraftPath(response.data.currentPath);
            setRootPath(response.data.root);
            onDirectoryChange?.(response.data.currentPath);
            setError(null);
        } catch (err) {
            setError(err.response?.data?.error || 'Failed to load files');
            console.error(err);
        }
    }, [onDirectoryChange]);

    useEffect(() => {
        // File listing is the external state this component synchronizes with.
        fetchFiles(currentPath);
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [refreshKey]);

    const handleNavigate = (path, type, file) => {
        if (type === 'directory') {
            fetchFiles(path);
        } else {
            onFileSelect(file);
        }
    };

    const handleUp = () => {
        const cleanPath = currentPath.replace(/\/$/, '');
        const parent = cleanPath.split('/').slice(0, -1).join('/') || '/';
        fetchFiles(parent);
    };

    return (
        <div className="file-explorer">
            <h3>File Explorer</h3>
            <div className="path-header">
                <button onClick={handleUp} disabled={currentPath === rootPath}>Up</button>
                <div className="path-input-container">
                    <input
                        type="text"
                        value={draftPath}
                        onChange={(e) => setDraftPath(e.target.value)}
                        onKeyDown={(e) => e.key === 'Enter' && fetchFiles(draftPath)}
                        className="path-input"
                    />
                    <button onClick={() => fetchFiles(draftPath)}>Go</button>
                </div>
            </div>
            {error && <div className="error">{error}</div>}
            <ul className="file-list">
                {files.map((file, index) => (
                    <li key={index}
                        className={`file-item ${file.type}`}
                        onClick={() => handleNavigate(file.path, file.type, file)}>
                        <span className="icon">{file.type === 'directory' ? 'DIR' : file.plotKind || 'FILE'}</span>
                        {file.name}
                    </li>
                ))}
            </ul>
        </div>
    );
};

export default FileExplorer;
