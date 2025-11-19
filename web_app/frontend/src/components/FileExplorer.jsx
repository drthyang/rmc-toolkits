import React, { useState, useEffect } from 'react';
import axios from 'axios';
import './FileExplorer.css';

const FileExplorer = ({ onFileSelect }) => {
    const [currentPath, setCurrentPath] = useState('./');
    const [files, setFiles] = useState([]);
    const [error, setError] = useState(null);

    useEffect(() => {
        fetchFiles(currentPath);
    }, [currentPath]);

    const fetchFiles = async (path) => {
        try {
            const response = await axios.get(`http://localhost:5000/api/files?dir=${path}`);
            setFiles(response.data.files);
            setError(null);
        } catch (err) {
            setError('Failed to load files');
            console.error(err);
        }
    };

    const handleNavigate = (path, type) => {
        if (type === 'directory') {
            setCurrentPath(path);
        } else {
            onFileSelect(path);
        }
    };

    const handleUp = () => {
        const parent = currentPath.split('/').slice(0, -1).join('/') || '/';
        setCurrentPath(parent);
    };

    return (
        <div className="file-explorer">
            <h3>File Explorer</h3>
            <div className="path-header">
                <button onClick={handleUp} disabled={currentPath === '/'}>⬆ Up</button>
                <span className="current-path">{currentPath}</span>
            </div>
            {error && <div className="error">{error}</div>}
            <ul className="file-list">
                {files.map((file, index) => (
                    <li key={index}
                        className={`file-item ${file.type}`}
                        onClick={() => handleNavigate(file.path, file.type)}>
                        <span className="icon">{file.type === 'directory' ? '📁' : '📄'}</span>
                        {file.name}
                    </li>
                ))}
            </ul>
        </div>
    );
};

export default FileExplorer;
