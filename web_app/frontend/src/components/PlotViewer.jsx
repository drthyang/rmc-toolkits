import React, { useState, useEffect } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import './PlotViewer.css';

const PlotViewer = ({ file, onConverted }) => {
    const [imageUrl, setImageUrl] = useState(null);
    const [metadata, setMetadata] = useState(null);
    const [loading, setLoading] = useState(false);
    const [converting, setConverting] = useState(false);
    const [error, setError] = useState(null);
    const [message, setMessage] = useState(null);
    const filePath = file?.path;

    useEffect(() => {
        if (filePath && file.plotKind) {
            fetchPlot(filePath);
        } else {
            setImageUrl(null);
            setMetadata(null);
            setError(null);
            setMessage(null);
        }
    }, [filePath, file?.plotKind]);

    const fetchPlot = async (path) => {
        setLoading(true);
        setError(null);
        setMessage(null);
        setImageUrl(null);
        setMetadata(null);
        try {
            const [plotResponse, metadataResponse] = await Promise.all([
                axios.get(`${API_BASE_URL}/api/plot`, {
                    params: { path },
                    responseType: 'blob'
                }),
                axios.get(`${API_BASE_URL}/api/plot/metadata`, {
                    params: { path }
                })
            ]);
            const url = URL.createObjectURL(plotResponse.data);
            setImageUrl(url);
            setMetadata(metadataResponse.data);
        } catch (err) {
            setError(err.response?.data?.error || 'Failed to generate plot. Ensure the file is a valid RMC output.');
            console.error(err);
        } finally {
            setLoading(false);
        }
    };

    const convertRmc6f = async () => {
        setConverting(true);
        setError(null);
        setMessage(null);
        try {
            const response = await axios.post(`${API_BASE_URL}/api/convert/frac`, {
                path: filePath,
                overwrite: true
            });
            setMessage(`Generated ${response.data.name}`);
            onConverted?.();
        } catch (err) {
            setError(err.response?.data?.error || 'Failed to generate Frac_coord file.');
            console.error(err);
        } finally {
            setConverting(false);
        }
    };

    useEffect(() => {
        return () => {
            if (imageUrl) {
                URL.revokeObjectURL(imageUrl);
            }
        };
    }, [imageUrl]);

    if (!filePath) {
        return <div className="plot-viewer empty">Select a file</div>;
    }

    return (
        <div className="plot-viewer">
            <div className="viewer-toolbar">
                <div>
                    <h2>{file.name}</h2>
                    <p>{file.path}</p>
                </div>
                {filePath?.endsWith('.rmc6f') && (
                    <button className="conversion-button" onClick={convertRmc6f} disabled={converting}>
                        {converting ? 'Generating...' : 'Generate Frac_coord'}
                    </button>
                )}
            </div>
            <div className="viewer-body">
                {loading && <div className="loading">Generating plot...</div>}
                {message && <div className="message">{message}</div>}
                {error && <div className="error">{error}</div>}
                {metadata && (
                    <div className="plot-meta">
                        <span>{metadata.title}</span>
                        {Object.entries(metadata.metrics || {}).map(([key, value]) => (
                            <span key={key}>{key}: {Number(value).toPrecision(5)}</span>
                        ))}
                    </div>
                )}
                {imageUrl && (
                    <div className="plot-container">
                        <img src={imageUrl} alt="RMC Plot" />
                    </div>
                )}
                {!loading && !imageUrl && !error && !filePath?.endsWith('.rmc6f') && (
                    <div className="empty-panel">No plot preview is available for this file.</div>
                )}
                {!loading && !imageUrl && filePath?.endsWith('.rmc6f') && (
                    <div className="empty-panel">Structure files can be converted for KDE and 3D workflows.</div>
                )}
            </div>
        </div>
    );
};

export default PlotViewer;
