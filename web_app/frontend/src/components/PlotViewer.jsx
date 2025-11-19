import React, { useState, useEffect } from 'react';
import axios from 'axios';
import './PlotViewer.css';

const PlotViewer = ({ filePath }) => {
    const [imageUrl, setImageUrl] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);

    useEffect(() => {
        if (filePath) {
            fetchPlot(filePath);
        }
    }, [filePath]);

    const fetchPlot = async (path) => {
        setLoading(true);
        setError(null);
        setImageUrl(null);
        try {
            const response = await axios.get(`http://localhost:5000/api/plot?path=${encodeURIComponent(path)}`, {
                responseType: 'blob'
            });
            const url = URL.createObjectURL(response.data);
            setImageUrl(url);
        } catch (err) {
            setError('Failed to generate plot. Ensure the file is a valid RMC output.');
            console.error(err);
        } finally {
            setLoading(false);
        }
    };

    if (!filePath) {
        return <div className="plot-viewer empty">Select a file to view plot</div>;
    }

    return (
        <div className="plot-viewer">
            {loading && <div className="loading">Generating Plot...</div>}
            {error && <div className="error">{error}</div>}
            {imageUrl && (
                <div className="plot-container">
                    <img src={imageUrl} alt="RMC Plot" />
                </div>
            )}
        </div>
    );
};

export default PlotViewer;
