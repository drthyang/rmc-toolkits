import React, { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];

const Dashboard = ({ directory }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);

    useEffect(() => {
        const fetchDashboard = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/files`, {
                    params: { dir: directory || '.' }
                });
                const loadedFiles = response.data.files || [];
                setFiles(loadedFiles);
                const plotFiles = loadedFiles.filter((file) => file.plotKind);
                const metadataEntries = await Promise.all(
                    plotFiles.map(async (file) => {
                        try {
                            const meta = await axios.get(`${API_BASE_URL}/api/plot/metadata`, {
                                params: { path: file.path }
                            });
                            return [file.path, meta.data];
                        } catch {
                            return [file.path, null];
                        }
                    })
                );
                setMetadata(Object.fromEntries(metadataEntries));
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to load dashboard data');
            } finally {
                setLoading(false);
            }
        };

        fetchDashboard();
    }, [directory]);

    const plotFiles = useMemo(() => {
        return files
            .filter((file) => file.plotKind)
            .sort((a, b) => plotOrder.indexOf(a.plotKind) - plotOrder.indexOf(b.plotKind));
    }, [files]);

    const counts = useMemo(() => {
        return files.reduce((acc, file) => {
            const key = file.type === 'directory' ? 'directories' : file.plotKind || 'other';
            acc[key] = (acc[key] || 0) + 1;
            return acc;
        }, {});
    }, [files]);

    return (
        <section className="dashboard-page">
            <div className="dashboard-toolbar">
                <div>
                    <h2>Run Dashboard</h2>
                    <p>{directory}</p>
                </div>
                {loading && <span className="status-pill">Loading</span>}
            </div>

            {error && <div className="dashboard-error">{error}</div>}

            <div className="metric-strip">
                <div className="metric-tile">
                    <span>Plots</span>
                    <strong>{plotFiles.length}</strong>
                </div>
                <div className="metric-tile">
                    <span>RMC Files</span>
                    <strong>{counts.other || 0}</strong>
                </div>
                <div className="metric-tile">
                    <span>Folders</span>
                    <strong>{counts.directories || 0}</strong>
                </div>
                <div className="metric-tile">
                    <span>Rwp Values</span>
                    <strong>{Object.values(metadata).filter((meta) => meta?.metrics?.rwp !== undefined).length}</strong>
                </div>
            </div>

            <div className="plot-grid">
                {plotFiles.map((file) => {
                    const meta = metadata[file.path];
                    return (
                        <article className="plot-card" key={file.path}>
                            <div className="plot-card-header">
                                <div>
                                    <h3>{meta?.title || file.name}</h3>
                                    <p>{file.plotKind}</p>
                                </div>
                                {meta?.metrics?.rwp !== undefined && (
                                    <span className="rwp-chip">Rwp {Number(meta.metrics.rwp).toPrecision(4)}</span>
                                )}
                            </div>
                            <img
                                src={`${API_BASE_URL}/api/plot?path=${encodeURIComponent(file.path)}`}
                                alt={file.name}
                                loading="lazy"
                            />
                        </article>
                    );
                })}
            </div>

            {!loading && plotFiles.length === 0 && (
                <div className="empty-state">Open a run folder such as data to populate the dashboard.</div>
            )}
        </section>
    );
};

export default Dashboard;
