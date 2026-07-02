const SUPPORTED_NAMES = new Set(['scale_ft.gr', 'scale_ft.sq', 'scale_ft_rmc.fq', 'stog_input.dat']);

export const isStaticMode = () => {
    if (import.meta.env.VITE_STATIC_MODE === 'true') return true;
    if (import.meta.env.VITE_STATIC_MODE === 'false') return false;
    if (import.meta.env.DEV && !import.meta.env.VITE_API_BASE_URL) return true;
    return window.location.hostname.endsWith('github.io');
};

// File System Access API (window.showDirectoryPicker) powers static-mode Live Data.
// Chromium-only (Chrome/Edge/Arc/Opera); Firefox/Safari fall back to <input webkitdirectory>.
export const supportsFileSystemAccess = () => (
    import.meta.env.DEV && !import.meta.env.VITE_API_BASE_URL
        ? false
        : typeof window !== 'undefined' && typeof window.showDirectoryPicker === 'function'
);

export const WATCH_INTERVAL_MS = 3000;

// Stable fingerprint of a run's files; a change signals the dashboard should reload.
export const fileSignature = (items) => items
    .map((file) => `${file.path}:${file.modified ?? ''}:${file.size ?? ''}:${file.plotKind ?? ''}`)
    .sort()
    .join('|');

export const detectPlotKind = (name) => {
    if (/-EXAFS-.+_Q_OUTPUT\.csv$/.test(name)) return 'exafs_q';
    if (/-EXAFS-.+_R_OUTPUT\.csv$/.test(name)) return 'exafs_r';
    if (/_FT_XFQ\d+\.csv$/.test(name)) return 'xpdf';
    if (name.includes('PDF') && name.endsWith('.csv')) {
        return name.includes('PDFpartials') ? 'pdf_partials' : 'npdf';
    }
    if (name.endsWith('_FQ1.csv')) return 'xray_sq';
    if (name.endsWith('_SQ1.csv')) return 'neutron_sq';
    if (name.endsWith('_bragg.csv')) return 'bragg';
    if (/-\d{2,}\.log$/.test(name)) return 'r_value';
    if (['scale_ft.gr', 'scale_ft.sq', 'scale_ft_rmc.fq'].includes(name)) return 'stog';
    return null;
};

const isSupportedFile = (name) => (
    name.endsWith('.csv')
    || name.endsWith('.log')
    || name.endsWith('.rmc6f')
    || name.startsWith('Frac')
    || SUPPORTED_NAMES.has(name)
);

const basename = (path) => path.split('/').pop() || path;
const dirname = (path) => path.includes('/') ? path.split('/').slice(0, -1).join('/') : '';

const runStemFromOutputName = (name) => {
    const patterns = [
        [0, /^(.+)-\d{2,}\.log$/],
        [1, /^(.+)-EXAFS-.+_[QR]_OUTPUT\.csv$/],
        [1, /^(.+)_FT_XFQ\d+\.csv$/],
        [1, /^(.+)_[FS]Q\d+\.csv$/],
        [1, /^(.+)_bragg\.csv$/],
        [1, /^(.+)_PDF(?:partials|\d+)?\.csv$/],
        [2, /^Frac_coord_(.+)\.txt$/]
    ];
    for (const [priority, pattern] of patterns) {
        const match = name.match(pattern);
        if (match) return { priority, stem: match[1] };
    }
    return null;
};

const chooseStructureFile = (files) => {
    const rmc6fFiles = files.filter((file) => file.name.endsWith('.rmc6f'));
    if (!rmc6fFiles.length) return null;
    const rmc6fByLocationAndStem = new Map(
        rmc6fFiles.map((file) => [`${dirname(file.path)}/${file.name.replace(/\.rmc6f$/, '')}`, file])
    );
    const outputStems = files
        .map((file) => {
            const stemMatch = runStemFromOutputName(file.name);
            if (!stemMatch) return null;
            return {
                ...stemMatch,
                directory: dirname(file.path),
                sortName: file.name.toLowerCase()
            };
        })
        .filter(Boolean)
        .sort((a, b) => a.priority - b.priority || a.sortName.localeCompare(b.sortName));

    for (const output of outputStems) {
        const match = rmc6fByLocationAndStem.get(`${output.directory}/${output.stem}`);
        if (match) return match;
    }
    return rmc6fFiles[0];
};

const parseNumberRows = (lines, startIndex = 0, separator = /\s+/) => {
    const rows = [];
    lines.slice(startIndex).forEach((line) => {
        const values = line.trim().split(separator).filter(Boolean);
        if (!values.length) return;
        const parsed = values.map(Number);
        if (parsed.every(Number.isFinite)) rows.push(parsed);
    });
    return rows;
};

const transpose = (rows) => rows[0].map((_, index) => rows.map((row) => row[index]));

const readRmcCsv = (text, name) => {
    const lines = text.split(/\r?\n/).filter((line) => line.trim());
    if (!lines.length) throw new Error(`${name} is empty`);
    const labels = lines[0].split(',').map((label) => label.trim());
    const rows = lines.slice(1).map((line, index) => {
        const values = line.split(',').map((value) => value.trim()).filter(Boolean);
        if (!values.length) return null;
        if (values.length !== labels.length) {
            throw new Error(`${name} line ${index + 2} has ${values.length} values; expected ${labels.length}`);
        }
        return values.map(Number);
    }).filter(Boolean);
    if (!rows.length) throw new Error(`${name} does not contain numeric rows`);
    return { labels, data: transpose(rows) };
};

const csvValues = (line) => line.split(',').map((value) => value.trim()).filter(Boolean);

const numericCsvValues = (line) => {
    const values = csvValues(line);
    if (!values.length) return null;
    const parsed = values.map(Number);
    return parsed.every(Number.isFinite) ? parsed : null;
};

const readExafsCsv = (text, name) => {
    const lines = text.split(/\r?\n/);
    const dataStart = lines.findIndex((line) => numericCsvValues(line));
    if (dataStart <= 0) {
        throw new Error(`${name} does not contain an EXAFS column header and numeric rows`);
    }

    const labels = csvValues(lines[dataStart - 1]);
    const rows = lines.slice(dataStart).map((line, index) => {
        const values = numericCsvValues(line);
        if (!values) return null;
        if (values.length !== labels.length) {
            throw new Error(`${name} line ${dataStart + index + 1} has ${values.length} values; expected ${labels.length}`);
        }
        return values;
    }).filter(Boolean);
    if (!rows.length) throw new Error(`${name} does not contain numeric rows`);
    return { labels, data: transpose(rows) };
};

const readChi = (text) => {
    const chiR = [];
    text.split(/\r?\n/).slice(2).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (parts.length >= 2) {
            const value = Number(parts[parts.length - 1]);
            if (Number.isFinite(value)) chiR.push(value);
        }
    });
    return chiR;
};

const readStog = (text, name) => {
    const rows = parseNumberRows(text.split(/\r?\n/), 2);
    if (!rows.length) throw new Error(`${name} does not contain STOG numeric rows`);
    return transpose(rows);
};

const rwp = (x, observed, fitted) => {
    void x;
    const denom = observed.reduce((sum, value) => sum + value * value, 0);
    if (!denom) return 0;
    const residual = fitted.reduce((sum, value, index) => sum + (value - observed[index]) ** 2, 0);
    return Math.sqrt(residual / denom);
};

const pdfIndex = (name) => {
    const match = name.match(/PDF(\d+)\.csv$/);
    return match ? Number(match[1]) : 0;
};

const cleanAxisLabel = (label) => {
    const normalized = label.trim();
    if (normalized === 'Q') return 'Q (Å^{-1})';
    if (normalized === 'r' || normalized === 'R') return 'r (Å)';
    return normalized
        .replace('(A^-1)', '(Å^{-1})')
        .replace('(A^{-1})', '(Å^{-1})')
        .replace('(A)', '(Å)');
};

// Bragg x-axis from the data's column header: time-of-flight CSVs label the first
// column "Flight time (us)", shown as ToF in microseconds (the conventional neutron
// unit, used as-is). Everything else (e.g. "Q or theta") stays a Q axis.
const braggAxis = (header) => (
    /tof|flight|time/.test((header || '').toLowerCase()) ? 'ToF (µs)' : 'Q (Å^{-1})'
);

export const plotMetadataFromFile = (file) => {
    const kind = file.plotKind;
    if (kind === 'xpdf') return { kind, title: 'xPDF', metrics: file.plotData?.metrics || {} };
    if (kind === 'exafs_q') return { kind, title: 'EXAFS Q-space', metrics: file.plotData?.metrics || {} };
    if (kind === 'exafs_r') return { kind, title: 'EXAFS R-space', metrics: file.plotData?.metrics || {} };
    if (kind === 'npdf') return { kind, title: file.name.replace(/\.[^.]+$/, '').split('_').pop(), metrics: file.plotData?.metrics || {} };
    if (kind === 'pdf_partials') return { kind, title: file.name.replace(/\.[^.]+$/, '').split('_').pop(), metrics: file.plotData?.metrics || {} };
    if (kind === 'xray_sq') return { kind, title: 'S(Q) (x-ray)', metrics: file.plotData?.metrics || {} };
    if (kind === 'neutron_sq') return { kind, title: 'S(Q) (neutron)', metrics: file.plotData?.metrics || {} };
    if (kind === 'bragg') return { kind, title: 'BRAGG', metrics: file.plotData?.metrics || {} };
    if (kind === 'r_value') return { kind, title: 'R-value', metrics: file.plotData?.metrics || {} };
    if (kind === 'stog') return { kind, title: file.name, metrics: file.plotData?.metrics || {} };
    return null;
};

export const plotDataFromText = (file) => {
    const kind = file.plotKind;
    if (!kind) return null;

    if (kind === 'r_value') {
        const yValues = readChi(file.text);
        if (!yValues.length) throw new Error(`${file.name} does not contain chi values`);
        return {
            kind,
            title: 'R-value',
            metrics: { final_chi_r: yValues[yValues.length - 1] },
            xLabel: 'Time steps',
            yLabel: 'log(χ)',
            series: [{
                label: 'R',
                x: yValues.map((_, index) => index),
                y: yValues.map((value) => Math.log(Math.max(value, 1e-12)))
            }]
        };
    }

    if (kind === 'stog') {
        const data = readStog(file.text, file.name);
        return {
            kind,
            title: file.name,
            metrics: {},
            xLabel: file.name.endsWith('.gr') ? 'r (Å)' : 'Q (Å^{-1})',
            yLabel: file.name.endsWith('.gr') ? 'G(r)' : 'S(Q)',
            series: [{ label: file.name, x: data[0], y: data[1] }]
        };
    }

    const csv = ['exafs_q', 'exafs_r'].includes(kind)
        ? readExafsCsv(file.text, file.name)
        : readRmcCsv(file.text, file.name);
    const metrics = {};
    if (['xpdf', 'npdf', 'xray_sq', 'neutron_sq', 'bragg'].includes(kind) && csv.data.length >= 3) {
        metrics.rwp = rwp(csv.data[0], csv.data[1], csv.data[2]);
    }
    if (kind === 'npdf') metrics.pdf_index = pdfIndex(file.name);

    let xLabel = csv.labels[0] || 'x';
    let yLabel = 'data';
    if (kind === 'exafs_q') {
        xLabel = 'k (Å^{-1})';
        yLabel = 'χ(k) k²';
    } else if (kind === 'exafs_r') {
        xLabel = 'r (Å)';
        yLabel = 'FT[χ(k) k²]';
    } else if (['xpdf', 'npdf', 'pdf_partials'].includes(kind)) {
        xLabel = 'r (Å)';
        yLabel = 'G(r)';
    } else if (['xray_sq', 'neutron_sq'].includes(kind)) {
        xLabel = 'Q (Å^{-1})';
        yLabel = 'S(Q)';
    } else if (kind === 'bragg') {
        xLabel = braggAxis(csv.labels[0]);
        yLabel = 'Intensity';
    } else xLabel = cleanAxisLabel(xLabel);

    const title = plotMetadataFromFile({ ...file, plotData: { metrics } })?.title || file.name;
    return {
        kind,
        title,
        metrics,
        xLabel,
        yLabel,
        series: csv.labels.slice(1).map((label, index) => ({
            label: label.trim() || `Series ${index + 1}`,
            x: csv.data[0],
            y: csv.data[index + 1]
        })).filter((series) => series.y)
    };
};

const readCellVectors = (text) => {
    const lines = text.split(/\r?\n/);
    let latticeVectors = null;
    let supercell = null;
    lines.forEach((line, index) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Supercell') supercell = parts.slice(-3).map(Number);
        if (parts[0] === 'Lattice') {
            latticeVectors = [lines[index + 1], lines[index + 2], lines[index + 3]]
                .map((row) => row.trim().split(/\s+/).map(Number));
        }
    });
    if (!latticeVectors || !supercell) throw new Error('Missing lattice or supercell metadata');
    return { latticeVectors, supercell };
};

// Circular mean of an angle-like quantity in [0,1): averages the box copies of a
// site's within-cell fraction so a boundary-wrapping site (≈0 ≡ 1) lands on the
// true position, and thermal displacement in a single snapshot averages out.
const TWO_PI = 2 * Math.PI;

export const structureFromRmc6f = (file, maxPoints = 100) => {
    const { latticeVectors, supercell } = readCellVectors(file.text);
    const counts = {};
    const atomIndices = {};
    const atoms = [];
    const rnAcc = new Map();   // referenceNumber -> { element, sc:[3], ss:[3] } for the circular-mean basis
    let inAtoms = false;
    file.text.split(/\r?\n/).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Atoms:') {
            inAtoms = true;
            return;
        }
        if (!inAtoms || parts.length < 10) return;
        const referenceNumber = Number(parts[6]);
        const element = parts[1];
        const coords = parts.slice(3, 6).map(Number);
        const cellIndices = parts.slice(7, 10).map(Number);
        if (![referenceNumber, ...coords, ...cellIndices].every(Number.isFinite)) return;
        counts[element] = (counts[element] || 0) + 1;
        if (!atomIndices[element]) atomIndices[element] = new Set();
        atomIndices[element].add(referenceNumber);
        atoms.push({ element, referenceNumber, coords, cellIndices });
        // Accumulate this atom's within-cell fraction into its reference-number site.
        let acc = rnAcc.get(referenceNumber);
        if (!acc) { acc = { element, sc: [0, 0, 0], ss: [0, 0, 0] }; rnAcc.set(referenceNumber, acc); }
        for (let i = 0; i < 3; i++) {
            const wf = ((coords[i] * supercell[i]) % 1 + 1) % 1;   // within-unit-cell fraction
            acc.sc[i] += Math.cos(TWO_PI * wf); acc.ss[i] += Math.sin(TWO_PI * wf);
        }
    });

    // One representative site per reference number (its circular-mean fraction) —
    // the (element, fractional) basis the symmetry finder consumes.
    const basis = [...rnAcc.entries()]
        .sort(([a], [b]) => a - b)
        .map(([referenceNumber, acc]) => ({
            el: acc.element,
            referenceNumber,
            frac: [0, 1, 2].map((i) => {
                const a = Math.atan2(acc.ss[i], acc.sc[i]) / TWO_PI;
                return a - Math.floor(a);
            })
        }));

    const stride = Math.max(1, Math.ceil(atoms.length / maxPoints));
    const points = atoms.filter((_, index) => index % stride === 0).slice(0, maxPoints).map((atom) => {
        const reduced = atom.coords.map((value, index) => value - atom.cellIndices[index] / supercell[index]);
        const unitCell = reduced.map((value, index) => ((value * supercell[index]) % 1 + 1) % 1);
        return {
            element: atom.element,
            referenceNumber: atom.referenceNumber,
            boxX: atom.coords[0],
            boxY: atom.coords[1],
            boxZ: atom.coords[2],
            x: unitCell[0],
            y: unitCell[1],
            z: unitCell[2]
        };
    });

    return {
        source: file.path,
        totalAtoms: atoms.length,
        sampledAtoms: points.length,
        sampleStride: stride,
        elements: Object.keys(counts).sort(),
        elementCounts: counts,
        atomIndices: Object.fromEntries(Object.entries(atomIndices).map(([key, value]) => [key, [...value].sort((a, b) => a - b)])),
        supercell,
        latticeVectors,
        basis,
        points
    };
};

export const readAndParseLocalPlotFile = async (file) => {
    if (!file.sourceFile) {
        throw new Error(`${file.name} is not backed by a browser file`);
    }
    const text = await file.sourceFile.text();
    return plotDataFromText({ ...file, text });
};

// Build a run object from { path, file } pairs. Shared by the <input webkitdirectory>
// path (buildLocalRun) and the File System Access path (buildLocalRunFromHandle).
const makeRunFromEntries = (entries) => {
    const files = entries
        .filter(({ path }) => isSupportedFile(basename(path)))
        .map(({ path, file }) => ({
            name: basename(path),
            path,
            type: 'file',
            plotKind: detectPlotKind(basename(path)),
            sourceFile: file,
            size: file.size,
            modified: file.lastModified
        }));
    if (!files.length) {
        throw new Error(`No supported RMCProfile files found in ${entries.length} selected files`);
    }

    const rmc6f = chooseStructureFile(files);
    const directoryRoot = files
        .map((file) => file.path)
        .find((path) => path.includes('/'))
        ?.split('/')[0];

    return {
        name: directoryRoot || 'Local files',
        files,
        structure: null,
        structureFile: rmc6f || null,
        structureError: rmc6f ? 'Structure data loads when needed' : 'No model structure detected',
        diagnostics: {
            selectedFileCount: entries.length,
            supportedFileCount: files.length,
            plotFileCount: files.filter((file) => file.plotKind).length,
            hasStructureFile: Boolean(rmc6f)
        }
    };
};

export const buildLocalRun = async (fileList) => {
    const entries = [...fileList].map((file) => ({
        path: file.webkitRelativePath || file.name,
        file
    }));
    return makeRunFromEntries(entries);
};

// Recursively walk a FileSystemDirectoryHandle into { path, file } pairs. getFile()
// returns a fresh File each call, so re-walking the same handle reflects on-disk changes.
export const collectHandleEntries = async (dirHandle, prefix = '') => {
    const entries = [];
    for await (const [name, entry] of dirHandle.entries()) {
        const path = prefix ? `${prefix}/${name}` : name;
        if (entry.kind === 'file') {
            entries.push({ path, file: await entry.getFile() });
        } else if (entry.kind === 'directory') {
            entries.push(...await collectHandleEntries(entry, path));
        }
    }
    return entries;
};

export const buildLocalRunFromHandle = async (dirHandle) => {
    // Prefix with the folder name so paths match the webkitdirectory layout
    // (root folder as the first segment) and directoryRoot resolves to the folder.
    const entries = await collectHandleEntries(dirHandle, dirHandle.name);
    return makeRunFromEntries(entries);
};
