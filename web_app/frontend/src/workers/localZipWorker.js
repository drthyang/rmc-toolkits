import { strFromU8, unzipSync } from 'fflate';
import { buildLocalRunFromEntries } from '../browserData';

self.onmessage = (event) => {
    const { id, buffer, name } = event.data;
    try {
        const unzipped = unzipSync(new Uint8Array(buffer));
        const entries = Object.entries(unzipped)
            .filter(([path]) => !path.endsWith('/'))
            .map(([path, bytes]) => ({
                name: path.split('/').pop() || path,
                path,
                text: strFromU8(bytes)
            }));
        const run = buildLocalRunFromEntries(entries);
        self.postMessage({
            id,
            result: {
                ...run,
                name: run.name === 'Local files' ? name.replace(/\.zip$/i, '') : run.name
            }
        });
    } catch (error) {
        self.postMessage({
            id,
            error: error.message || 'Could not import ZIP file'
        });
    }
};
