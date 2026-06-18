import { structureFromRmc6f } from '../browserData';

self.onmessage = async (event) => {
    const { id, file, maxPoints } = event.data;
    try {
        if (!file?.sourceFile) {
            throw new Error('No browser structure file available');
        }
        const text = await file.sourceFile.text();
        const result = structureFromRmc6f({ ...file, text }, maxPoints);
        self.postMessage({ id, result });
    } catch (error) {
        self.postMessage({
            id,
            error: error.message || 'Browser structure parser failed'
        });
    }
};
