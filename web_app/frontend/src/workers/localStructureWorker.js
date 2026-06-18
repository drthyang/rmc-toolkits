import { structureFromRmc6f } from '../browserData';

self.onmessage = (event) => {
    const { id, file, maxPoints } = event.data;
    try {
        if (!file?.text) {
            throw new Error('No structure file data available');
        }
        self.postMessage({
            id,
            result: structureFromRmc6f(file, maxPoints)
        });
    } catch (error) {
        self.postMessage({
            id,
            error: error.message || 'Browser structure parser failed'
        });
    }
};
