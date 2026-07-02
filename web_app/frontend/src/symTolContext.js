import { createContext } from 'react';

// Shared symmetry-detection tolerance so the "Detected SG" ladder selection is
// kept when switching between the Dashboard and KDE/3D pages (both keep their
// ModelSummary mounted). Value is a [symTol, setSymTol] pair (useState shape);
// null when no provider is present, so ModelSummary can fall back to local state.
export const SymTolContext = createContext(null);
