// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { createContext, useCallback, useMemo, useState } from 'react';

// Shared symmetry-detection tolerance so the "Detected SG" ladder selection is
// kept when switching between the Dashboard and KDE/3D pages (both keep their
// ModelSummary mounted). Value is a [symTol, setSymTol] pair (useState shape);
// null when no provider is present, so ModelSummary can fall back to local state.
export const SymTolContext = createContext(null);

/**
 * The shared tolerance of the open run, as a [value, setValue] pair (useState
 * shape; setValue also takes an updater). It starts at `initial` and returns
 * to it whenever `runKey` changes — another run opened, so a pick made for one
 * run never sets the next run's card or the AI Copilot's context — even when
 * that runKey was open before. It is kept across page switches and across a
 * Live Data save or a Load of the same run (the same runKey).
 */
export const useRunTolerance = (runKey, initial) => {
    const [picked, setPicked] = useState(() => ({ runKey, value: initial }));
    // Another run: back to the default in this render, before any page or the
    // Copilot reads the last run's pick.
    if (picked.runKey !== runKey) setPicked({ runKey, value: initial });
    const value = picked.runKey === runKey ? picked.value : initial;
    const setValue = useCallback((next) => setPicked((current) => ({
        runKey: current.runKey,
        value: typeof next === 'function' ? next(current.value) : next
    })), []);
    return useMemo(() => [value, setValue], [value, setValue]);
};
