// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Join the truthy class names: cx('ui-card', clip && 'ui-card--clip', className).
// Returns undefined (no class attribute) when nothing is left.
const cx = (...parts) => {
    const joined = parts.filter(Boolean).join(' ');
    return joined || undefined;
};

export default cx;
