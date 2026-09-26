// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';

// The app footer shared by the workspace pages. The look is the kit's
// `ui-footer` class; `tight` trims its top margin on the fill-height pages.
const AppFooter = ({ tight = false }) => (
    <footer className={tight ? 'ui-footer ui-footer--tight' : 'ui-footer'}>
        &copy; 2026 Tsung-Han Yang &middot;{' '}
        <a
            href="https://github.com/drthyang/rmc-toolkits/blob/main/LICENSE"
            target="_blank"
            rel="noreferrer"
        >
            AGPLv3
        </a>
        {' '}&middot;{' '}
        <a
            href="https://github.com/drthyang/rmc-toolkits#readme"
            target="_blank"
            rel="noreferrer"
        >
            About & documentation
        </a>
    </footer>
);

export default AppFooter;
