// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import cx from './cx';

/**
 * Workspace page shell: the scrolling, padded area below the app header.
 *
 * @param {string}  [as='section'] - element to render.
 * @param {boolean} [column]   - flex column (fill-height pages).
 * @param {boolean} [mobile=true] - tighter padding at ≤760px.
 * @param {boolean} [wide]     - tighter padding on ≥1500px 16:10 screens.
 * @param {boolean} [pbSm]     - smaller bottom padding.
 * @param {boolean} [focusAll] - focus ring on every button on the page.
 */
const Page = ({ as = 'section', column, mobile = true, wide, pbSm, focusAll, className, ...rest }) => {
    const Tag = as;
    return (
        <Tag
            className={cx(
                'ui-page',
                column && 'ui-page--column',
                mobile && 'ui-page--mobile',
                wide && 'ui-page--wide',
                pbSm && 'ui-page--pb-sm',
                focusAll && 'ui-page--focus-all',
                className
            )}
            {...rest}
        />
    );
};

export default Page;
