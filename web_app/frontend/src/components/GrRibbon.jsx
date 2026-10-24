// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Dashboard's G(r) error ribbon, drawn: the σ toggle for the card header
// and the caption row under the plot, whose ? help (and title) is the tooltip
// with the transform, the reproduction residuals and their tolerance. What
// they say comes from ../grRibbon.js.

import React from 'react';
import { ribbonCaption, ribbonChoice } from '../grRibbon';
import { CardNote, Segmented, SegmentedButton } from '../ui';
import InfoBadge from '../ui/InfoBadge';

/** The σ toggle for the card header: only when the file's σ and the propagated one both exist. */
export const RibbonToggle = ({ sources, choice, onChoose }) => {
    if (!sources?.fileSigma || !sources?.propagated) return null;
    const which = ribbonChoice(sources, choice);
    return (
        <Segmented role="group" aria-label="Error ribbon">
            <SegmentedButton type="button" active={which === 'file'} aria-pressed={which === 'file'} onClick={() => onChoose('file')}
                title="The G(r) file's own σ column">
                σ file
            </SegmentedButton>
            <SegmentedButton type="button" active={which === 'propagated'} aria-pressed={which === 'propagated'} onClick={() => onChoose('propagated')}
                title="σ propagated from F(Q)'s σ through RMCProfile's transform">
                σ propagated
            </SegmentedButton>
        </Segmented>
    );
};

/** The caption row under an xPDF plot: what the ribbon is, or why there is none. */
export const RibbonNote = ({ sources, choice }) => {
    if (!sources) return null;
    const { which, text, details } = ribbonCaption(sources, choice);
    return (
        <CardNote className="gr-ribbon-note" data-ribbon={which} title={details.join('\n') || undefined}>
            {text}
            {details.length > 0 && (
                <>
                    {' '}
                    <InfoBadge label="About the error ribbon" align="end" side="above">
                        {details.map((line) => <span key={line} className="gr-ribbon-note__line">{line}</span>)}
                    </InfoBadge>
                </>
            )}
        </CardNote>
    );
};
