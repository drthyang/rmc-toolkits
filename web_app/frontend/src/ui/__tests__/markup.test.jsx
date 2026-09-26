// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The kit components must render exactly the markup of the call sites they
// replaced: the same elements, attributes and roles — a component never adds a
// role, an aria attribute, a `type` or a wrapper the call site did not have.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import {
    Banner, Card, CardHeader, CardMeta, CardNote, CardTitle, Chip, Control, ControlGroup, ControlsBar,
    EmptyState, Hint, IconButton, Page, Pill, PrimaryButton, Segmented, SegmentedButton, Stat, StatCard,
    StatRail, Switch, ToolButton,
} from '..';

const html = (node) => renderToStaticMarkup(node);

describe('Page', () => {
    it('renders a section with the page classes (mobile padding by default)', () => {
        expect(html(<Page>x</Page>)).toBe('<section class="ui-page ui-page--mobile">x</section>');
    });
    it('maps the flags and passes the rest through', () => {
        expect(html(<Page as="div" column wide pbSm focusAll mobile={false} className="hook" data-x="1" />))
            .toBe('<div class="ui-page ui-page--column ui-page--wide ui-page--pb-sm ui-page--focus-all hook" data-x="1"></div>');
    });
});

describe('Card', () => {
    it('renders the surface with its modifiers', () => {
        expect(html(<Card roundEnds className="pca-viewport" />))
            .toBe('<div class="ui-card ui-card--round-ends pca-viewport"></div>');
        expect(html(<Card as="article" clip lift data-figure-card="" />))
            .toBe('<article class="ui-card ui-card--clip ui-card--lift" data-figure-card=""></article>');
        expect(html(<Card as="details" pad="bar" />)).toBe('<details class="ui-card ui-card--pad-bar"></details>');
        expect(html(<Card clip note />)).toBe('<div class="ui-card ui-card--clip ui-card--note"></div>');
    });

    it('CardHeader slots render label + actions inside an h3', () => {
        expect(html(<CardHeader title="T" help={<i />} actions={<b />} />))
            .toBe('<h3 class="ui-card__header"><span class="ui-card__label">T<i></i></span><span class="ui-card__actions"><b></b></span></h3>');
    });

    it('CardHeader renders a meta node as a direct child and no actions span when absent', () => {
        expect(html(<CardHeader fixed title="P" meta={<CardMeta>m</CardMeta>} />))
            .toBe('<h3 class="ui-card__header ui-card__header--fixed"><span class="ui-card__label">P</span><span class="ui-card__meta">m</span></h3>');
        expect(html(<CardHeader title="A" />)).toBe('<h3 class="ui-card__header"><span class="ui-card__label">A</span></h3>');
    });

    it('CardHeader children form keeps the caller markup (div for a non-heading bar)', () => {
        expect(html(<CardHeader as="div" wrap><span>S</span></CardHeader>))
            .toBe('<div class="ui-card__header ui-card__header--wrap"><span>S</span></div>');
    });

    it('CardTitle, CardMeta and CardNote', () => {
        expect(html(<CardTitle>t</CardTitle>)).toBe('<h3 class="ui-card__title">t</h3>');
        expect(html(<CardMeta fixed>m</CardMeta>)).toBe('<span class="ui-card__meta ui-card__meta--fixed">m</span>');
        expect(html(<CardNote emph role="status">n</CardNote>))
            .toBe('<div class="ui-card__note ui-card__note--emph" role="status">n</div>');
    });
});

describe('Controls', () => {
    it('ControlsBar variants', () => {
        expect(html(<ControlsBar />)).toBe('<div class="ui-controls"></div>');
        expect(html(<ControlsBar variant="dense" />)).toBe('<div class="ui-controls ui-controls--dense"></div>');
        expect(html(<ControlsBar variant="stacked" sub />)).toBe('<div class="ui-controls ui-controls--stacked ui-controls--sub"></div>');
        expect(html(<ControlsBar variant="stacked" footer />)).toBe('<div class="ui-controls ui-controls--stacked ui-controls--footer"></div>');
    });

    it('ControlGroup is a labelled group', () => {
        expect(html(<ControlGroup label="Site and sampling" />))
            .toBe('<div class="ui-control-group" role="group" aria-label="Site and sampling"></div>');
    });

    it('Control renders label, widget and value in order', () => {
        expect(html(<Control label="Box" value="3.0σ"><input type="range" /></Control>))
            .toBe('<label class="ui-control"><span class="ui-control-label">Box</span><input type="range"/><span class="ui-control-value">3.0σ</span></label>');
        expect(html(<Control label="Slice" value="0.50" valueWide />))
            .toContain('<span class="ui-control-value ui-control-value--wide">0.50</span>');
        expect(html(<Control as="div" label="Shell" />)).toBe('<div class="ui-control"><span class="ui-control-label">Shell</span></div>');
    });

    it('Switch: label, checkbox, aria-hidden track; bare has no label', () => {
        expect(html(<Switch label="Wireframe" checked={false} onChange={() => {}} inputProps={{ 'aria-label': 'Show ellipsoid wireframe' }} />))
            .toBe('<label class="ui-control ui-switch"><span class="ui-control-label">Wireframe</span><input type="checkbox" aria-label="Show ellipsoid wireframe"/><i class="ui-switch__track" aria-hidden="true"></i></label>');
        expect(html(<Switch bare checked onChange={() => {}} />))
            .toBe('<label class="ui-control ui-switch ui-switch--bare"><input type="checkbox" checked=""/><i class="ui-switch__track" aria-hidden="true"></i></label>');
    });
});

describe('Segmented', () => {
    it('frame / overlay / nav containers', () => {
        expect(html(<Segmented role="group" aria-label="Reference frame" />))
            .toBe('<div class="ui-seg ui-seg--frame" role="group" aria-label="Reference frame"></div>');
        expect(html(<Segmented as="nav" variant="nav" className="page-tabs" aria-label="Workspace pages" />))
            .toBe('<nav class="ui-seg ui-seg--nav page-tabs" aria-label="Workspace pages"></nav>');
    });

    it('SegmentedButton adds no type and only the state class', () => {
        expect(html(<SegmentedButton active>PC</SegmentedButton>)).toBe('<button class="is-active">PC</button>');
        expect(html(<SegmentedButton>Dashboard</SegmentedButton>)).toBe('<button>Dashboard</button>');
        expect(html(<SegmentedButton type="button" overlay warm active>a</SegmentedButton>))
            .toBe('<button class="ui-seg__btn ui-seg__btn--warm is-active" type="button">a</button>');
    });
});

describe('Buttons', () => {
    it('Pill looks', () => {
        expect(html(<Pill aria-expanded>Hide</Pill>)).toBe('<button type="button" class="ui-pill" aria-expanded="true">Hide</button>');
        expect(html(<Pill tint>Reset zoom</Pill>)).toBe('<button type="button" class="ui-pill-tint">Reset zoom</button>');
        expect(html(<Pill size="md" active>Advanced</Pill>)).toBe('<button type="button" class="ui-pill-md is-active">Advanced</button>');
    });

    it('ToolButton, IconButton, PrimaryButton', () => {
        expect(html(<ToolButton axes active aria-pressed>abc</ToolButton>))
            .toBe('<button type="button" class="ui-tool-btn ui-tool-btn--axes is-active" aria-pressed="true">abc</button>');
        expect(html(<IconButton variant="remove" aria-label="Hide x">×</IconButton>))
            .toBe('<button type="button" class="ui-icon-btn ui-icon-btn--remove" aria-label="Hide x">×</button>');
        expect(html(<PrimaryButton type="button" className="geom-compute">Compute</PrimaryButton>))
            .toBe('<button class="ui-btn-primary geom-compute" type="button">Compute</button>');
        expect(html(<PrimaryButton type="button" outlined>Go</PrimaryButton>))
            .toBe('<button class="ui-btn-primary ui-btn-primary--outlined" type="button">Go</button>');
    });
});

describe('Chip', () => {
    it('tones and flags', () => {
        expect(html(<Chip tone="success" strong>Rwp 0.1</Chip>)).toBe('<span class="ui-chip ui-chip--strong ui-chip--success">Rwp 0.1</span>');
        expect(html(<Chip center truncate title="t">f</Chip>)).toBe('<span class="ui-chip ui-chip--center ui-chip--truncate" title="t">f</span>');
    });
});

describe('Stats', () => {
    it('StatRail is a card section with a title cell and a dl', () => {
        expect(html(<StatRail aria-label="Model information" heading="Model information" headingProps={{ title: 'f.rmc6f' }}><Stat label="Cell">1</Stat></StatRail>))
            .toBe('<section class="ui-card ui-stat-rail" aria-label="Model information"><h2 class="ui-stat-rail__title" title="f.rmc6f">Model information</h2>'
                + '<dl class="ui-stat-rail__stats"><div class="ui-stat"><dt>Cell</dt><dd>1</dd></div></dl></section>');
    });

    it('Stat passes role, dt and dd attributes through', () => {
        expect(html(<Stat end role="status" dtProps={{ className: 'sym-ladder-dt' }} ddProps={{ title: 'x' }} label="L">v</Stat>))
            .toBe('<div class="ui-stat ui-stat--end" role="status"><dt class="sym-ladder-dt">L</dt><dd title="x">v</dd></div>');
    });

    it('StatCard renders label, value, sub', () => {
        expect(html(<StatCard tone="warn" label="L" value="V" sub="S" />))
            .toBe('<div class="ui-stat-card is-warn"><span class="ui-stat-card__label">L</span><span class="ui-stat-card__value">V</span><span class="ui-stat-card__sub">S</span></div>');
    });
});

describe('Feedback', () => {
    it('Banner keeps the caller role; onDismiss adds the close button', () => {
        expect(html(<Banner as="p" tone="caution" role="status">w</Banner>))
            .toBe('<p class="ui-banner ui-banner--caution" role="status">w</p>');
        expect(html(<Banner tone="danger" flush role="alert" onDismiss={() => {}}>boom</Banner>))
            .toBe('<div class="ui-banner ui-banner--danger ui-banner--flush ui-banner--dismissible" role="alert"><span>boom</span>'
                + '<button type="button" class="ui-icon-btn ui-icon-btn--close" aria-label="Close notification" title="Close">×</button></div>');
        expect(html(<Banner tone="danger" gapLg>e</Banner>)).toBe('<div class="ui-banner ui-banner--danger ui-banner--gap-lg">e</div>');
        expect(html(<Banner as="p" tone="danger" sm>e</Banner>)).toBe('<p class="ui-banner ui-banner--danger ui-banner--sm">e</p>');
        expect(html(<Banner tone="danger-light" inline>e</Banner>)).toBe('<div class="ui-banner ui-banner--danger-light ui-banner--inline">e</div>');
    });

    it('Hint and EmptyState', () => {
        expect(html(<Hint>h</Hint>)).toBe('<p class="ui-hint">h</p>');
        expect(html(<EmptyState fill>e</EmptyState>)).toBe('<div class="ui-empty ui-empty--fill">e</div>');
    });
});
