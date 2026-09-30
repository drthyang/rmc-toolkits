// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The kit components must render exactly the markup of the call sites they
// replaced: the same elements, attributes and roles — a component never adds a
// role, an aria attribute, a `type` or a wrapper the call site did not have.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import {
    Banner, BondDash, Card, CardHeader, CardMeta, CardNote, CardTitle, Chip, Control, ControlGroup, ControlsBar,
    ElementChip, EmptyState, Hint, IconButton, Kpi, KpiRail, Page, Pill, PrimaryButton, SaveMenu, Segmented,
    SegmentedButton, Stat, StatCard, StatRail, Switch, ToolButton, UnitField,
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

    it('ControlsBar as a form', () => {
        expect(html(<ControlsBar as="form" className="geom-controls" noValidate />))
            .toBe('<form class="ui-controls geom-controls" noValidate=""></form>');
    });

    it('UnitField: inputProps on the input, the unit inside the border, invalid state', () => {
        expect(html(<UnitField unit="Å" inputProps={{ type: 'number', value: '2.00', readOnly: true, 'aria-label': 'Min' }} />))
            .toBe('<span class="ui-unit-field"><input class="ui-unit-field__input" type="number" readOnly="" aria-label="Min" value="2.00"/>'
                + '<span class="ui-unit-field__unit">Å</span></span>');
        expect(html(<UnitField unit="°" invalid className="x" data-f="1" />))
            .toBe('<span class="ui-unit-field is-invalid x" data-f="1"><input class="ui-unit-field__input"/><span class="ui-unit-field__unit">°</span></span>');
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

    it('PrimaryButton run: states only with run, the label stays the text', () => {
        expect(html(<PrimaryButton type="submit" run>Compute</PrimaryButton>))
            .toBe('<button class="ui-btn-primary ui-btn-primary--run" type="submit">Compute</button>');
        expect(html(<PrimaryButton type="submit" run busy stale aria-busy>Update</PrimaryButton>))
            .toBe('<button class="ui-btn-primary ui-btn-primary--run is-busy is-stale" type="submit" aria-busy="true">Update</button>');
        expect(html(<PrimaryButton busy stale>Go</PrimaryButton>)).toBe('<button class="ui-btn-primary">Go</button>');
    });
});

describe('Chip', () => {
    it('tones and flags', () => {
        expect(html(<Chip tone="success" strong>Rwp 0.1</Chip>)).toBe('<span class="ui-chip ui-chip--strong ui-chip--success">Rwp 0.1</span>');
        expect(html(<Chip center truncate title="t">f</Chip>)).toBe('<span class="ui-chip ui-chip--center ui-chip--truncate" title="t">f</span>');
    });

    it('ElementChip: element color on --chip, a dot, the central ring', () => {
        expect(html(<ElementChip color="#00A087">Se</ElementChip>))
            .toBe('<span class="ui-element-chip" style="--chip:#00A087"><i class="ui-element-chip__dot" aria-hidden="true"></i>Se</span>');
        expect(html(<ElementChip color="#F39B7F" central title="central atom">Ta</ElementChip>))
            .toBe('<span class="ui-element-chip ui-element-chip--central" style="--chip:#F39B7F" title="central atom">'
                + '<i class="ui-element-chip__dot" aria-hidden="true"></i>Ta</span>');
    });

    it('BondDash: bond color on --bond and hidden text, so a chain reads as text', () => {
        expect(html(<BondDash color="#1f6fd6" />))
            .toBe('<span class="ui-bond-dash" style="--bond:#1f6fd6"><span class="ui-visually-hidden">–</span></span>');
        expect(html(<BondDash lead text="" />)).toBe('<span class="ui-bond-dash ui-bond-dash--lead"><span class="ui-visually-hidden"></span></span>');
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

    it('KpiRail wraps its dl; Kpi shows value, unit and sub, or "—" and a kept sub line', () => {
        // A live region's role goes on the wrapper: on the <dl> itself it would
        // replace the description-list semantics (ARIA in HTML allows a dl only
        // group, list, none and presentation).
        expect(html(<KpiRail role="status" aria-live="polite" aria-label="R"><Kpi label="Angles" value="14.8" unit="per Ta" sub="236,431 angles" title="t" /></KpiRail>))
            .toBe('<div class="ui-kpis" role="status" aria-live="polite" aria-label="R"><dl class="ui-kpis__list">'
                + '<div class="ui-kpi" title="t"><dt>Angles</dt><dd>'
                + '<span class="ui-kpi__value">14.8<span class="ui-kpi__unit">\u2009per Ta</span></span>'
                + '<span class="ui-kpi__sub">236,431 angles</span></dd></div></dl></div>');
        expect(html(<Kpi label="Coordination" value={null} unit="per Ta" sub="x" />))
            .toBe('<div class="ui-kpi"><dt>Coordination</dt><dd><span class="ui-kpi__value is-empty">—</span>'
                + '<span class="ui-kpi__sub">\u00a0</span></dd></div>');
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

describe('SaveMenu', () => {
    it('renders the trigger (single format: no menu semantics)', () => {
        expect(html(<SaveMenu onSave={() => {}} />))
            .toBe('<div class="ui-save"><button type="button" class="ui-save__trigger" title="Save figure">'
                + '<span class="ui-save__icon" aria-hidden="true">⤓</span>Save</button></div>');
    });

    it('announces a menu for several formats and keeps the accent passthrough', () => {
        const markup = html(<SaveMenu onSave={() => {}} className="ui-save--accent" label="Save all figures"
            options={[{ id: 'png', label: 'PNG' }, { id: 'svg', label: 'SVG' }]} />);
        expect(markup).toContain('<div class="ui-save ui-save--accent">');
        expect(markup).toContain('aria-haspopup="menu"');
        expect(markup).toContain('aria-expanded="false"');
        expect(markup).toContain('Save all figures');
    });

    it('shows the busy label and disables the trigger', () => {
        expect(html(<SaveMenu onSave={() => {}} busy />)).toContain('class="ui-save__trigger" disabled=""');
        expect(html(<SaveMenu onSave={() => {}} busy />)).toContain('Saving…');
    });
});
