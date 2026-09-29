# UI kit (`src/ui`)

The one shared look for every workspace page: a stylesheet of `ui-` classes
([`ui.css`](ui.css)) plus a handful of thin React components, built only on the
design tokens in [`../index.css`](../index.css).

- `ui.css` is imported **once**, from `main.jsx`, directly after `index.css`.
  It therefore loads before every component and page stylesheet.
- Kit `.jsx` files import **no CSS**. Import components from `src/ui`
  (`import { Card, CardHeader } from '../ui'`) or `InfoBadge` / `SaveMenu`
  from their own files.
- Page CSS **places** things; the kit **styles** them. A page may key layout
  off a kit class (`.geom-layout .ui-card { display: flex }`) or wire data to a
  kit token (`.ui-file-chip.kind-bragg { --kind-color: … }`), but never changes
  a kit class's look.

The kit was introduced as a zero-visual-change refactor (every page was moved
onto it pixel for pixel). Near-duplicate looks that existed at that point were
kept as separate variants on purpose — see [Not unified yet](#not-unified-yet).

## Inventory

### Components

| Component | Renders | Props (besides `className`, `...rest`) | Use when |
|---|---|---|---|
| `Page` | `<section class="ui-page …">` | `as`, `column`, `mobile` (default `true`), `wide`, `pbSm`, `focusAll` | the scrolling root of a workspace page |
| `Card` | `<div class="ui-card …">` | `as`, `clip`, `lift`, `roundEnds`, `pad` (`'plot'`/`'bar'`), `note` | any panel surface |
| `CardHeader` | `<h3 class="ui-card__header">` label + meta/actions | `as` (`'div'` for a non-heading bar), `wrap`, `fixed`, slots `title`, `help`, `meta`, `actions` — or `children` | the bar header of a card |
| `CardTitle` | `<h3 class="ui-card__title">` | `as` | title inside a flush header |
| `CardMeta` | `<span class="ui-card__meta">` | `fixed` | a small count / run-meta readout in a header |
| `CardNote` | `<div class="ui-card__note">` | `emph` | a note row under a card's canvas |
| `ControlsBar` | `<div class="ui-controls …">` | `variant` (`'default'`/`'dense'`/`'stacked'`), `sub`, `footer` | the controls bar above a page's cards |
| `ControlGroup` | `<div class="ui-control-group" role="group">` | `label` (→ `aria-label`) | related controls that wrap as a unit |
| `Control` | `<label class="ui-control">` micro-label + widget + value | `as` (`'div'` when the row holds two interactive elements), `label`, `value`, `valueWide` | one labeled control row — inside a `ControlsBar` only |
| `Switch` | `<label class="ui-control ui-switch">` … checkbox + track | `label`, `bare`, `checked`, `onChange`, `inputProps` | a boolean option — inside a `ControlsBar` only |
| `Segmented` | `<div class="ui-seg ui-seg--{variant}">` | `as` (`'nav'`), `variant` (`'frame'`/`'overlay'`/`'nav'`) | a set of mutually exclusive buttons |
| `SegmentedButton` | `<button>` (adds no `type`) | `active`, `overlay`, `warm` | one segment |
| `Pill` | `<button type="button" class="ui-pill">` | `size` (`'md'`), `tint`, `active` | small pill buttons (Show/Hide, Reset zoom) |
| `ToolButton` | `<button type="button" class="ui-tool-btn">` | `axes`, `active` | card-header tools (Reset view, a b c) |
| `IconButton` | `<button type="button" class="ui-icon-btn ui-icon-btn--{variant}">` | `variant` (`'close'`/`'remove'`) | round × buttons |
| `PrimaryButton` | `<button class="ui-btn-primary">` (caller passes `type`) | `outlined` | the one primary action of a bar |
| `Chip` | `<span class="ui-chip">` | `tone` (`'success'`/`'warn'`/`'danger'`), `strong`, `center`, `truncate` | small read-only pills (Rwp, file info) |
| `StatRail` | `<section class="ui-card ui-stat-rail">` title cell + `<dl>` | `heading`, `headingProps` | Model information / Detected SG / Triplet result |
| `Stat` | `<div class="ui-stat"><dt/><dd/></div>` | `label`, `end`, `dtProps`, `ddProps` | one column of a stat rail |
| `StatCard` | readout tile with a status edge | `tone` (`'good'`/`'warn'`/`'bad'`), `label`, `value`, `sub` | result readouts (Auto StoG) |
| `Banner` | `<div class="ui-banner ui-banner--{tone}">` | `as`, `tone` (`'danger'`/`'neutral'`/`'caution'`/`'danger-light'`), `sm`, `gapLg`, `flush`, `inline`, `onDismiss` | messages above or inside a card |
| `Hint` | `<p class="ui-hint">` | | what to do next (dashed box) |
| `EmptyState` | `<div class="ui-empty">` | `fill` | nothing to show yet |
| `InfoBadge` | `?` trigger + popover (`ui-info`) | `label`, `align` (`'start'`/`'end'`), `children` | a short explanation beside a label |
| `SaveMenu` | save trigger + format menu (`ui-save`, `ui-menu`) | `onSave`, `options`, `label`, `align`, `disabled`, `busy`, `className` (`'ui-save--accent'`) | figure export |

Every component appends `className` to its kit classes and spreads `...rest` on
its root, so `role`, `aria-*`, `title`, `style`, `data-*` and `ref` pass through.
Two exceptions take only the props listed: `InfoBadge` (no `className`, no
`...rest`) and `SaveMenu` (`className`, but no `...rest` — its root holds the
ref that closes the menu on an outside click).

`Control` and `Switch` are styled only inside a `ControlsBar`: their rules are
keyed `.ui-controls .ui-control…` (and `.ui-controls .ui-control.ui-switch…`),
so outside one the label loses its micro-label look and a `Switch` shows the
native checkbox beside its track. Making them standalone would change
specificity ties, so it waits for a visible-change pass.

### Class families (used directly where markup varies too much for a component)

| Family | Classes / variants | Use when |
|---|---|---|
| Page | `ui-page` `--column` `--mobile` `--wide` `--pb-sm` `--focus-all` | page root |
| Footer | `ui-footer` `--tight` | the app footer (content in `components/AppFooter.jsx`) |
| Card | `ui-card` `--clip` `--round-ends` `--lift` `--pad-plot` `--pad-bar` `--note` | surfaces |
| Bar header | `ui-card__header` `--wrap` `--fixed`; `ui-card__label`, `__actions`, `__cluster`, `__meta` (`--fixed`), `__readout` | card title bars |
| Flush header | `ui-card__header-flush` `--padded`; `ui-card__heading`, `__title`, `__source` (`__source-item`: one name in a list), `__subtitle`, `__header-actions` | chart cards whose plot continues below the title |
| Inset header | `ui-card__header-inset` (styles its `h3` and `span`) | compact plot cards (Auto StoG) |
| Card notes | `ui-card__note` (`--emph`), `ui-card__caption`, `ui-card__section` | rows and dividers inside a card |
| Legend | `ui-legend`, `__item`, `__swatch`, `__note`, `__warning`, `__group`, `__subitem`, `__credit` | a color key under a canvas |
| Stage | `ui-stage` `--glow` `--glow-soft` `--glow-faint` `--divided` `--orbit` | canvas / WebGL backgrounds (`--orbit` gives the grab cursor) |
| Overlays | `ui-overlay-badge` (`.is-error`), `ui-overlay-controls` `--left` `--right` | status and controls over a canvas |
| Table | `ui-table` `--labels` `--strong-heads` `--abbr`; `td.is-highlight`, `tr.is-dim`, `__center`, `__quiet`, `__note`, `__dot`; `ui-table-block`, `-title`, `-scroll` | dense numeric tables |
| Tag | `ui-tag` (`.is-clean`, `.is-flagged`), `__count` | a status callout row |
| Controls | `ui-controls` `--dense` `--stacked` `--sub` `--footer`; `ui-control-group`; `ui-cluster` `--grow` `--end`, `ui-cluster-label`; `ui-control`, `ui-control-label`, `ui-control-value` (`--wide`); `ui-color-dot` | controls bars |
| Form widgets | `select.ui-select` (`--ring`), `select.ui-select-native`, `input.ui-input`, `input.ui-input-compact`, `ui-pair` + `input.ui-input-strong`, `input.ui-range` (`--lg`), `ui-field` `--formula` `--wide` `--select` | inputs |
| Switches | `ui-switch` (`--bare`, `__track`), `ui-switch-outline` (`--button`), `ui-chip-toggle` (`.is-on`) | boolean options |
| Dropdown | `ui-dropdown`, `__button`, `__list` (`button.is-selected`) | a custom listbox |
| Buttons | `ui-btn-primary` (`--outlined`), `ui-btn-brand` (`.is-active`), `ui-pill`, `ui-pill-tint`, `ui-pill-md` (`.is-active`), `ui-tool-btn` (`--axes`, `.is-active`), `ui-icon-btn` `--close` `--remove` | actions |
| Segmented | `ui-seg` `--frame` `--overlay` `--nav`; `ui-seg__label`, `ui-seg__btn` (`--warm`, `.is-active`); frame/nav buttons take `.is-active` | exclusive choices |
| Chips | `ui-chip` `--strong` `--center` `--truncate` `--success` `--warn` `--danger`; `ui-file-chip` (`.is-hidden`), `__kind`, `__name` | read-only pills |
| Stats | `ui-stack`; `ui-stat-rail`, `__title`, `__source`, `__stats`, `__line`; `ui-stat` (`--end`), `__sub`; `ui-stat-card` (`.is-good/-warn/-bad`), `__label`, `__value`, `__sub`; `ui-inline-stats`, `ui-inline-stat` (`.is-flagged`), `__null` | numbers with labels |
| Feedback | `ui-banner` `--danger` `--neutral` `--caution` `--inline` `--danger-light` `--sm` `--gap-lg` `--flush` `--dismissible`; `ui-status` (`.is-error`); `ui-hint`; `ui-empty` (`--fill`); `ui-placeholder`; `ui-loading` `--sm` `--error` | messages and empty states |
| Floating | `ui-save` (`--accent`), `__trigger`, `__icon`; `ui-menu` `--right` `--left`, `__item`; `ui-info`, `__trigger`, `__popover` `--start` `--end` | menus and popovers |
| Field bar | `ui-fieldbar` (`--readonly`), `__value`, `__ghost` | the app header's run-folder field |
| Forms | `ui-fieldset`, `__fields`; `ui-dropzone` (`.is-drag`), `__hint` | grouped parameters, uploads |
| Utility | `ui-visually-hidden` | hidden but accessible |

## Tokens

All tokens live in `index.css`: the two theme blocks (surfaces, text, accent,
status, shadows, focus ring) and one theme-invariant block (brand, ink on
accent, extra shadows, status and plot-kind colors, chart ink, the select
chevron, radii `--radius-*`, durations `--dur-*`, control heights `--h-*`, and
the type scale `--fs-62` … `--fs-115`, named by hundredths of a rem).

Device tokens: `--viewport-h` (the visible viewport height, `100dvh` where
supported — use it for every height budget instead of `100vh`) and the
safe-area insets `--safe-top` / `--safe-right` / `--safe-bottom` /
`--safe-left` (0 except on notched / rounded screens with
`viewport-fit=cover`). On touch screens (`pointer: coarse`) the `--h-*`
control heights are larger.

The root type size is fluid (15px up to a 1080p-class window, 17px at 1440p,
21px at 2160p; see `index.css`), so every rem grows on 2K / 4K monitors. Size
anything that should grow with the UI in rem — a px size stays small on a 4K
screen.

In `ui.css` every color, shadow, radius, duration, control height and font size
is a `var()`. Literals are allowed only for spacing (padding / gap / margin),
font-weight, letter-spacing, line-height, `color-mix()` percentages, relative
`em` sizes, and the geometry of a single widget (range track and thumbs, switch
knob offsets; rem where it should scale with the UI). A new value becomes a
token first.

Media blocks in `ui.css` adapt the kit to handhelds: at ≤ 760px card headers
wrap their actions under the title and `InfoBadge` popovers open as a
full-width sheet under their trigger, and on touch phones / portrait
tablets (`hover: none` and `pointer: coarse`, ≤ 1040 wide or ≤ 540 tall) the
kit's fields use 16px text, below which iOS Safari zooms the page on focus.

## Rules

- **No page-specific look.** Page CSS may only place things — grid templates
  and areas, flex sizing, `order`, outer margins between blocks, width/height
  clamps, canvas sizing — keyed by page hook classes (`pca-layout`,
  `geom-layout`, `orient-*`, `analysis-layout`, `r-value-card`, …). It may wire
  domain data to tokens (e.g. `.ui-file-chip.kind-bragg { --kind-color: … }`).
  Borders, colors, radii, shadows, typography and control geometry come from
  the kit only, apart from the exceptions listed next.
- **Domain visualizations keep their look in their component CSS**, token-only:
  the InteractivePlot SVG marks, legend and tooltip (`InteractivePlot.css`), the
  ModelSummary tolerance ladder (`ModelSummary.css`), the OrientationView
  colorbar, axis-view panel and tooltip (`OrientationView.css`).
- **The other look rules outside the kit** are these, and only these:
  - the app shell in `App.css`, token-only: the `.app-container` background,
    the `.app-header` bar (bottom border, background) and the brand lockup
    (`.brand-mark` radius, fill, glow and ink; `.brand-copy h1` and its `span`
    type);
  - the statistics-column dividers in `PcaKdePage.css` (`.pca-stats-col`
    `::before` hairline, and its `border-top` when stacked), token-only,
    because they are positioned from that page's grid gap;
  - the llm module's own stylesheets (`src/llm/components/*.css`; the module
    is kept extractable);
  - `FileExplorer.css` and `PlotViewer.css`, whose components nothing renders
    (not moved onto the kit).
- **State classes** (`is-active`, `is-on`, `is-hidden`, `is-error`, …) are
  always compounded with a block class, never a bare `.is-hidden {}`. In
  `ui.css` that block is a `ui-` class; outside it a state rule compounds with
  its component's own class — today `.workspace-page.is-hidden` (`App.css`),
  `.sym-brick.is-active` (the ModelSummary ladder); the llm module scopes its
  states to its own `llm-*` classes.
- **Specificity is part of the look.** `index.css` element rules interact with
  the kit, notably `button:hover:not(:disabled)` (0,2,1: panel-raised
  background + strong border) and `button:focus-visible` (0,1,1: the focus
  ring). Several looks depend on who wins those ties — e.g. `.ui-seg--nav
  button` (0,1,1) suppresses the ring because `ui.css` loads after `index.css`,
  and single-class buttons (`ui-pill`, `ui-save__trigger`, `ui-tool-btn`, …)
  still receive the global hover background. So: keep selector shapes when
  editing, and never use `@layer`, `:where()` or `!important` in `ui.css`.
  (The one `!important` in the app, `background: transparent` on the
  InteractivePlot legend's dashed swatch in `InteractivePlot.css`, predates
  the kit.)
- **Order inside `ui.css`**: blocks in catalog order (page, footer, card,
  headers, notes, legend, stage, table, tag, controls, form widgets, switches,
  dropdown, buttons, segmented, chips, stats, feedback, floating, field bar /
  forms / utilities); modifiers after their base; media queries after the
  base. One deliberate exception: `ui-table--labels` precedes the base table
  rules (its label column takes the base `muted` color).
- **`src/llm` uses kit class strings only**, never kit JS — it must stay
  extractable (see `src/llm/README.md`).

## Adding a variant

1. Reuse an existing variant if it is pixel-identical.
2. Otherwise add `ui-<block>--<name>` directly after the block's base — or a
   separate class when the base's hover/active rules would repaint it (that is
   why `ui-pill`, `ui-pill-tint` and `ui-pill-md` are three classes).
3. New values become tokens in `index.css`.
4. Add it to the inventory above, with the page that uses it.
5. Never override a kit look from page CSS.

## Adding a component

Only where markup repeats. Keep it thin: `className` appended, `...rest` spread
on the root, `as` when the element varies, slots or `children`. A component
must not add roles, aria attributes, `type` or wrapper elements that its call
sites lacked (callers that differ use the children form or the classes). Add a
case to `__tests__/markup.test.jsx`.

## Verifying a change

- `npm test` (includes the kit markup tests), `npx eslint src`, `npm run build`.
- Layout changes: load the Demo run and walk every page at the browser content
  sizes of the target devices (window chrome subtracted): iPhone 17 Pro
  402 × 874 and landscape 874 × 402, iPhone 17 Pro Max 440 × 956, iPad Air 11"
  820 × 1110, iPad Pro 13" 1032 × 1310 and 1376 × 970, MacBook Air 13"
  1470 × 840, MacBook Pro 14" 1512 × 870, MacBook Pro 16" 1728 × 1000, a
  1366 × 650 laptop, 1080p 1920 × 960, 1440p 2560 × 1310 and 4K 3840 × 2030
  (emulate touch for the phones and iPads). Check that the header keeps its
  tier, no page scrolls sideways, and the fill-height pages end at the footer.
- Compare computed styles before/after on every page and breakpoint (the tag +
  index path of each element → its computed style, `::before`/`::after`
  included) — a class rename with identical styling must produce an identical
  dump.
- Screenshot diffs of headless Chromium are not perfectly deterministic here: the page loads
  Inter from Google Fonts (`display=swap`), and text advances occasionally come out 1/64 px
  different between runs of the *same* build, which flips the sub-pixel phase of a few glyphs
  (≈0.06 % of a phone screenshot, a handful of glyphs). Re-shoot before chasing such a diff; a
  real CSS change reproduces.
- A second noise class does reproduce, so re-shooting alone does not rule it out: the
  anti-aliasing of a rounded corner can flip between two values depending on paint history.
  Seen on the brand mark's top corners in the header (≈6 pixels, max channel delta ≈46, e.g.
  rgb(211,223,248) vs rgb(165,190,246)) after the Live Data switch was toggled on and off; it
  persisted through a whole browser context and repeated across runs of each build, while
  geometry and computed styles were identical, and a later run of `main` rendered the other
  value. Before treating a few-pixel corner diff as a regression, shoot `main` against itself
  (A vs A′) in the same scenario — if `main` flips too, it is noise.
- Walk the states by hand (DevTools :hover / :focus-visible / :active):
  nav tabs (hover, active-tab hover, no focus ring), pills and the save trigger
  (incl. accent), info trigger, tool buttons, overlay and frame segments,
  the dropdown trigger, legend pills (ring only on Auto StoG), the primary
  buttons, × buttons (incl. a hidden file chip), the Demo button in both
  states, the field bar `:focus-within`, the header switch (label and button
  forms), control switches, range thumbs, card lift, InfoBadge popovers.

## Not unified yet

Kept as separate variants when the kit was introduced, to change nothing on
screen; candidates for a later, visible cleanup: the micro-label sizes
(.62–.68rem, .05–.09em tracking), the five pill looks, three primary buttons,
three card-header styles, the banner looks, styled vs native selects and the
UA-styled Auto StoG export input, frame vs overlay vs nav segments, solid vs
outline switch vs chip toggle, the accent save trigger losing its tint on hover,
the active tab going transparent on hover, the dropdown chevron vanishing on
hover, legend pills without a focus ring (except on Auto StoG), the `--fs-*`
scale, chart ink that ignores the theme, and the legacy hook names (`pca-*`,
`orient-*`, `geom-*`).
