// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The UI kit's React components. The look lives in ./ui.css, which main.jsx
// imports once (right after index.css); the components import no CSS.
// See ./README.md for the inventory and the rules.
export { default as Page } from './Page';
export { Card, CardHeader, CardTitle, CardMeta, CardNote } from './Card';
export { ControlsBar, ControlGroup, Control, Switch, UnitField } from './Controls';
export { Segmented, SegmentedButton } from './Segmented';
export { Pill, ToolButton, IconButton, PrimaryButton } from './Buttons';
export { default as Chip, ElementChip, BondDash } from './Chip';
export { StatRail, Stat, StatCard, KpiRail, Kpi } from './Stats';
export { Banner, Hint, EmptyState } from './Feedback';
export { default as InfoBadge } from './InfoBadge';
export { default as SaveMenu } from './SaveMenu';
export { IssueStoreProvider, IssueScope, IssueList, IssueWatcher, PageIssues } from './Issues';
export { useIssue, useIssueSet, useReportIssue, useResolveIssue } from './issueStore';
