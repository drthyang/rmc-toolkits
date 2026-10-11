// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Where the Copilot sits and how big: the docked width keeps the pages at
// least PAGE_MIN wide, the sheet stays between its limits, the remembered
// state is validated, and the keyboard and pointer resize land on the limits.

import { describe, expect, it } from 'vitest';
import {
    DOCK_MAX, DOCK_MIN, DOCK_QUERY, PAGE_MIN, SHEET_DEFAULT, SHEET_MAX, SHEET_MIN,
    clampDockWidth, clampSheet, defaultDockWidth, dragDockWidth, dragSheet, keyResize, maxDockWidth,
    readDockState, writeDockState
} from '../components/dockLayout';

const memoryStorage = () => {
    const data = new Map();
    return {
        getItem: (key) => (data.has(key) ? data.get(key) : null),
        setItem: (key, value) => data.set(key, String(value))
    };
};

describe('dock layout', () => {
    it('docks from 1080 × 541 up: the narrowest panel beside pages of 760px', () => {
        expect(DOCK_QUERY).toBe('(min-width: 1080px) and (min-height: 541px)');
    });

    it('defaults to 28 % of the window within 23–30rem', () => {
        expect(defaultDockWidth(1440)).toBe(403);
        expect(defaultDockWidth(1920)).toBe(450);
        expect(defaultDockWidth(1200)).toBe(345);
        // On a 4K monitor the root type size is 21px.
        expect(defaultDockWidth(3840, 21)).toBe(630);
    });

    it('keeps the pages at least PAGE_MIN wide, and the dock within its limits', () => {
        expect(PAGE_MIN).toBe(760);
        expect(maxDockWidth(1440)).toBe(680);
        expect(maxDockWidth(1200)).toBe(440);
        expect(maxDockWidth(2560)).toBe(DOCK_MAX);
        expect(maxDockWidth(900)).toBe(DOCK_MIN);
        expect(clampDockWidth(1000, 1440)).toBe(680);
        expect(clampDockWidth(100, 1440)).toBe(DOCK_MIN);
        expect(clampDockWidth(Number.NaN, 1440)).toBe(403);
        expect(clampDockWidth(410.6, 1440)).toBe(411);
    });

    it('keeps the sheet between its limits', () => {
        expect(clampSheet(0.9)).toBe(SHEET_MAX);
        expect(clampSheet(0.1)).toBe(SHEET_MIN);
        expect(clampSheet(undefined)).toBe(SHEET_DEFAULT);
        expect(clampSheet(0.5)).toBe(0.5);
    });

    it('reads back what it wrote, and ignores anything malformed', () => {
        const storage = memoryStorage();
        expect(readDockState(storage)).toEqual({ open: false, width: null, sheet: SHEET_DEFAULT });
        writeDockState(storage, { open: true, width: 512.4 });
        writeDockState(storage, { sheet: 0.6 });
        expect(readDockState(storage)).toEqual({ open: true, width: 512, sheet: 0.6 });
        storage.setItem('rmc-llm-dock-v1', '{not json');
        expect(readDockState(storage)).toEqual({ open: false, width: null, sheet: SHEET_DEFAULT });
        storage.setItem('rmc-llm-dock-v1', JSON.stringify({ open: 'yes', width: -3, sheet: 7 }));
        expect(readDockState(storage)).toEqual({ open: false, width: null, sheet: SHEET_MAX });
        expect(readDockState(null)).toEqual({ open: false, width: null, sheet: SHEET_DEFAULT });
        const failing = { getItem: () => null, setItem: () => { throw new Error('quota'); } };
        expect(writeDockState(failing, { open: true })).toBeNull();
    });

    it('resizes from the keyboard: ← widens the dock, ↑ raises the sheet, Home/End go to the limits', () => {
        expect(keyResize('docked', 'ArrowLeft', 400, 1440)).toBe(424);
        expect(keyResize('docked', 'ArrowRight', 400, 1440)).toBe(376);
        expect(keyResize('docked', 'Home', 400, 1440)).toBe(DOCK_MIN);
        expect(keyResize('docked', 'End', 400, 1440)).toBe(680);
        expect(keyResize('docked', 'ArrowUp', 400, 1440)).toBeNull();
        expect(keyResize('sheet', 'ArrowUp', 0.46, 1024)).toBe(0.51);
        expect(keyResize('sheet', 'ArrowDown', 0.3, 1024)).toBe(SHEET_MIN);
        expect(keyResize('sheet', 'End', 0.3, 1024)).toBe(SHEET_MAX);
        expect(keyResize('sheet', 'Tab', 0.3, 1024)).toBeNull();
    });

    it('follows the pointer: the dock grows leftward, the sheet\'s top tracks it', () => {
        expect(dragDockWidth(400, 1000, 900, 1440)).toBe(500);
        expect(dragDockWidth(400, 1000, 1200, 1440)).toBe(DOCK_MIN);
        expect(dragSheet(512, 1024)).toBe(0.5);
        expect(dragSheet(10, 1024)).toBe(SHEET_MAX);
    });
});

describe('a sheet that makes room for a focused card', () => {
    it('lowers the sheet until the card fits above it, within its limits, never raising it', async () => {
        const { sheetForCard } = await import('../components/dockLayout');
        // iPad landscape: 820 tall, the pages' header 125 px, a 420 px card.
        expect(sheetForCard({ viewportHeight: 820, topOffset: 125, cardHeight: 420, current: 0.46 })).toBe(0.316);
        // A small card: the sheet stays.
        expect(sheetForCard({ viewportHeight: 820, topOffset: 125, cardHeight: 200, current: 0.46 })).toBe(0.46);
        // A tall sheet comes down to half the window at least.
        expect(sheetForCard({ viewportHeight: 820, topOffset: 125, cardHeight: 100, current: 0.7 })).toBe(0.5);
        // A card taller than the room: the lowest sheet, the card's top shows.
        expect(sheetForCard({ viewportHeight: 820, topOffset: 125, cardHeight: 900, current: 0.46 })).toBe(SHEET_MIN);
        expect(sheetForCard({ viewportHeight: 0, cardHeight: 10, current: 0.6 })).toBe(0.5);
    });
});
