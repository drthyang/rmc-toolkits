// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// rho0 self-consistency: the worker and the page report WHY estimateRho0 did not
// converge — a trial density the auto-scale could not fit (estimate.stopped) is
// not "the amplitudes disagree at every density" (scaling_cli.py parity).
// The page shows a one-line summary with the full message behind its help.

import { describe, expect, it } from 'vitest';
import {
  rho0NonConvergenceError,
  rho0NonConvergenceMessage,
  rho0NonConvergenceSummary,
} from '../workers/autoScale';

describe('rho0NonConvergenceMessage', () => {
  it('names the trial density that could not be fitted', () => {
    const stopped = 'autoscale failed at rho0 = 0.66: autoscale: could not locate the first shell';
    const message = rho0NonConvergenceMessage({ concordance: 9.87, stopped });
    expect(message).toContain(stopped);
    expect(message).not.toContain('disagree at every density');
    expect(message).toContain('9.87');
  });

  it('keeps the discordance explanation when the iteration ran out', () => {
    const message = rho0NonConvergenceMessage({ concordance: 0.4321, stopped: null });
    expect(message).toContain('disagree at every density');
    expect(message).toContain('0.432');
  });

  it('has a one-line summary, carried beside the unchanged message', () => {
    const estimate = { concordance: 0.4321, stopped: null };
    expect(rho0NonConvergenceSummary(estimate))
      .toBe('ρ₀ did not converge (concordance 0.432) — set ρ₀, or use the FZ amplitude.');
    const error = rho0NonConvergenceError(estimate);
    expect(error).toBeInstanceOf(Error);
    expect(error.message).toBe(rho0NonConvergenceMessage(estimate));
    expect(error.summary).toBe(rho0NonConvergenceSummary(estimate));
  });
});
