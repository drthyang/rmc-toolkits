// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// readStogXy's row grammar (0.6.0 audit, stog-b review): the same inputs as
// tests/test_stog_b_readers.py::NumericTokenTests, so both engines are pinned
// to the same rows. JS \d (no u flag) is ASCII-only, and tokens split on the
// ECMAScript \s set — rmc_toolkits.parsers.read_stog_xy mirrors both.

import { describe, expect, it } from 'vitest';
import { readStogXy } from '../workers/autoScale';

const rows = (text) => readStogXy(text).map((column) => Array.from(column));

describe('readStogXy numeric-token grammar', () => {
  it('reads Fortran D exponents and skips 1_0 / 0x10 rows', () => {
    const text = '3\ntitle\n 0.5D+00 1.1D+00\n 0.51d0 1.25E0\n 5.2E-01 1.3\n 1_0 2\n 0x10 3\n';
    expect(rows(text)).toEqual([[0.5, 0.51, 0.52], [1.1, 1.25, 1.3]]);
  });

  it('treats only ASCII digits as numeric', () => {
    const text = ' 0.5 1.0\n 0.6 1.1\n ٠.7 1.2\n 0.8 ١.5\n 0.9 １.0\n';
    expect(rows(text)).toEqual([[0.5, 0.6], [1.0, 1.1]]);
  });

  it('splits tokens on the ECMAScript whitespace set', () => {
    const text = ' 0.5 1.0\n 0.6 1.1\n 0.7 1.2\n 0.8﻿1.3\n'
      + ' 0.9\u001f1.4\n 1.0\u00851.5\n 1.1\u001c1.6\n';
    expect(rows(text)).toEqual([[0.5, 0.6, 0.7, 0.8], [1.0, 1.1, 1.2, 1.3]]);
  });
});
