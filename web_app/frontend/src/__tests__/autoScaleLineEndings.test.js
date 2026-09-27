// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// CR-only (classic Mac) line endings (0.6.0 audit, stog-b): Python reads text with
// universal newlines, so a bare '\r' ends a line; the JS readers split only on
// /\r?\n/ and saw one line ("no numeric rows" / "1 non-empty lines"). All three
// readers now accept LF, CRLF and CR alike.

import { describe, expect, it } from 'vitest';
import { readDatHeader, readStogInp, readStogXy } from '../workers/autoScale';

const INP = [
  '1', 'data.dat', '1.0 28.0', '-9 0.1', '0', 'scale.fq', 'scale.gr', '50', '5000', 'N',
  '0.063049', '0', 'N', 'Y', '1.0', 'scale_ft.sq', 'scale_ft.gr', '0.015407',
  'rmc.fq', 'rmc.gr', 'rmc.dr', '2.48 2.65 3.1',
];
const XY = ['        3', 'Q S(Q)', ' 0.50 1.10 0.01', ' 0.51 1.25 0.01', ' 0.52 1.30 0.01'];
const DAT = ['TITLE :: FeCoSn 199K', 'NUMBER_DENSITY :: 0.057329 Angstrom^(-3)', 'MINIMUM_DISTANCES :: 2.4 2.2', ' 0.5 1.0'];

describe('readers accept LF, CRLF and CR-only line endings', () => {
  ['\n', '\r\n', '\r'].forEach((eol) => {
    const label = JSON.stringify(eol);
    it(`readStogInp with ${label}`, () => {
      const inp = readStogInp(INP.join(eol) + eol);
      expect(inp.dataFile).toBe('data.dat');
      expect(inp.a).toBeCloseTo(10, 12);
      expect(inp.peakRmax).toBeCloseTo(3.1, 12);
    });

    it(`readStogXy with ${label}`, () => {
      const columns = readStogXy(XY.join(eol));
      expect(columns.length).toBe(3);
      expect(Array.from(columns[1])).toEqual([1.1, 1.25, 1.3]);
    });

    it(`readDatHeader with ${label}`, () => {
      const header = readDatHeader(DAT.join(eol));
      expect(header.title).toBe('FeCoSn 199K');
      expect(header.numberDensity).toBeCloseTo(0.057329, 12);
      expect(header.minDistance).toBeCloseTo(2.2, 12);
    });
  });
});
