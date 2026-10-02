// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/cifWriter.js
//
// CIF 1.1 text for the symmetry-averaged structure (averageStructure.js). Formatting only:
// every number was computed upstream. The file is ASCII throughout, as CIF 1.1 requires,
// so comments spell Å, ×, ≥ and δ out. The symmetry operations are always written out. On
// ITA's origin (the usual case) they are ITA's own and the H–M symbol, with ':2' for origin
// choice 2 and ':H' for hexagonal R axes, describes the file by itself; otherwise the origin
// is the .rmc6f model's moved onto the 1/48 grid, and the listed operations are authoritative.

const DIGITS = { xyz: 6, cell: 5, angle: 4, U: 5, occupancy: 4 };

/** Fixed-point text with negative zero suppressed. */
export const fixed = (value, digits) => {
  const text = Number(value).toFixed(digits);
  return /^-0\.?0*$/.test(text) ? text.slice(1) : text;
};

const gcd = (a, b) => (b ? gcd(b, a % b) : a);

/** A translation component as a reduced fraction ('1/2', '3/4'), '' for 0, else decimals. */
export function fractionText(x) {
  const y = x - Math.floor(x);
  if (y < 1e-6 || y > 1 - 1e-6) return '';
  for (let d = 2; d <= 48; d++) {
    const k = Math.round(y * d);
    if (Math.abs(y * d - k) < 1e-6 * d) {
      const g = gcd(k, d);
      return `${k / g}/${d / g}`;
    }
  }
  return fixed(y, DIGITS.xyz);
}

/** A symmetry operation {R | t} as a coordinate triplet: 'x, y, z', '-y+1/2, x-y, z+1/3'. */
export function symopText(R, t) {
  return [0, 1, 2].map((i) => {
    let term = '';
    R[i].forEach((c, j) => {
      if (!c) return;
      const sign = c < 0 ? '-' : '+';
      const magnitude = Math.abs(c) === 1 ? '' : String(Math.abs(c));
      term += `${sign}${magnitude}${'xyz'[j]}`;
    });
    const tau = fractionText(t[i]);
    if (tau) term += `+${tau}`;
    return term.replace(/^\+/, '') || '0';
  }).join(', ');
}

/**
 * A Hermann–Mauguin symbol as CIF spells it, one blank between the lattice letter and each
 * position, subscripts inline: 'P2_1/c' → 'P 21/c', 'Fd-3m' → 'F d -3 m',
 * 'I4_1/amd' → 'I 41/a m d', 'P6_3mc' → 'P 63 m c'.
 */
export function cifSpaceGroupName(symbol) {
  const lattice = symbol[0];
  const parts = symbol.slice(1).match(/-?\d(?:_\d)?(?:\/[a-z])?|[a-z]/g) || [];
  return [lattice, ...parts.map((p) => p.replace('_', ''))].join(' ');
}

// A CIF value: quoted when it is empty or holds a blank or a quote-sensitive character.
const value = (text) => (/^[^\s'"_#$;[\]][^\s]*$/.test(text) ? text : `'${text.replace(/'/g, '')}'`);

const ascii = (text) => String(text)
  .replace(/Å/g, 'A').replace(/×/g, 'x').replace(/≥/g, '>=').replace(/[–—]/g, '-')
  .replace(/[^\x20-\x7e]/g, '?');

/** A small rational as text ('2', '1/2', '3/2'), decimals when it is not one. */
export function rationalText(c) {
  for (let d = 1; d <= 48; d++) {
    const k = Math.round(c * d);
    if (Math.abs(c * d - k) < 1e-6 * d) {
      const g = gcd(Math.abs(k), d) || 1;
      return d / g === 1 ? String(k / g) : `${k / g}/${d / g}`;
    }
  }
  return fixed(c, 4);
}

// The output basis as text, column j of Q being the new vector j: "a' = a+b, b' = -a+b, …".
const basisText = (Q) => ['a', 'b', 'c'].map((name, j) => {
  let out = '';
  [0, 1, 2].forEach((i) => {
    const c = Q[i][j];
    if (Math.abs(c) < 1e-9) return;
    const magnitude = Math.abs(Math.abs(c) - 1) < 1e-9 ? '' : rationalText(Math.abs(c));
    out += `${c < 0 ? '-' : '+'}${magnitude}${'abc'[i]}`;
  });
  return `${name}' = ${out.replace(/^\+/, '')}`;
}).join(', ');

const formulaText = (counts) => Object.entries(counts)
  .sort(([p], [q]) => (p < q ? -1 : 1))
  .map(([el, n]) => {
    const rounded = Math.abs(n - Math.round(n)) < 1e-3 ? String(Math.round(n)) : fixed(n, 3).replace(/0+$/, '');
    return rounded === '1' ? el : `${el}${rounded}`;
  })
  .join(' ');

/**
 * The CIF text of a symmetryAveragedStructure() model.
 * @param {object} model
 * @param {{ dataName?:string, date?:string, tolRange?:[number, number] }} [options]
 *   dataName — the data block name (sanitized here); date — ISO date for the audit record;
 *   tolRange — the ladder range over which the group holds (Å), for the header comment.
 */
export function writeCif(model, { dataName = 'rmc_average', date = null, tolRange = null } = {}) {
  const { spaceGroup, cell, operations, sites, formula, provenance: p } = model;
  const block = String(dataName).replace(/[^A-Za-z0-9_.-]+/g, '_') || 'rmc_average';
  const source = p.source ? String(p.source).split('/').pop() : 'the .rmc6f model';
  const lines = [];
  const comment = (text = '') => lines.push(`# ${ascii(text)}`.trimEnd());

  comment('Symmetry-averaged RMC average structure, written by rmc-toolkits.');
  comment();
  const box = p.supercell ? `, supercell ${p.supercell.join('x')}` : '';
  const atoms = Number.isFinite(p.atoms) ? `${p.atoms} atoms, ` : '';
  comment(`Source: ${source} (${atoms}${p.sites} reference sites${box}).`);
  const range = tolRange ? `holds from ${fixed(tolRange[0], 3)} to ${fixed(tolRange[1], 3)} A on the tolerance ladder; ` : '';
  const residual = Number.isFinite(p.maxResidual) ? `, worst residual ${fixed(p.maxResidual, 3)} A` : '';
  comment(`Space group: ${spaceGroup.label}${spaceGroup.number ? ` (No. ${spaceGroup.number})` : ''}; ${range}`
    + `orbits taken at ${fixed(p.tol, 3)} A.`);
  comment(`  ${p.nSpace} operation${p.nSpace === 1 ? '' : 's'} in the .rmc6f unit cell${residual}.`);
  comment(spaceGroup.standard
    ? `Cell: the standard cell the group is named in: ${basisText(p.Q)} (a, b, c: the .rmc6f unit cell).`
    : 'Cell: the .rmc6f unit cell (no standard setting was found for this group, so no symbol or number is given).');
  comment(`Origin: x(CIF) = inv(Q) (x(.rmc6f) + d), d = (${p.originShift.map((v) => fixed(v, 5)).join(', ')})`);
  if (p.itaOrigin) {
    comment(`  in .rmc6f cell fractions (${fixed(p.originShiftA, 4)} A): the ITA standard origin`
      + `${p.originChoice ? ` (origin choice ${p.originChoice})` : ''}, of the equivalent ones the nearest the .rmc6f origin.`);
  } else {
    comment(`  in .rmc6f cell fractions (${fixed(p.originShiftA, 4)} A)${p.niceOrigin ? ', the smallest shift that puts every'
      : '; no nearby origin puts the'} translation${p.niceOrigin ? ' on the 1/48 grid.' : 's on the 1/48 grid, so they are decimals.'}`);
  }
  comment('Positions: each site\'s arithmetic mean over its copies in the box, averaged over its orbit');
  comment('  through the exact operations (special positions are exact; free coordinates are measured).');
  comment(`  Largest move of a site mean onto its symmetrized position: ${fixed(p.maxShiftA, 4)} A (rms ${fixed(p.rmsShiftA, 4)} A).`);
  if (p.unmatched > 0) comment(`  ${p.unmatched} site(s) no operation reached did not enter the average.`);
  comment('ADPs: U_ij = mean-square displacement of every atom of the orbit about its symmetrized');
  comment('  position (static + thermal, in A^2); it includes any distortion the picked group averages away.');
  comment('Occupancy: atoms of each element / (sites in the orbit x unit cells in the box).');
  comment(p.itaOrigin
    ? 'The symmetry operations below are those of the ITA standard setting.'
    : 'The symmetry operations below are authoritative; the origin is not the ITA standard one.');
  lines.push('');

  lines.push(`data_${block}`);
  lines.push(`_audit_creation_method                 ${value('rmc-toolkits symmetry-averaged RMC structure')}`);
  if (date) lines.push(`_audit_creation_date                   ${date}`);
  lines.push(`_chemical_formula_sum                  '${formulaText(formula.counts)}'`);
  lines.push(`_cell_formula_units_Z                  ${formula.Z}`);
  lines.push('');
  lines.push(`_cell_length_a                         ${fixed(cell.a, DIGITS.cell)}`);
  lines.push(`_cell_length_b                         ${fixed(cell.b, DIGITS.cell)}`);
  lines.push(`_cell_length_c                         ${fixed(cell.c, DIGITS.cell)}`);
  lines.push(`_cell_angle_alpha                      ${fixed(cell.alpha, DIGITS.angle)}`);
  lines.push(`_cell_angle_beta                       ${fixed(cell.beta, DIGITS.angle)}`);
  lines.push(`_cell_angle_gamma                      ${fixed(cell.gamma, DIGITS.angle)}`);
  lines.push(`_cell_volume                           ${fixed(cell.volume, 3)}`);
  lines.push('');
  if (spaceGroup.system) lines.push(`_space_group_crystal_system             ${spaceGroup.system}`);
  if (spaceGroup.symbol) {
    const setting = spaceGroup.centring === 'R' ? ' :H' : p.originChoice ? ` :${p.originChoice}` : '';
    const name = cifSpaceGroupName(spaceGroup.symbol) + setting;
    lines.push(`_space_group_name_H-M_alt              '${name}'`);
    lines.push(`_symmetry_space_group_name_H-M         '${name}'`);
  }
  if (spaceGroup.number) {
    lines.push(`_space_group_IT_number                 ${spaceGroup.number}`);
    lines.push(`_symmetry_Int_Tables_number            ${spaceGroup.number}`);
  }
  lines.push('');
  lines.push('loop_');
  lines.push('_space_group_symop_id');
  lines.push('_space_group_symop_operation_xyz');
  operations.forEach(({ R, t }, i) => lines.push(`${String(i + 1).padStart(4)}  '${symopText(R, t)}'`));
  lines.push('');

  // One row per element of each site: a mixed site is several rows at one position.
  const labelCount = {};
  const rows = [];
  for (const site of sites) {
    for (const { element, occupancy } of site.elements) {
      labelCount[element] = (labelCount[element] || 0) + 1;
      rows.push({ site, element, occupancy, label: `${element}${labelCount[element]}` });
    }
  }
  const anisotropic = (site) => site.Ueq > 1e-8;
  lines.push('loop_');
  for (const tag of ['label', 'type_symbol', 'symmetry_multiplicity', 'Wyckoff_symbol', 'fract_x', 'fract_y', 'fract_z',
    'occupancy', 'U_iso_or_equiv', 'adp_type']) lines.push(`_atom_site_${tag}`);
  for (const { site, element, occupancy, label } of rows) {
    lines.push([
      label.padEnd(6),
      element.padEnd(3),
      String(site.multiplicity).padStart(3),
      site.wyckoff ?? '?',
      ...site.x.map((v) => fixed(v, DIGITS.xyz).padStart(9)),
      fixed(occupancy, DIGITS.occupancy),
      fixed(site.Ueq, DIGITS.U),
      anisotropic(site) ? 'Uani' : 'Uiso',
    ].join(' '));
  }
  const aniso = rows.filter(({ site }) => anisotropic(site));
  if (aniso.length) {
    lines.push('');
    lines.push('loop_');
    for (const tag of ['label', 'U_11', 'U_22', 'U_33', 'U_12', 'U_13', 'U_23']) lines.push(`_atom_site_aniso_${tag}`);
    for (const { site, label } of aniso) {
      const U = site.U;
      lines.push([label.padEnd(6), ...[U[0][0], U[1][1], U[2][2], U[0][1], U[0][2], U[1][2]]
        .map((v) => fixed(v, DIGITS.U).padStart(9))].join(' '));
    }
  }
  lines.push('');
  return lines.join('\n');
}
