# The judgment subset: a spec for its held-out cases

Plan item 1.6, part ii. This file is what the author of the held-out judgment cases works from, and
all they read about the dev judgment cases: the four situations, the kind and option each dev case
uses, what makes it a judgment case, and the knobs that set its difficulty. It holds no parameter of
a dev case, no dev result and no arm score, and it must stay that way:
`validation/__tests__/judgment.test.js` scans it for the dev cases' ids and values.

## Why the subset exists

- **The ceiling.** The engines alone (arm A+) score close to 1 on the paired cases. Most of a case
  score is made of facts a fixed rule gets right: the symmetry fields and the spurious-claim test.
  On those cases a model can at best tie A+.
- **The subset.** Cases where A+'s fixed joining rules (`validation/armA.js`, `runArmSweep`) are
  ambiguous or wrong, so a model that reads the tools' results with judgment has room to add. One
  situation, the user's choice, is a control on the score (below).
- **Dev now, held out later.** The dev cases exist (`judgment: true`). Reported numbers need
  held-out matches, written by you without the dev results, joining the test set as a new pinned
  version on the user's decision. (Done at protocol 2026-10-19.11: the five matches are in
  `TEST_CASES`, and the dev cases are no longer `devOnly`.)

## What not to read

- `validation/judgmentCases.js`: the dev cases' parameters.
- `validation/baseline/arm-a-dev.json`: the dev records, judgment cases included.
- The judgment pins in `validation/__tests__/arms.test.js` (its last five `BAR` rows and the test
  "the judgment subset: where A+'s fixed joining rules fail").
- The PR that added the dev cases, and any arm record or transcript on them.

Read `validation/cases.js` for the generators: the kinds, their options and the truths they write.

## How to write a match

- **One match per dev case: five.** Each has `judgment: true`, the same `situation`, and the same
  kind and option. Change the elements, directions, magnitudes or structures.
- **The same physics, a comparable difficulty.** Set the knobs below from the physics and from the
  engines' own rules, never from a dev value.
- **Seeds by a rule fixed before you look.** The rule depends only on your case's own draws or an
  engine statistic of your case, never on an arm's score. The precedent is protocol 2026-10-19.2's
  seed rule (docs/algorithms/copilot-validation.md).
- **Check the situation's property on your own configuration.** Use the engines and the partial
  g(r) directly. Do not score arm A, arm A+ or a model on your cases before they are frozen.
- **Coincidences.** The lead compares your elements, directions and magnitudes with the dev cases'
  without showing them to you, and asks you to change any that coincide.
- **Where they go.** `TEST_CASES` in `validation/cases.js`. `TEST_SET_SHA256` in
  `validation/__tests__/split.test.js` changes, with a dated protocol note in
  docs/algorithms/copilot-validation.md. The split test's judgment rules then pair each dev judgment
  case with its match.

## The four situations

### 1. Two features at once (`two_features`)

- **The rule at issue.** A+ lets one finding explain another reading. A split coordination shell
  around a centre that is itself a split site is taken as the off-centring's own shell
  (`runArmSweep`: `if (found.has(shell.centre)) continue`). With two features on one centre, that
  rule keeps one of them.
- **The kind.** `jahn_teller` with the option `offCentre: { axis, delta }`:
  - isolated Jahn–Teller octahedra, 4 short + 2 long bonds, the long axis random over x, y and z;
  - the centre also sits ±delta off along one fixed axis, the sign random per atom.
- **The truth.** `jahn_teller` and `off_centre` on the centre, both planted and scored; the kind's
  tolerated partial descriptions (the ligand's `off_centre` and `cluster_distortion`).
- **What a good answer does.** It sees that an off-centred centre splits its shell 1 + 4 + 1 about
  the mean distance, and that this cannot put four ligands below the mean and two above it. Both
  distortions are there. The evidence:
  - the shell's peaks against the average-structure distance;
  - the coordination number of each peak;
  - the ligands' +u/−u asymmetry along their bonds.
- **Knobs.**
  - The centre's shift against its thermal width: large enough that the centre is a split site by
    the one split rule (`isSplitSite`, `workers/clusters.js`; two equal wells reach the kurtosis bar
    at a shift equal to the thermal rms).
  - The Jahn–Teller split (long − short) against the shift: clearly larger, so that the shell's two
    peaks are the elongation's 4 + 2, not the off-centring's satellites.
  - The cell: large enough that the octahedra stay isolated.

### 2. A signal at the edge of significance (`marginal`)

- **The rule at issue.** A fixed 3σ bar. At 2.5–4σ, whether a planted signal clears it is a coin
  flip. An honest answer gives the significance and says how confident it is.
- **The kind.** `zone_boundary` with `weak: true` and `decoy: null`:
  - an ordered antiphase displacement of the B cation at a zone-boundary q;
  - its amplitude below half the B thermal rms (the generator refuses more), so no single-site
    reading shows it.
- **The truth.** `modulation` on B at q. No axis: no tool resolves it at this amplitude.
- **Two cases at one amplitude.** One seed realised on each side of the bar.
- **What the score sees.** The verdict has no confidence field, so a claim in the VERDICT counts
  however hedged the prose is. A claim needs the element and q (within the scoring's wavevector
  tolerance). The app's answer check on a superstructure claim needs an ordered point (3σ or more),
  so a claim below the bar gets a finding beside its score. That tension is the point.
- **Knobs.**
  - **The expected significance: 3σ.** Set it from the zone-boundary power after the family-wise
    correction over every zone-boundary point and element tested (`zoneBoundarySigma`,
    `workers/modulations.js`). Take one B site of n copies, thermal rms sigma per axis, amplitude
    delta along one axis. Along that axis the power is a noncentral χ² of one degree of freedom
    with non-centrality n·delta²/sigma². The two other axes are central. So
    E[P] ≈ (n·delta²/sigma² + 3) / (3 + delta²/sigma²). The realised significance scatters by about
    one σ around its expectation.
  - **The seeds.** A pair, by a rule that does not depend on which side an arm lands: in a fixed
    sequence, the first seed whose realised significance at q lies in [3, 4], and the first in
    [2.5, 3).
  - **The q.** Any zone-boundary point the supercell holds (its axes even).

### 3. Two tools read differently (`two_readings`)

- **The rule at issue.** A+ calls a site split only from its ellipsoid, by the app's one split rule
  (`isSplitSite`, `workers/clusters.js`): PC1 resolved, and PC1's excess kurtosis at or below the
  split bar and clearly below zero. Wells arranged cubically average to a round ellipsoid. No axis
  is resolved, nothing is called split, and the split coordination shell the wells cause is left
  unexplained.
- **The kind.** `offcentre` with `wells: 8` (`axis` a <111> direction): each B atom sits delta off
  its site along one of the eight <111> directions at random (eight-site disorder).
- **The truth.** `off_centre` on B, its direction <111> in any cubic image
  (`direction_mode: 'cubic'`), its magnitude delta. Nothing is tolerated, as in the kind.
- **The readings that differ.**
  - The ellipsoid and the orientation tensor: isotropic.
  - The Mardia kurtosis (`list_sites`): strongly negative.
  - The direction map: anisotropic at its default tiling (the map test). Its eight <111> lobes are
    resolved by the peak test only on a finer tiling (a higher `frequency`).
  - The B–X shell: split 3 + 3 by the off-centring itself. It is not a Jahn–Teller distortion: a
    `jahn_teller` claim on B is spurious.
- **Knobs.**
  - The shift against the B thermal rms: a lobe's angular width is about their ratio.
  - The copies per site (the supercell): the default tiling has about 12 points per cell, so more
    copies resolve the lobes at a finer default.
  - The 3 + 3 split, 2·delta/√3, against the shell's width.

### 4. A choice only the user can make (`user_choice`)

- **The choice.** A window.
  - The Workbench joins a split shell whose peaks close at an ideal coordination number. It asks
    the user when an analysis' window cuts that shell (`src/llm/copilot/clarify.js`).
  - The simulated user skips, so the default (the whole shell) stands. A good answer states that
    default and what the other window would give. The truth follows the default.
- **The control.** Since protocol 2026-10-19.8, A+ takes its cluster windows from the same shell
  reading as the Workbench (`src/llm/context/shells.js`, `validation/armA.js`), so its fixed rule
  and the default coincide. On the score this situation is a control; the judgment shows in the
  answer only. The dev case was designed against the earlier rule (a cut at the first trough of the
  M–M g(r)). Match it as it is: it pins the question and its default, and it would show a reader
  that leaves the one shell reading again.
- **The kind.** `elongated_clusters`, the subset's one kind of its own (`validation/cases.js`):
  - finite M6 octahedra around the body centre, an anion at the origin;
  - each octahedron elongated along x, y or z at random, its two apical atoms pushed out by
    `elongation` Å;
  - locally 4 short (equatorial, `edge`) and 8 long edges; the average octahedron is regular.
- **The truth.** `cluster_distortion` on M. Tolerated: `off_centre` and `jahn_teller` on M (each
  vertex sits on one of two positions along its own axis, and each M's own M–M shell is split).
  Symmetry: Pm-3m.
- **Knobs: constraints the configuration must meet.**
  - The M–M shell is split across a deep trough: its first minimum is below
    `SHELL_TROUGH_FRACTION` of the first peak (`src/llm/context/shells.js`). This sets the
    elongation against the edge's thermal width.
  - The joined shell closes at 4 M per M (an octahedron), within `SPLIT_SHELL_REACH_A` of the first
    peak.
  - The average octahedron's edge lies beyond the first trough, so the first peak alone holds no
    edge of the average structure.
  - The cell is large enough that the long edges decay into an empty gap before the nearest M–M
    contact between clusters.

## How the report reads the subset

- `run.js report` gives each split with judgment cases a table: each column's mean score per
  situation, over the subset and over the rest.
- The comparison beside the case score (plan item 1.6, part i) lists their components and
  disagreements as for any case.
- Until the held-out matches exist, every judgment result is a dev result.
