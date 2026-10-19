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
  version on the user's decision. (Done at protocol 2026-10-19.11: the five matches are in the
  test set, sealed since the holdout gate, and the dev cases are no longer `devOnly`. At protocol
  2026-10-19.19 the marginal situation's knob changed for the engine's calibrated null, and the dev
  pair was re-set from it; its two matches are to be re-set from the same knob, below.)

## What not to read

- `validation/judgmentCases.js`: the dev cases' parameters.
- `validation/baseline/arm-a-dev.json`: the dev records, judgment cases included.
- The judgment pins in `validation/__tests__/arms.test.js` (its last five `BAR` rows and the test
  "the judgment subset: where A+'s fixed joining rules fail"), and the tests that read the dev
  marginal pair's values: `judgment.test.js`'s marginal-pair and edge-of-significance tests, and
  the weak-order tests in `src/workers/__tests__/modulations.test.js` and
  `src/__tests__/copilotTools.test.js`.
- The PR that added the dev cases, and any arm record or transcript on them.
- `validation/findings/decisions/2026-10-19-marginal-reset-dev.json` and its PR: the dev marginal
  pair's amplitude, seeds and results.

Read `validation/cases.js` for the generators: the kinds, their options and the truths they write.

## How to write a match

- **One match per dev case: five.** Each has `judgment: true`, the same `situation`, and the same
  kind and option. Change the elements, directions, magnitudes or structures.
- **The same physics, a comparable difficulty.** Set the knobs below from the physics and from the
  engines' own rules, never from a dev value.
- **Seeds by a rule fixed before you look.** The rule depends only on your case's own draws or an
  engine statistic of your case, never on an arm's score. The precedent is protocol 2026-10-19.2's
  seed rule (docs/algorithms/copilot-validation.md; its numbers are sealed).
- **Check the situation's property on your own configuration.** Use the engines and the partial
  g(r) directly. Do not score arm A, arm A+ or a model on your cases before they are frozen.
- **Coincidences.** The lead compares your elements, directions and magnitudes with the dev cases'
  without showing them to you, and asks you to change any that coincide.
- **Where they go.** The sealed test set: `TEST_CASES` in the held-out module that only
  `validation/holdoutLoader.js` imports (docs/algorithms/copilot-validation.md, "The sealed test
  set"), and the id in `TEST_CASE_IDS` (`validation/cases.js`). `TEST_SET_SHA256` in
  `validation/__tests__/split.test.js` changes, with a dated protocol note in
  docs/algorithms/copilot-validation.md that states no value. The split test's judgment rules then
  pair each dev judgment case with its match.

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
- **The truth.** `modulation` on B at q. No axis: no single-site reading resolves it at this amplitude.
  (Since protocol 2026-10-19.12 `analyze_modulations` gives an ordered point's polarisation, the
  direction of the mode; the truth does not score it.)
- **Two cases at one amplitude.** One seed realised on each side of the bar.
- **What the score sees.** The verdict has no confidence field, so a claim in the VERDICT counts
  however hedged the prose is. A claim needs the element and q (within the scoring's wavevector
  tolerance). The app's answer check on a superstructure claim needs an ordered point (3σ or more),
  so a claim below the bar gets a finding beside its score. That tension is the point.
- **Re-set at protocol 2026-10-19.19.** Until then this section took the modulation engine's null
  before finding F1 (protocol 2026-10-19.16), which held a site to one component. The engine now
  tests against independent cells, and one isotropic site has three components: an amplitude set
  from the earlier recipe reads far above the bar, so the situation is no longer marginal. A match
  made from that recipe is re-set from the one below.
- **Knobs.**
  - **The expected significance: 3σ.** The median of the realised significance at q over
    configurations drawn with the case's parameters (fresh seeds): half of them clear the bar, which
    is what an expected significance means here. The realised significance is B's σ at q in the
    default `analyze_modulations` call (every element), the engine's unrounded value. Compute the
    expectation from the engine's null, as follows (the engine: `workers/modulations.js`, with
    `rmc_toolkits/modulations.py` its reference).
    1. **The power at q.** One B site of n copies (n is the supercell's cell count), thermal rms σ
       along each axis (the kind fixes it: `THERMAL.B` in `validation/cases.js`), displaced ±δ
       along one axis in the antiphase pattern of a zone-boundary q. The engine transforms the
       site's displacements u(c) from its mean position over the cell grid,
       U(q) = Σ_c u(c) e^(−2πi q·c), a real sum at a zone-boundary point. Each component of U is
       normal with variance nσ² (n independent thermal terms). The two components across the
       displacement have mean 0, and the one along it has mean nδ. The power is
       P = |U|² / (n tr C), where tr C = 3σ² + δ² is the site's total displacement variance, the
       static ±δ included. So
       `P = σ²/(3σ² + δ²) · X`, `X ~ χ'²₃(Z)`, `Z = n δ²/σ²`.
       X is a noncentral χ² with three degrees of freedom (two central components, one noncentral)
       and noncentrality Z. The mean of P is (Z + 3)/(3 + δ²/σ²).
    2. **The null.** The independent-cell null places the configuration's own cells at random, which
       keeps every correlation inside a cell and destroys every correlation between cells. At a
       zone-boundary point the element's power under it is `P₀ = Σ_k λ_k χ²₁`, a sum of independent
       χ²₁ terms. Its weights λ_k are the eigenvalues of the element's cell covariance
       N/(N − 1) Σ_c y(c) y(c)ᵀ, normalised as the power (`independentCellNull`): N is the cell
       count, y(c) holds the displacements of the element's sites in cell c, and the weights sum to
       N/(N − 1), the mean power. For the B site of step 1,
       `λ = n/(n − 1) · (σ², σ², σ² + δ²) / (3σ² + δ²)`: three components, one per axis the site
       moves along, where the earlier recipe took one.
    3. **The family-wise significance of a power P.** The p-value is p = Pr(P₀ ≥ P)
       (`nullTail(P, λ, true)`: a Lugannani–Rice saddlepoint, 0 to 8 % above the exact tail). Then
       σ = Φ⁻¹((1 − p)^F) (`zoneBoundarySigma(P, λ, F)`). F, the Šidák family, is the number of
       zone-boundary points the supercell holds (7 when every axis is even) times the number of
       elements in the configuration with displacements (three in an ABX₃ case), whichever element
       a call names.
    4. **The expected significance** is the σ of step 3 at the median power,
       `P_med = σ²/(3σ² + δ²) · median χ'²₃(Z)`, since σ increases with P. The median is
       `scipy.stats.ncx2.median(3, Z)` in Python. In JavaScript it is the root x of
       Σ_j e^(−Z/2) (Z/2)^j / j! · P(3/2 + j, x/2) = ½, where P is the lower regularized gamma
       (`regularizedGamma(a, x).lower`, `workers/orientation.js`). Solve σ(P_med) = 3 for δ.
       - At a fixed F the result depends on Z, and only weakly on δ/σ: for F = 21, σ(P_med) = 3 at
         Z ≈ 20. The amplitude therefore scales as σ/√n, so the supercell sets it.
       - Do not take the σ at the mean power: it reads about 0.14σ higher than at the median.
       - The formula leaves out how tr C and λ vary from one configuration to the next; the engine
         check below includes it.
    5. **The ratio bar.** An ordered point also needs a power of at least 3 × the median power over
       q. Near an expected 3σ that bar is passed by a wide margin, because an isotropic site's median
       power is about 0.9.
  - **The scatter.** At an expected 3σ the realised significance has a standard deviation of about
    1.4σ, with a longer tail below than above. About 30 % of configurations land in [3, 4] and
    about 14 % in [2.5, 3).
  - **Check it with the engine** before taking seeds:
    - Build K configurations of your case on seeds that no case uses and that lie outside your seed
      sequence (`buildCase(item, { partials: false })` skips the g(r)).
    - Run the engine as the default call does, every element:
      `displacementModulations(siteDisplacementsFromRmc6f(rmc6f))` (`workers/modulations.js`,
      `workers/pcaKde.js`). Read B's σ at q.
    - Their median should match your expected significance to within its sampling error, about
      1.25 × 1.4σ/√K.
    - About half of them should reach 3σ, within the binomial interval.
    - B's `nullDof` should be close to 3, one isotropic site's three components.
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
