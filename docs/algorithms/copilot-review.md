# AI Copilot: tool and skill review (2026-10-14)

A review of the Copilot's analysis tools, answer checks and skills, made with code reading and the
synthetic configurations of the [validation protocol](copilot-validation.md): supercells with
planted local structure, where the right answer is known. No model was run.

## Tools

"Known answer" means a test that plants a structure and checks the tool recovers it, not a
regression pin on the demo run.

| Tool | What it computes | Known-answer test | Gaps |
| --- | --- | --- | --- |
| `list_sites` (PCA) | Per site: U_iso, principal rms, rms1/rms3, Mardia kurtosis | `test_pca_kde.py`: known covariance and anisotropy, kurtosis flags a non-Gaussian cloud; validation: two-site splits recovered (direction ≤ 10°, δ within 2 %) | Mardia kurtosis only; the per-axis kurtosis that shows a split along PC1 is not in the tool's rows |
| `analyze_displacement_directions` | Sphere histogram per site; calibrated peak, map, antipodal and Bingham tests | `test_orientation.py`: isotropic cloud flat, lobed cloud peaks on its axis | **L3:** "significant" includes the map test, and any anisotropic Gaussian ellipsoid passes it (thermal-only O site: peak 2.3σ, map 7.2σ), so the direction answer check lets "preferred direction" through on plain thermal anisotropy |
| `compute_bond_angles` | A–B–C angles, sin θ-corrected histogram, CN, ideal-polyhedron split; refuses a window that cuts a shell | `test_triplets.py`: ideal geometries, Se–Nb–Se octahedra | No `facts`, so angle and CN claims are not checked |
| `find_gr_peaks` + the Workbench's window questions | First shell from the partial g(r): window, peaks, CN, cut sensitivity | `clarify.test.js`: split shells, tolerances | **L1:** a Jahn–Teller shell whose two bond lengths are well resolved (4 × 1.93 + 2 × 2.25 Å) is cut at the deep trough: `first_shell` is 1.70–2.10 Å with CN 4, and the long bonds look like "the next shell" (only a "not separated" warning). Following the coordination-shells skill gives a four-coordinate Cu |
| `analyze_clusters` | Finite same-element clusters vs a shuffle null: edge ladder, breathing, E/T2, gyration axes | `test_clusters.py`: trigonal, static, dimers, wrap at fraction 0; validation: rigid dimers at x = 0 (breathing −87 %, z = −37), breathing triangles (+49 %, z = 22) | **L2:** the null shuffles each site's displacements whole, so collective motion lowers the shape statistics: breathing clusters come out "28 % LESS distorted", and the gyration-axis test reads breathing tetrahedra as a shape preference (a spurious "distortion" in the A+ sweep). For split sites (JT) the null already holds the distortion's amplitude. No `facts` |
| `analyze_modulations` | Per-element displacement power at the zone-boundary points vs the configuration's own background | `test_modulations.py`: planted doubling, under a stronger long wave, independent atoms; validation: (0,0,½) and (½,½,½) found, decoy long waves not called order | Fixed below. Still blind to incommensurate or short-range order between cells |
| `detect_symmetry` + context ladder | Space group at a tolerance; point-group facts; noise floor | the symmetry suites (230 groups, spglib-pinned); validation: F-43m, P-43m, P2₁3, P4/mmm, P4mm and R3m found | The card's 0.2 Å default misses a 0.06–0.07 Å polar shift (P4mm, R3m read as Pm-3m); the context's ladder and noise floor show it |
| Answer checks (`answerAudit.js`) | Regex claims vs tool facts: convergence, chirality, centrosymmetry, polar, inversion lost, direction, superstructure | `answerAudit.test.js` | They cover positive symmetry, convergence and direction claims only. A missed feature, or ordered order read as disorder, is caught by no check; the validation's scoring measures those |

### Fixed in this review (each with a regression test)

1. **Superstructure false positives.** A zone-boundary point counted as order when it ranked in
   the top 1 % of the modes, with no correction for the 7 points × elements tested. Real-valued
   zone-boundary sums also have half the degrees of freedom, and so a heavier tail. The rule found
   "order" in 6 of 12 configurations without any, the thermal null included. Now a point needs
   ratio ≥ 3 and ≥ 3σ against independent displacements: χ² with at least m degrees of freedom
   for m sites, which holds for any site anisotropy, Šidák-corrected. Planted order stays above 8σ.
   (`zoneBoundarySigma`, in the modulation engines since 2026-10-19: workers/modulations.js,
   modulations.py, and their tests)
2. **The cluster check denied a distortion its own axis test found.** For JT octahedra it printed
   "LESS distorted … no correlated shape distortion on top" next to "more clusters elongated along
   <100>". The sentence now names the favoured axis.
3. **The polar answer check** that the local-symmetry skill promised did not exist. Added.
4. **Superstructure answer check.** `analyze_modulations` now returns `facts`, and a claimed
   superstructure is flagged when every modulation analysis found none.
5. **A JSON key is not a claim.** `"chiral": false` in a structured block was read as calling the
   group chiral.

### Larger findings (not changed; tool and skill changes need the user's approval)

- **L1, split shells cut at a deep trough** (above). Proposal: treat a peak within ~0.5 Å of the
  first as the same shell when the summed CN reaches an ideal polyhedral count (4, 6, 8, 12), and
  let the Workbench ask. Touches `firstShell`, `clarify.js` and the coordination-shells skill.
- **L2, a cluster null that removes collective motion.** Proposal: shuffle after taking out each
  cluster's mean radial scale (breathing), report breathing on its own, and say "split sites:
  the null already holds the amplitude" when the sites' clouds are bimodal.
- **L3, "preferred direction" means a peak, not any anisotropy.** Proposal: `facts.direction.significant`
  from the peak test only, and the map test as "anisotropic". It changes a documented rule.
- **L4, missing tools.** The coherent shift of a symmetry lowering has no magnitude tool: the orbit
  positions are in `detect_symmetry`, but the shift must be worked out by hand. Nothing tests
  correlations between cells beyond the zone-boundary points.
- **L5, unchecked claims.** Bond angles and clusters return no `facts`. A cluster rule would
  misfire today because of L2.

### Implemented in validation batch 1 (approved by the user on 2026-10-14; 2026-10-17)

Every tool proposal above (the skill proposals below go in a PR of their own, which waits for the
user's review), each recorded as a finding in
`web_app/frontend/validation/findings/batch-1.jsonl` with its commit and tests:

| Item | Change | Commit |
| --- | --- | --- |
| L1 | `firstShell` joins a later peak within 0.5 Å when the summed CN is 4, 6, 8 or 12 (either direction); the first peak alone is `first_peak_only`, offered and asked about (clarify.js) | 0a7fdf0 |
| L2 | `clusters.py` / `clusters.js` `shapeNull`: each cluster keeps its breathing, only the breathing-free displacements are shuffled, their radial component scaled by 1/√(1 − β_m); shape rows read against it, size rows against independent atoms; split member sites named | 6b308f3 |
| L3 | `facts.direction.significant` from the peak test; `anisotropic` from the map test | b867260 |
| L4 | `measure_symmetry_lowering` (`symmetryShift.js`); `analyze_modulations` tests general wavevectors | 56cd672 |
| L5 | coordination facts from `compute_bond_angles` and `find_gr_peaks`; cluster and coordination answer checks | 4481a86 |

Two consequences the proposals did not spell out: the per-cluster breathing fit takes out a share
β_m of each member's own radial noise, which the shape null puts back (without it independent atoms
read ~6 % too regular); and against the shape null a trigonal push of one vertex keeps its pattern
(trigonality, the <111> axis family) but no longer its excess edge spread, which the vertex sites'
own two positions hold.

## Skills

All four are drafts awaiting the user's interview. Changes to them need the user's approval; these
are proposals; batch 1 implements them in a separate PR, with the missing local-order skill as a fifth draft.

| Skill | Every decision tied to a tool output? | Criteria explicit and physical? | Contradictions and gaps |
| --- | --- | --- | --- |
| rmc-fit-assessment | Yes: `assess_convergence`, `compare_fit`, `find_gr_peaks`, the context's weights | Convergence by the watchdog; a good Rw, sampling and overfitting are left as judgment calls | Not exercised by the synthetic cases (no fit data). Box-scale waves are named, but nothing says when they void a modulation claim |
| coordination-shells | Yes: `find_gr_peaks` first shell, CN, refusal, sin θ histogram | The 60 % modal rule and the CN comparison are explicit | **Contradiction (L1):** "use `first_shell.window_A`" and "a shell with two peaks is one split shell" disagree when the trough between the peaks is deep: the JT octahedron reads as CN 4 |
| displacement-anisotropy | Yes: `list_sites`, directions, `analyze_clusters` | The σ thresholds are; "name the family only when the angle is small" has no number (proposal: 15°, the scoring's) | **Gap:** a split site (negative kurtosis) is read as disorder. Ordered antiphase displacements give the same single-site statistics, so `analyze_modulations` must come first. The zone-boundary case is built on exactly this |
| local-symmetry | Yes: ladder, `noise_floor_A`, point-group facts | Explicit; noise floor = 3√2 × rms/√copies | The polar check it promised was missing (fixed). It leaves "which tolerance is *the* group" open. The symmetry-lowering cases show the 0.2 Å card misleads; proposal: report the lowest rung above the noise floor |
| *(missing)* local order | — | — | No skill says how to conclude "no order": which analyses must all come back empty (directions per site, clusters, modulations, the ladder above the noise floor). Proposal: a short skill, with the user's approval |
