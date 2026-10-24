# Order signatures: where each kind of order shows

Measured with the Copilot's tools on the validation's dev cases (synthetic supercells whose order
is known; 2026-10-17, after validation batch 1). Use them to recognise a pattern, never as another
run's numbers.

| Order | Shows in | Hides from |
|---|---|---|
| Two-site disorder (each B atom ±δ along an axis, at random) | `list_sites`: negative kurtosis, prolate; directions: a peak along the axis; `find_gr_peaks`: the B–X shell split by the off-centring | `analyze_modulations` (no correlation between cells); the ladder (the average site stays) |
| Ordered antiphase displacements (a superstructure) | `analyze_modulations`: the zone-boundary point, ratio ≥ 3 and ≥ 3σ | `list_sites` and directions: the same split cloud as two-site disorder; the ladder |
| Order at a general wavevector (incommensurate, short-range between cells) | `analyze_modulations`: a strongest mode that is neither zone boundary nor long wave, ratio ≥ 3 and ≥ 3σ | the average structure; single-site statistics |
| Breathing clusters (size fluctuates collectively, shape regular) | `analyze_clusters`: breathing above independent atoms | every page; the cluster's shape at its own size |
| Rigid clusters (atoms of a cluster move together) | `analyze_clusters`: shape spread or bond-length spread below independent atoms | every page |
| Distorted polyhedra, disordered orientation (Jahn–Teller) | `find_gr_peaks`: a split shell (first peak alone a lower CN); split X sites; `analyze_clusters` on the X polyhedron: an axis-family pattern, split sites named | the average structure (a regular average polyhedron) |
| A coherent lowering of the average structure | the ladder above the noise floor; `measure_symmetry_lowering`: the shift | the Detected SG card at 0.2 Å when the shift is smaller; displacement directions (the cloud is centred on the shifted mean) |

## The dev cases

- **Two-site disorder** (Nb ±0.20 Å along [001] in KNbO₃): Nb kurtosis −0.43, rms 0.211 / 0.062 /
  0.059 Å; directions peak 17σ along [0 0 −1], 9.1× isotropic; Nb–O shell peaks 1.835 and
  2.015 Å, CN 6; no zone-boundary order.
- **Ordered antiphase** (Ti ±0.10 Å along c, q = (0, 0, ½)): Ti kurtosis −0.27, rms 0.115 / 0.059 /
  0.058 Å, as for disorder; `analyze_modulations`: (0, 0, ½) at 498× the median, > 8σ. Sr's
  strongest mode, a long wave at (0.125, 0, 0), 100× the median, is box-scale, not order.
- **Breathing Mo₃ triangles** (size sd 4 %): breathing +49 % (z = 22); shape at the clusters' own
  size as for independent atoms (edge spread −1 %, z = −0.6).
- **Rigid V₂ dimers**: bond-length spread −87 % (z = −37).
- **Jahn–Teller CuF₆** (4 × 1.93 + 2 × 2.25 Å, long axis at random): Cu–F shell CN 6, first peak
  alone 4; F sites prolate (rms 0.16 / 0.05 / 0.05 Å); F₆ clusters: more elongated along <100>
  than by chance (z = 12), split member sites named.
- **Coherent lowering** (every Ta 0.07 Å along [001]): ladder P4mm above the 0.026 Å noise floor,
  the card Pm-3m at 0.2 Å; `measure_symmetry_lowering` Ta 0.065 Å along [0 0 1].
- **Thermal motion only** (BaTiO₃ composition): every kurtosis within ±0.15; O map 7.2σ with peak
  2.3σ (anisotropic thermal motion); no zone-boundary or general-wavevector order; Pm-3m above the
  noise floor. The answer: no order.
