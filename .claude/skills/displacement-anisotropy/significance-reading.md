# Reading the significances

Every significance the direction and cluster tools report is a calibrated one-sided normal
deviate (σ), with tails floored at 1e-300; ≥ 3σ counts as significant. A local, uncorrected
z-score is never one of them.

| Quantity | What it tests | Null |
|---|---|---|
| `peak.significance_sigma` | the most populated sphere cell, against the expected count | exact Poisson tail, Šidák-corrected over all cells |
| the map significance | the whole map's departure from isotropy | Pearson X² against a gamma matched to its exact mean, variance and skewness (undefined below ~0.1 expected coincident pairs) |
| `anisotropy_significance_sigma` | the orientation tensor's departure from isotropy | Bingham χ²₅ |
| +u/−u asymmetry z | more displacements one way along an axis than the other | each antipodal cell pair Bin(T, ½) |

**Enhancement** is the density over isotropic (1 = isotropic). The sphere is tiled by hexagons and
exactly 12 pentagons, every count divided by its cell's solid angle.

**Mardia kurtosis** (b₂ − 15)/5 is rotation- and affine-invariant: 0 for a Gaussian cloud, negative
for a flat or split site, positive for heavy tails. Per-axis kurtosis means something only on axes
the PCA resolves (an eigenvalue gap above 3 standard errors).

**Cluster statistics** (`analyze_clusters`) count as a finding only when both the z and the effect
size clear their thresholds (z ≥ 3 and a relative change ≥ 5 %; T2 trigonality ≥ 0.02; axis
enhancement ≥ 1.1×).

## The Demo (GaTa₄Se₈, 250 K)

- Most anharmonic site: ref 23 (Ta), kurtosis 5.62; largest U_iso 0.0175 Å².
- Ga, pooled: no preferred direction (peak −0.89σ, map 0.56σ).
- Ta₄ clusters: 23 % less distorted than independent atoms (z = −64), breathing +27 %; <111>
  compression preferred (19 % against 14 %, z = 11); no <001> signal.
