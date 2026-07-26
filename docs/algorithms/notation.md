# Notation and conventions

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

Shared symbols, function conventions, coordinate frames, and the meaning of *reference-grade* vs *visualization-grade* as used throughout this reference.

## Contents

- [Symbols, collisions, and conventions](#symbols-collisions-and-conventions)
  - [1. Consolidated symbol table](#1-consolidated-symbol-table)
  - [2. Symbol collisions: read this before trusting a single letter](#2-symbol-collisions-read-this-before-trusting-a-single-letter)
  - [3. Conventions](#3-conventions)

---

## Symbols, collisions, and conventions

The Algorithms reference is fourteen self-contained sections, spread over six pages and written by
different passes over the same codebase. Each section carries its own local notation block; this
page is the union of those blocks, reconciled. **Where two sections use the same letter for
different things, both meanings are listed separately with their scope named.** Nothing is silently
merged — see
[Symbol collisions](#2-symbol-collisions-read-this-before-trusting-a-single-letter) for the index.

Section references use the numbers below.

| § | Section |
|---|---|
| §1 | Auto StoG — Step 0: inputs, composition constants, and data conditioning |
| §2 | Auto StoG — the transform layer: Fourier pair, Lorch, filter, Keen conversions |
| §3 | Auto StoG — the auto-scaling engine: level sweep, closed-form fit, self-consistency |
| §4 | Auto StoG — page workflow, outputs, and the written file family |
| §5 | Run Dashboard — file detection, parsing, and fit-quality metrics |
| §6 | Plot rendering, interaction, and figure export |
| §7 | Model summary and the Detected SG symmetry finder |
| §8 | Structure page — KDE density slices |
| §9 | Structure page — Slab In Cell projection and the 3D unit-cell view |
| §10 | PCA Ellipsoid page — displacement clouds, ADP tensors, and the separable 3D KDE |
| §11 | AI Assistant — run context construction, pair correlations, convergence heuristics |
| §12 | Principal axes in the crystallographic frame |
| §13 | Displacement Directions — the orientation-distribution engine |
| §14 | Displacement Directions — the sphere view, axis views, and site picker |

---

### 1. Consolidated symbol table

Ordered by topic: **scattering functions → structure and lattice → statistics, KDE and fit quality
→ directions and orientation → rendering geometry**. "Units" is the unit the *code* carries, not
the unit a textbook would use, wherever the two differ.

#### 1a. Scattering functions and the Auto StoG fit (§1–§4)

| Symbol | Meaning | Units | Where used |
|---|---|---|---|
| $Q$ | momentum transfer, $4\pi\sin\theta/\lambda$ | Å⁻¹ | §1–§5 |
| $Q_0$ | the **first** point of the supplied/cropped $Q$ grid, `q[0]` — lower limit of the numerical integral | Å⁻¹ | §2, §3 |
| $Q_{N-1}$ | the **last** point of the cropped grid, `q[-1]` — truncation edge and Lorch width. §2 calls this $Q_\mathrm{max}$; §3 renames it to keep it apart from `config.qmax` | Å⁻¹ | §2, §3 |
| $Q_\mathrm{min}, Q_\mathrm{max}$ | the requested crop window from the config (`config.qmax` is **not** $Q_{N-1}$) | Å⁻¹ | §1, §3 |
| $r$ | real-space pair separation | Å | §1–§4, §11 |
| $r_0$ | closest interatomic approach (first-shell onset) | Å | §1–§4 |
| $r_\mathrm{cut}$ | Fourier-filter cutoff radius | Å | §1, §2 |
| $r_\mathrm{fit,lo}$, $r_\mathrm{fit,hi}$ | bounds of the C2 (low-$r$) fit window | Å | §1, §3 |
| $S(Q)$ | Faber–Ziman normalized total-scattering structure factor, $\to 1$ at high $Q$ (Keen Eq. 19) | dimensionless | §1–§4 |
| $S_\mathrm{meas}(Q)$ | the measured, mis-scaled structure factor as read from file | dimensionless | §3, §4 |
| $S_\mathrm{corr}(Q) = a S_\mathrm{meas}(Q) + b$ | the corrected structure factor | dimensionless | §1, §3, §4 |
| $\sigma(Q)$ | per-point uncertainty of $S(Q)$ (optional third data column) | same as $S$ | §1, §4 |
| $F(Q) = Q[S(Q)-1]$ | PDF-community / pystog $F(Q)$; the transform argument | Å⁻¹ | §2, §3, §4 |
| $F_K(Q) = \langle b\rangle^2[S(Q)-1]$ | **Keen's** $F(Q)$ (Eq. 9); written to `*_rmc.fq` | barn | §2, §3, §4 |
| $g(r)$ | pair distribution function; $\equiv 0$ below $r_0$, $\to 1$ at large $r$ | dimensionless | §2, §4, §11 |
| $G_\mathrm{PDF}(r) = 4\pi\rho_0 r[g(r)-1]$ | PDFFIT-style $G(r)$, the transform intermediate | Å⁻² | §2, §3, §4 |
| $G_K(r) = \langle b\rangle^2[g(r)-1]$ | **Keen's** $G(r)$ (Eq. 10/16); flat at $-\langle b\rangle^2$ below $r_0$ | barn | §2, §3, §4 |
| $D(r) = 4\pi\rho_0 r\,G_K(r)$ | Keen's $D(r)$ (Eq. 29) | barn·Å⁻² | §2, §3, §4 |
| $c_i$ | atomic fraction of element $i$, $\sum_i c_i = 1$ | — | §1, §2, §3 |
| $b_i$ | bound coherent neutron scattering length of element $i$ (Sears/NIST) | fm | §1, §2, §3 |
| $\langle b\rangle^2 = (\sum_i c_i b_i)^2$ | the classic `stog.inp` "Faber–Ziman coefficient"; sets the normalization and the low-$r$ level | barn (stored fm², divided by 100) | §1–§4 |
| $\langle b^2\rangle = \sum_i c_i b_i^2$ | Laue/self term; sets the $Q\to0$ limit (Keen Eq. 14). **A different number** whenever there is >1 element | barn | §1–§4 |
| $\rho_0$ | atomic number density | Å⁻³ | §1–§4 |
| $\rho_m$ | **mass** density (§1 only — see the collision with the orientation $\rho_m$) | g/cm³ | §1 |
| $M$ | **molar mass** of one formula unit (§1 only) | g/mol | §1 |
| $n$ | atoms per formula unit (§1 only) | — | §1 |
| $a$ | fitted multiplicative scale of the affine correction | dimensionless | §1–§4 |
| $b$ | fitted additive offset of the affine correction | dimensionless | §1–§4 |
| $L$ | the measured high-$Q$ level (asymptote) of $S_\mathrm{meas}$, from the level sweep | dimensionless | §3, §4 |
| $a_L = \pi/Q_{N-1}$ | **Lorch length** — real-space resolution scale of the window. The source calls this `a` too (local to `low_q_correction_basis`); §2 writes $a_L$ to keep it apart from the fitted $a$ | Å | §2, §3 |
| $M(Q) = \sin(\pi Q/Q_{N-1})\big/(\pi Q/Q_{N-1})$ | the Lorch modification window (`lorch_window` / `lorchWindow`), $M(0)\equiv1$ | dimensionless | §2, §3 |
| $T[\cdot]$ | the forward sine transform $F(Q) \to G_\mathrm{PDF}(r)$ (carries the $2/\pi$; the backward direction does not) | operator | §2, §3 |
| $s_0 = F(Q_0)/Q_0 + 1$ | the measured $S(Q_0)$ (`s0`) | dimensionless | §2 |
| $s_0^\mathrm{target}$ | assumed $S(0)$ for the low-$Q$ extrapolation (`s0_target` / `effective_s0_target`) | dimensionless | §1, §2, §3 |
| $f_1, f_2$ | the two closed-form moment integrals of the low-$Q$ correction (`f1`, `f2`) | Å⁻³, Å⁻² (non-Lorch) | §2 |
| $K(r)$ | the numerical "offset kernel", the transform of $F \equiv Q$ (`g_one` / `gOne`) | Å⁻² | §2 |
| $\delta S(Q)$, $\delta(Q) = S_\mathrm{ft}(Q) - 1$ | the Fourier-filter correction term (`delta_sq`; classic `ft.dat` shifted), held fixed inside each fit iteration | dimensionless | §1, §3 |
| **C1** | the *high-$Q$ residual block* of the least-squares fit: $S_\mathrm{corr}\to1$ on $Q\ge Q_\mathrm{tail,min}$ | — | §1, §3, §4 |
| **C2** | the *low-$r$ residual block*: $g_\mathrm{corr}\to0$ on $[r_\mathrm{fit,lo}, r_\mathrm{fit,hi}]$ | — | §1, §3, §4 |
| $f_\mathrm{tail}$ | tail fraction of the $Q$ window that defines C1 (default 0.15) | — | §1 |
| $Q_\mathrm{tail,min}$ | start of the C1 window | Å⁻¹ | §1, §3 |
| $a_\mathrm{density}$ | amplitude implied by the low-$r$ density limit (the default criterion) | dimensionless | §1, §3, §4 |
| $a_\mathrm{fz}$ | amplitude implied by the $Q\to0$ Faber–Ziman limit (independent criterion) | dimensionless | §1, §3, §4 |
| concordance | $a_\mathrm{fz}/a_\mathrm{density}$ — equals 1 when the two criteria agree | — | §1, §3, §4 |
| $w_i$ | least-squares row weight ($1/\max(\sigma_i,10^{-12})$ normalized to unit mean on the tail; Huber IRLS re-weighting) | — | §3, §4 |

Unit convention throughout §1–§4: **1 barn = 100 fm²**. Scattering lengths are stored in fm;
every coefficient that reaches the engine is in barns (`faberZiman()` in `autoScale.js`,
`faber_ziman()` in `rmc_toolkits/scattering.py`, both dividing fm² by 100). For a normalized x-ray
$S(Q)$ the convention is $\langle b\rangle^2 = 1$ and
$\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$.

#### 1b. Structure, lattice and frames (§7–§10, §12)

| Symbol | Meaning | Units | Where used |
|---|---|---|---|
| $\mathsf{L}$ (rows $\mathbf L_i$) | $3\times3$ matrix whose **rows** are the RMC simulation-box (supercell) lattice vectors, from the `.rmc6f` `Lattice vectors` block | Å | §7, §8, §9, §10, §12, §13 |
| $\mathbf N = (N_1,N_2,N_3)$ | supercell tiling / repeat counts declared in the `.rmc6f` header | integers | §7–§10, §12 |
| $\mathsf{A}$ (= $A_\mathrm{conv}$), rows $\mathbf a,\mathbf b,\mathbf c$ | unit-cell (conventional-cell) matrix; row $i$ is $\mathbf L_i/\max(N_i,1)$ | Å | §7, §8, §9, §10, §13 |
| $M$ | §12's name for the **same** unit-cell matrix (rows $\mathbf a_1,\mathbf a_2,\mathbf a_3$) | Å | §12 |
| $\mathbf a,\mathbf b,\mathbf c$ ($\mathbf a_1,\mathbf a_2,\mathbf a_3$) | the three unit-cell edge vectors, in the shared Cartesian basis | Å | §9, §12, §13, §14 |
| $\ell_i = \lVert\mathbf a_i\rVert$, $\ell_\mathrm{max}$ | unit-cell edge lengths (`cellLengths`, `cellEdgeA`) and their maximum | Å | §9, §12 |
| $\mathbf b_i = \mathbf A_i/\max(\ell_1,\ell_2,\ell_3,10^{-9})$ | the **normalized** cell basis (longest edge exactly 1) used only by the three.js scene (`basis`) | normalized cell units | §9 |
| $\mathbf c$ | body centre of the *normalized* cell, $\tfrac12\sum_i\mathbf b_i$ (`center`) — never a cube corner | normalized cell units | §9 |
| $G = \mathsf A\mathsf A^{\mathsf T}$ | metric tensor of the conventional cell | Å² | §7 |
| $\mathbf f = (f_1,f_2,f_3)$ | atom coordinate **as stored in `.rmc6f`** — a fraction of the *box*, not of one unit cell | — | §8, §9, §10, §12 |
| $\mathbf x = (x_1,x_2,x_3)$ | atom folded into **one unit cell**, $x_i\in[0,1)$ | — | §7, §8, §9 |
| $(n_1,n_2,n_3)$, $\mathbf c_n$ | per-atom supercell (box-copy) indices from the `.rmc6f` | integers | §9, §10 |
| $\delta\mathbf f_n$ | folded, mean-subtracted fractional offset — a fractional *displacement* (supercell fraction) | — | §10 |
| $\mathbf u_n$, $\Delta\mathbf r_i$, $\mathbf u_m$ | Cartesian displacement of an atom from its site mean | Å | §10, §12, §13 |
| $R$ | integer rotation part of a symmetry operation, acting on fractional **columns** | — | §7 |
| $\mathbf t$ | fractional translation part of a symmetry operation | — | §7 |
| $\tau$ (`tol`, `symTol`) | Cartesian atomic-position tolerance for symmetry acceptance | Å | §7 |
| $\varrho$ (`residual`) | worst-site Cartesian mapping error of one symmetry operation | Å | §7 |
| $\epsilon_G$ (`metricTol`) | *relative* tolerance on the metric-preservation test (default $10^{-2}$) | — | §7 |
| $\mathbf h = (h,k,l)$ | user-chosen slice normal, entered as three numbers in *fractional index* space | — | §8, §9 |
| $\hat{\mathbf h}$ | $\mathbf h/\lVert\mathbf h\rVert_2$ | — | §8, §9 |
| $\hat{\mathbf u},\hat{\mathbf v}$ | in-plane axes, orthonormal **in fractional-component space** (not in Å) | — | §8, §9 |
| $\mathbf e_u,\mathbf e_v,\mathbf e_h$ | Å-space images of $\hat{\mathbf u},\hat{\mathbf v},\hat{\mathbf h}$, i.e. $\sum_i\hat u_i\mathbf A_i$ etc. | Å | §9 |
| $d = \hat{\mathbf h}\cdot\mathbf x$ | depth of an atom along the slice normal | — | §8, §9 |
| $[d_\mathrm{min},d_\mathrm{max}]$, $\Delta_d$ | projection range and depth span of the unit **cube** along $\hat{\mathbf h}$ | — | §8, §9 |
| $\tilde d = (d-d_\mathrm{min})/\Delta_d$ | normalized depth (`pointDepth`) | — | §9 |
| $z_c$, $\Delta z$ | the "Slice" and "Thickness" sliders — **fractions of $\Delta_d$**, not Å, not cell edges (see §3b of this page) | — | §8, §9 |
| $\delta = \Delta z\,\Delta_d$; $d_\mathrm{start},d_\mathrm{end}$ | band depth thickness and its cell-clipped edges | — | §9 |
| $d_{hkl}$ | interplanar spacing of the plane family $\mathbf h$; only appears in the (never-executed) Å conversion of §8 | Å | §8 |
| $\mathbf p_i = (u_i,v_i)$ | in-plane projection $(\mathbf x_i\!\cdot\!\hat{\mathbf u},\ \mathbf x_i\!\cdot\!\hat{\mathbf v})$ of an atom | — | §8, §9 |
| $\mathbf R = u\mathbf e_u + v\mathbf e_v + d\mathbf e_h$ | the atom's true Å position in the slab | Å | §9 |
| $\mathbf c_j$, $j=1\ldots8$ | the eight unit-cube corners (`CUBE_CORNERS`) | fractional | §9 |
| $\mathbf P,\mathbf Q$ | the two **Å-space** vectors handed to `makePlane` (§9 only — not $Q$ the momentum transfer) | Å | §9 |
| $\mathbf p,\mathbf q$ | their flattened 2D images | Å | §9 |

#### 1c. Statistics, KDE, and fit-quality metrics (§5, §8, §10, §11)

| Symbol | Meaning | Units | Where used |
|---|---|---|---|
| $N_\mathrm{src}$ | unique source atoms contributing to the slab (`slabCount`) | count | §8 |
| $N_\mathrm{img}$ | slab rows including periodic images | count | §8 |
| $n$ | fit points actually handed to the estimator (`fitCount`; capped at 6000) | count | §8 |
| $m$ | periodic-image margin used when tiling neighbour images | fractional | §8 |
| $\kappa = N_\mathrm{img}/N_\mathrm{src}$ | periodic renormalization factor of the 2-D KDE | — | §8 |
| $\mathbf C$ | $2\times2$ sample covariance of the 2-D fit points | fractional² | §8 |
| $f$ | KDE bandwidth factor (`bw` in §8, `factor` — Scott/Silverman — in §10) | dimensionless | §8, §10 |
| $\mathbf H$ | kernel bandwidth (covariance) matrix: $f^2\mathbf C$ in §8; $\mathrm{diag}(h_a^2)$ in §10 | fractional² (§8), Å² (§10) | §8, §10 |
| $h_a = f\sigma_a$ | per-axis kernel bandwidth of the separable 3-D KDE (`bandwidth`) | Å | §10 |
| $\rho(u,v)$ | estimated 2-D density on the slice plane | per unit fractional area | §8, §9 |
| $v_\mathrm{min},v_\mathrm{max}$ | min/max of the density grid actually drawn (after the $\log_{10}$ toggle) | as $\rho$ | §8, §9 |
| $G$ | grid points per axis of the 3-D KDE volume (`grid`) | count | §10 |
| $\mathbf U$ (`covariance`) | $3\times3$ displacement covariance = the ADP tensor, in **Cartesian** $x,y,z$ (not the CIF $U^{ij}$ basis) | Å² | §10, §12, §13 |
| $C$ | §12's name for the **same** $3\times3$ cloud covariance | Å² | §12 |
| $\lambda_1\ge\lambda_2\ge\lambda_3$ | eigenvalues of $\mathbf U$, descending, clamped $\ge0$ (`eigenvalues`) | Å² | §10, §12 |
| $\mathsf P$ | principal axes as **rows**, orthonormal, right-handed (`axes`) | — | §10, §12 |
| $\sigma_a = \sqrt{\lambda_a}$ | principal RMS amplitude (`rms`) | Å | §10, §12 |
| $U_\mathrm{iso} = U_\mathrm{eq} = \tfrac13\operatorname{tr}\mathbf U$ | equivalent isotropic displacement parameter (`uIso`) | Å² | §10 |
| $B_\mathrm{iso} = 8\pi^2 U_\mathrm{eq}$ | isotropic $B$ factor (`bIso`) | Å² | §10 |
| $k(p)$ | ellipsoid scale factor at enclosed probability $p$ — the $\chi^2$ quantile (`probabilityScale`) | — | §10 |
| $\tau(p)$ | isosurface density threshold enclosing mass fraction $p$ | Å⁻³ | §10 |
| $m_2^{(a)}, m_4^{(a)}$ | raw ($1/n$) 2nd and 4th moments of the cloud projected on principal axis $a$ | Å², Å⁴ | §10 |
| $\kappa_a$, $\kappa_i$ | per-axis **excess kurtosis** $m_4^{(a)}/D_a - 3$, where $D_a$ is the guarded form of $\bigl(m_2^{(a)}\bigr)^2$ — floored differently in the two engines (§10 Step 12) (`excessKurtosis`) | — | §10, §12 |
| `nonGaussianity` $= \tfrac13\sum_a\kappa_a$ | the summary non-Gaussianity readout | — | §10, §11 |
| $R$ | the number the dashboard labels "Rwp": $\sqrt{\sum(y^{(3)}-y^{(2)})^2 / \sum (y^{(2)})^2}$, **unweighted**, denominator = the RMC-*calculated* column | — | §5 |
| $\chi_r$ ("$\chi$", "$\chi^2$") | RMCProfile's goodness metric as parsed from the **last whitespace column** of the `.log` — the code never checks the header | — | §5, §11 |
| $y_i = \ln(\max(v_i,10^{-12}))$ | the plotted R-value series — natural log, floored at $10^{-12}$ on the **interactive** paths only (`plot_data()` in `app.py`, `browserData.js`); the floor is **absent** in `plots.py`'s static PNG, which applies a bare `np.log` | ln-units | §5, §11 |
| $N$ | number of points in the series (§11) / number of data rows in a file (§5) | count | §5, §11 |
| $n$, $L$ | recent-window length $\max(\min(N,5),\lceil 0.2N\rceil)$ — §11 writes $n$, §5 writes $L$ | points | §5, §11 |
| $m$ | OLS tail slope of the $\ln\chi$ series against local index (`recentSlope`, `recent_slope_per_step`) | ln-metric per log row | §5, §11 |
| $\Delta_\mathrm{win} = m\,(n-1)$ | window delta (`recent_window_delta`) — the actual classification variable; **not** the same field as $m$ | ln-metric per window | §5, §11 |

#### 1d. Directions and orientation (§12, §13, §14)

| Symbol | Meaning | Units | Where used |
|---|---|---|---|
| $\hat{\mathbf e}_i$, $i=1,2,3$ | the $i$-th **principal axis** of the displacement cloud; row $i$ of `axes`; unit Cartesian, sign- and handedness-canonicalised | — | §12, §13, §14 |
| $P$ | $3\times3$ matrix whose **rows** are $\hat{\mathbf e}_1,\hat{\mathbf e}_2,\hat{\mathbf e}_3$ | — | §12, §14 |
| $\mathbf e_0,\mathbf e_1,\mathbf e_2$ | orthonormal **display** frame Gram–Schmidt'd from $\mathbf a,\mathbf b$ (**0-indexed**; $\mathbf c$ is never read) — *not* the principal axes | — | §12 |
| $\mathbf r$ | an arbitrary Cartesian vector (a *displacement*, not a position) | Å | §12 |
| $\mathbf f$ | a displacement in **unit-cell fractional** coordinates | — | §12 |
| $\mathbf q$ | a displacement expressed in the PCA frame, $\mathbf q = P\mathbf r$ | Å | §12 |
| $[u\,v\,w]$ | crystallographic direction indices of a principal axis, $\propto M^{-\mathsf T}\hat{\mathbf e}_i$, rescaled to $\max_k\lvert\cdot\rvert=1$ (`crystalDirection`) | — | §12 |
| $\theta_{ij}$ | angle between $\hat{\mathbf e}_i$ and the cell edge $\mathbf a_j$ | degrees | §12 |
| $L$ | drawn rod length of an axis triad (§12 Step 9a) | Å | §12 |
| $a_i = \lVert\Delta\mathbf r_i\rVert$ | displacement **amplitude** of point $i$ | Å | §13, §14 |
| $\mathbf u_i = \Delta\mathbf r_i/\lVert\Delta\mathbf r_i\rVert$ | unit direction on $S^2$ | — | §13, §14 |
| $N_\mathrm{tot}$, $N$ | total points (`totalPoints`) and points surviving the amplitude cut (`usedPoints`) | count | §13 |
| $t$ | effective amplitude threshold; keep $i \iff a_i > \max(t, 10^{-9}\,\text{Å})$ | Å | §13 |
| $q$ | `min_amplitude_quantile` $\in[0,1)$ (§13 only — not $\mathbf q$ the PCA-frame displacement) | — | §13 |
| $\mathbf A$ | PCA **rotation** matrix used when `frame='pca'`; rows = the three PCA axes in descending-$\lambda$ order (§13 only — not the unit-cell matrix $\mathsf A$) | — | §13, §14 |
| $\nu$ | geodesic **frequency** of the tiling; engine range $[1,64]$, UI max 24, UI default 10 | integer | §13, §14 |
| $C = 10\nu^2+2$ | number of Goldberg cells | count | §13, §14 |
| $\varphi = (1+\sqrt5)/2$ | golden ratio, generating the 12 icosahedron vertices | — | §13 |
| $V$, $F_\triangle = 20\nu^2$, $E_\triangle = 30\nu^2$ | geodesic vertex / triangle / edge counts | count | §13 |
| $\chi$ | **Euler characteristic** $V-E+F=2$ (§13 §4.2 only) | — | §13 |
| $\mathbf c_m$ | cell centre of cell $m$ (a geodesic vertex), unit vector | — | §13, §14 |
| $\Omega_m$ | exact solid angle of cell $m$ (spread $-36\%$ to $+18\%$ about $4\pi/C$ at the UI default $\nu=10$) | sr | §13, §14 |
| $d_m \in \{5,6\}$ | degree (neighbour count) of cell $m$ | — | §13 |
| $w_i$ | per-point histogram weight: $1$ (`count`), $a_i$ (`amplitude`), $a_i^2$ (`amplitude2`) | —, Å, Å² | §13 |
| $M_m$ | weighted mass accumulated in cell $m$ (`mass`) | as $w$ | §13 |
| $n_m$ | **raw** count in cell $m$ (`counts`) — never smoothed | count | §13 |
| $S_m$ | per-cell amplitude sum $\sum_i a_i$ | Å | §13 |
| $\rho_m = M_m/\big((\sum_{m'}M_{m'})\Omega_m\big)$ | direction density over solid angle; $\sum_m\rho_m\Omega_m = 1$ | sr⁻¹ | §13, §14 |
| $\mathcal E_m = 4\pi\rho_m$ | **enhancement** — the displayed scalar; exactly 1 everywhere for an isotropic site | — | §13, §14 |
| $e_m = N\Omega_m/4\pi$ | Poisson expectation under the isotropic null (`expected`) | count | §13 |
| $z_m = (n_m-e_m)/\sqrt{\max(e_m,10^{-12})}$ | per-cell Poisson $z$-score (`zScore`) — from raw counts, before smoothing | σ | §13, §14 |
| `significance` $=\sqrt{\tfrac1C\sum_m z_m^2}$ | whole-map RMS $z$; $\approx1$ means pure counting noise | σ | §13, §14 |
| $\alpha = 0.5$ | `SMOOTHING_ALPHA`, the diffusion coefficient of one neighbour-smoothing pass | — | §13 |
| $\langle a\rangle_m$, $A_c$ | per-cell mean amplitude (`cellMeanAmplitude`); §13 writes $\langle a\rangle_m$, §14 writes $A_c$ | Å | §13, §14 |
| $\langle a\rangle$, $\bar A$ | global mean amplitude over surviving points (`meanAmplitude`), never smoothed | Å | §13, §14 |
| $\phi_m$, $f_c$ | per-cell **relief** radial factor, $\mathrm{clamp}_{[0.3,2.2]}(1+\rho(A_c/\bar A - 1))$; §13 writes $\phi_m$, §14 writes $f_c$ | — | §13, §14 |
| $\rho$ | the **Amplitude-height (relief) slider**, $[0,1]$ step 0.05, default 0.5 (§14 only — not a density) | — | §14 |
| $r(v)$ | per-vertex radius = arithmetic mean of the $f_c$ of the cells incident on $v$ | — | §14 |
| $\mathcal A$ | antipodal (inversion) asymmetry $\tfrac1{2N}\sum_m\lvert n_m - n_{\mathrm{antipode}(m)}\rvert$ | — | §13, §14 |
| $\mathbf T$ | **orientation tensor** $\sum_i w_i\mathbf u_i\mathbf u_i^{\mathsf T}/\sum_i w_i$; $\operatorname{tr}\mathbf T = 1$; computed from points, not bins | — | §13 |
| $\lambda_1\ge\lambda_2\ge\lambda_3$ (of $\mathbf T$) | orientation eigenvalues, **dimensionless**, summing to 1 (`orientationEigenvalues`) — not the ADP eigenvalues | — | §13 |
| `orientationAnisotropy` $=3\lambda_1-1$ | Woodcock-style scalar; 0 uniform, 2 single axis | — | §13, §14 |
| $t(v)$ | colormap coordinate $\in[0,1]$ from `colorCoordinate` | — | §14 |
| $p = \mathrm{clip}_{[0,1]}(1/v_\mathrm{max})$ | colour coordinate of the **isotropic level** — the contrast pivot (§14 only) | — | §14 |
| $\gamma$ | the **Contrast** slider, range 0.5–3.0 step 0.1, default 1 | — | §14 |
| $v_\mathrm{max}$ | $\max_c \mathcal E_c$ over the tiling (`result.vmax`); the colour range is $[0,v_\mathrm{max}]$, **not** symmetric about 1 | — | §14 |

#### 1e. Rendering geometry (§6, §9)

| Symbol | Meaning | Units | Where used |
|---|---|---|---|
| $x_i,y_i$ | the $i$-th data pair of one series, in the file's own units | file units | §6 |
| $[x_0,x_1]$, $[y_0,y_1]$ | the current x- and y-**domains** | data units | §6 |
| $p_x,p_y$ | coordinates in the SVG **user space** = viewBox units (*not* CSS pixels) | user units | §6 |
| $L,R,T,B$ | left/right/top/bottom margins of the viewBox (§6 only) | user units | §6 |
| $W,H$; $W_p = W-L-R$, $H_p = H-T-B$ | viewBox size and plot-area size | user units | §6 |
| $s$ | screen scale factor, CSS px per user unit | px/user unit | §6 |
| $W_c,H_c$ | canvas panel width/height (§6 Step 15) | CSS px | §6 |
| $W_{px},H_{px}$ | slab/KDE canvas width/height (§9) | CSS px | §9 |
| $k$ | canvas scale factor of the plane mapper, $\min\big((W_{px}-36)/\Delta X,\ (H_{px}-36)/\Delta Y\big)$ | CSS px/Å | §9 |
| $o_x,o_y$ | centring offsets of the fitted content | CSS px | §9 |
| $X_\mathrm{min}\ldots Y_\mathrm{max}$, $\Delta X,\Delta Y$ | bounds and spans of the flattened plane | Å | §9 |
| $R_{sph}$ | camera framing radius (§9 Step 7.5) | normalized cell units | §9 |
| $k$ | photoelectron wavenumber (EXAFS files, §5 only) | Å⁻¹ | §5 |
| ToF | neutron time of flight (x-axis of some dashboard files) | µs | §5 |

---

### 2. Symbol collisions: read this before trusting a single letter

Every entry below is a letter whose meaning **changes between sections**. Both meanings are in the
table above; this is the index.

| Letter | Meaning A | Meaning B (and further) |
|---|---|---|
| $a$ | fitted multiplicative scale (§1–§4) | Lorch length $a_L=\pi/Q_{N-1}$, called `a` in `low_q_correction_basis` (§2); unit-cell edge vector $\mathbf a$ (§9, §12, §13); amplitude $a_i$ (§13); mean amplitude $\langle a\rangle$ (§13, §14) |
| $b$ | fitted additive offset (§1–§4) | scattering length $b_i$ (§1–§3); unit-cell edge vector $\mathbf b$ (§9, §12); normalized scene basis $\mathbf b_i$ (§9) |
| $C$ | number of Goldberg cells $10\nu^2+2$ (§13, §14) | $2\times2$ KDE fit covariance (§8); $3\times3$ cloud covariance = $\mathbf U$ (§12) |
| $\mathbf A$ / $\mathsf A$ | unit-cell / conventional-cell matrix, rows $a,b,c$ (§7, §8, §10) | PCA rotation matrix, rows = PC axes (§13, §14); triangle corner $\mathbf A$ in the geodesic subdivision (§13); per-cell mean amplitude $A_c$ (§14) |
| $\mathcal A$ vs $\mathbf A$ | antipodal asymmetry $\mathcal A$ (§13) | — script vs bold is the only distinction; do not read them as the same object |
| $L$ | measured high-$Q$ level of $S_\mathrm{meas}$ (§3, §4) | supercell lattice matrix $\mathsf L$ (§7–§10, §12, §13); left margin (§6); recent-window length (§5); rod length (§12) |
| $M$ | molar mass, g/mol (§1) | Lorch window $M(Q)$ (§2, §3); unit-cell matrix, rows $a,b,c$ (§12); per-cell weighted mass $M_m$ (§13) |
| $m$ | periodic-image margin (§8) | OLS tail slope of $\ln\chi$ (§5, §11); cell index (§13); box-copy index (§12); raw moments $m_2,m_4$ (§10) |
| $n$ | atoms per formula unit (§1) | KDE fit-point count (§8); recent-window length (§11); cell degree/count subscripts (§13) |
| $N$ | supercell repeat counts $N_i$ (§7–§10, §12) | data rows in a file (§5); series length (§11); surviving direction count (§13) |
| $\rho$ | atomic number density $\rho_0$ (§1–§4) | 2-D KDE density $\rho(u,v)$ (§8, §9); orientation density $\rho_m$ (§13); the Amplitude-height slider $\rho$ (§14) |
| $\rho_m$ | **mass density**, g/cm³ (§1) | **orientation density of cell $m$**, sr⁻¹ (§13) — same subscripted glyph, unrelated quantities |
| $\sigma$ | per-point uncertainty $\sigma(Q)$ (§1, §4) | principal RMS amplitude $\sigma_a$ (§10, §12) |
| $\lambda$ | ADP eigenvalues, Å² (§10, §12) | orientation-tensor eigenvalues, dimensionless, summing to 1 (§13); neutron wavelength inside $Q=4\pi\sin\theta/\lambda$ (§2) |
| $\kappa$ | periodic renormalization $N_\mathrm{img}/N_\mathrm{src}$ (§8) | per-axis excess kurtosis $\kappa_a$ (§10, §12) |
| $\chi$ | RMCProfile goodness metric (§5, §11) | $\chi^2$ quantile for probability ellipsoids (§10); Euler characteristic (§13) |
| $\mathbf T$ / $T$ | orientation tensor (§13) | forward transform operator $T[\cdot]$ (§2, §3); top margin (§6) |
| $\tau$ | symmetry position tolerance, Å (§7) | isosurface density threshold $\tau(p)$, Å⁻³ (§10) |
| $G$ | metric tensor $\mathsf A\mathsf A^{\mathsf T}$, Å² (§7) | KDE grid points per axis (§10); the functions $G_\mathrm{PDF}(r)$, $G_K(r)$ (§2–§4) |
| $\mathbf P$ / $P$ | principal axes as rows (§10, §12, §14) | one of the two Å-space vectors handed to `makePlane` (§9) |
| $\mathbf Q$ / $Q$ | momentum transfer (§1–§5) | the second Å-space `makePlane` vector $\mathbf Q$ (§9) |
| $\mathbf q$ / $q$ | displacement in the PCA frame, Å (§12) | amplitude quantile $q\in[0,1)$ (§13); flattened 2-D plane vector $\mathbf q$ (§9); spherical circumcentre $\mathbf q$ (§13) |
| $\hat{\mathbf e}_i$ vs $\mathbf e_0,\mathbf e_1,\mathbf e_2$ | hatted, 1-indexed = **principal axes** (§12, §13) | unhatted, 0-indexed = **orthonormal display frame from the cell** (§12). "The $\mathbf e_0$–$\mathbf e_1$ face" is a crystal-frame face, never a PC plane |
| $\hat{\mathbf u}$ vs $\mathbf u$ | in-plane slice axis (§8, §9) | Cartesian displacement $\mathbf u_n$ (§10, §12); unit direction $\mathbf u_i$ (§13); the $u$ of $[u\,v\,w]$ (§12) |
| $h$ | Miller index / slice-normal component (§8, §9) | per-axis kernel bandwidth $h_a$, Å (§10); hexagon count (§13) |
| $k$ | Miller index (§8, §9) | canvas scale factor, px/Å (§9); photoelectron wavenumber, Å⁻¹ (§5); ellipsoid probability scale $k(p)$ (§10) |
| $R$ | unweighted R-factor (§5) | integer rotation part of a symmetry operation (§7); atom Å position $\mathbf R$ (§9); right margin (§6) |
| $p$ | enclosed probability of an ellipsoid / isosurface (§10) | colour coordinate of the isotropic level (§14); SVG user-space coordinate $p_x,p_y$ (§6); in-plane projection $\mathbf p_i$ (§8); pentagon count (§13) |
| $f$ | KDE bandwidth factor (§8, §10) | fractional coordinate $\mathbf f$ (§8–§10, §12); relief factor $f_c$ (§14); $f_\mathrm{tail}$ (§1); moment integrals $f_1,f_2$ (§2) |
| $w_i$ | least-squares row weight (§3, §4) | per-point orientation histogram weight (§13) |
| $z$ | slice-centre slider $z_c$ (§8, §9) | Poisson $z$-score $z_m$ (§13) |
| $\delta$ | Fourier-filter correction $\delta S(Q)$ (§1, §3) | band depth thickness $\delta=\Delta z\,\Delta_d$ (§9); fractional offset $\delta\mathbf f_n$ (§10) |
| $\phi$ vs $\varphi$ | relief factor $\phi_m$ (§13) | golden ratio $\varphi$ (§13) — same section, different glyph; do not transcribe one for the other |

**Two more that are not collisions but naming inconsistencies** between sections describing the
*same* quantity, kept here so a reader moving between sections is not surprised:

- The unit-cell matrix is $\mathsf A$ in §7/§8/§10 and $M$ in §12; the cloud covariance is
  $\mathbf U$ in §10/§13 and $C$ in §12; the relief factor is $\phi_m$ in §13 and $f_c$ in §14;
  the per-cell mean amplitude is $\langle a\rangle_m$ in §13 and $A_c$ in §14. Same code fields
  (`unitCell`, `covariance`, `reliefFactors`, `cellMeanAmplitude`) in every case.
- $Q_\mathrm{max}$ means "last supplied data point" in §2 and "the configured crop ceiling
  `config.qmax`" in §1/§3. §3 writes $Q_{N-1}$ for the former precisely because the two are
  routinely different (crop at 30 Å⁻¹, data end at 29.40 Å⁻¹ → the Lorch window is built from
  29.40).

---

### 3. Conventions

#### 3a. Function names for the scattering quantities

These follow `docs/STOG_SCALING_PLAN.md` §1.1–§1.2, which cites D. A. Keen, *J. Appl. Cryst.*
**34**, 172–177 (2001), [doi:10.1107/S0021889800019993](https://doi.org/10.1107/S0021889800019993).
Equation numbers below are that paper's, **as cited by the repo's own docs and code**.

| Name in this document | Definition | Keen Eq. | Units |
|---|---|---|---|
| $S(Q)$ | $S(Q)-1 = F_K(Q)/\langle b\rangle^2$; $\to1$ at high $Q$ | (19) | dimensionless |
| $F(Q)$ | $Q[S(Q)-1]$ — the **pystog / PDF-community** $F(Q)$, and the argument of every sine transform in the code | — (not Keen's $F$) | Å⁻¹ |
| $F_K(Q)$ | $\langle b\rangle^2[S(Q)-1]$ = $\sum_{ij}c_ic_jb_ib_j[A_{ij}(Q)-1]$ — **Keen's** $F(Q)$; written to `scale_ft_rmc.fq` | (9) | barn |
| $g(r)$ | pair distribution function; $\equiv0$ below $r_0$, $\to1$ | — | dimensionless |
| $G_\mathrm{PDF}(r)$ | $4\pi\rho_0 r[g(r)-1]$ — PDFFIT-style $G(r)$, the transform intermediate | — | Å⁻² |
| $G_K(r)$ | $\langle b\rangle^2[g(r)-1]$ = $\sum_{ij}c_ic_jb_ib_j[g_{ij}(r)-1]$ — **Keen's** $G(r)$; written to `scale_ft_rmc.gr` | (10), low-$r$ level (15)/(16) | barn |
| $D(r)$ | $4\pi\rho_0 r\,G_K(r)$; written to `scale_ft_rmc.dr` | (29) | barn·Å⁻² |

Limits used as fit constraints (`STOG_SCALING_PLAN.md` §1.2):

| Limit | Value | Keen Eq. |
|---|---|---|
| $S(Q\to\infty)$ | 1 | (21) |
| $F_K(Q\to0)$ | $-\langle b^2\rangle + \eta$ ($\eta$ = compressibility term, ignorable for dense solids) | (14) |
| $S(Q\to0)$ | $1 - \langle b^2\rangle/\langle b\rangle^2$ — zero only for a monatomic sample, slightly negative otherwise | (21) |
| $G_K(r<r_0)$ | $-\langle b\rangle^2$ (flat) | (15) |
| $D(r<r_0)$ | $-4\pi\rho_0\langle b\rangle^2 r$ (straight line through the origin) | (29)+(15) |

**$\langle b\rangle^2$ vs $\langle b^2\rangle$ — never interchange them.**
$\langle b\rangle^2 = (\sum_i c_i b_i)^2$ is the classic `stog.inp` "Faber–Ziman coefficient": it
sets the normalization and the flat low-$r$ level of $G_K$. $\langle b^2\rangle = \sum_i c_i b_i^2$
is the Laue/self term: it sets the $Q\to0$ limit. They are equal **only for a single element**.
For Mn₃Sn ($b_\mathrm{Mn}=-3.73$ fm, $b_\mathrm{Sn}=6.225$ fm) they differ by a factor of about 13.
Whichever unit system you give them in fixes the units of $F_K$, $G_K$ and $D$; the repo uses
barns everywhere (1 barn = 100 fm²).

**Scaling sign convention.** The correction is modelled as
$S_\mathrm{corr}(Q) = a\,S_\mathrm{meas}(Q) + b$ — **multiply**. The classic `stog_new` *divides*
(`S/yscale + yoffset`) while pystog multiplies; the repo's API uses the multiply form everywhere
and the `stog.inp` reader converts on the way in
(`STOG_SCALING_PLAN.md` §1.4, and §1 Step 4 "the divide-vs-multiply footgun").

**Fourier pair.** The pair is asymmetric — the $r\to Q$ direction carries no $2/\pi$:

$$G_\mathrm{PDF}(r) = \frac{2}{\pi}\int F(Q)\sin(Qr)\,dQ, \qquad
F(Q) = \int G_\mathrm{PDF}(r)\sin(Qr)\,dr .$$

Both are trapezoid-rule sums on the supplied grids — no resampling, no analytic quadrature.

#### 3b. Coordinates and frames

**Fractional vs Cartesian.** An `.rmc6f` file stores each atom as a fraction of the **simulation
box** ($\mathbf f$), together with per-atom supercell indices. Nearly every page's first step is to
*fold* that into a fraction of **one unit cell** ($\mathbf x \in [0,1)^3$) using
$\mathsf A_{i\cdot} = \mathsf L_{i\cdot}/N_i$. Confusing the two is the single most common way to be
off by a factor of $N_i$: $\mathbf f$ is a box fraction, $\mathbf x$ is a cell fraction, and a
*displacement* $\delta\mathbf f$ read out of the file is still a box fraction (§10 Step 1).

**One shared Cartesian basis.** Everything Cartesian in §10–§14 lives in the orthonormal ångström
basis in which the `.rmc6f` lattice vectors are written. Displacement clouds are mapped into it
through $\mathsf L$; the unit-cell vectors are derived from the same $\mathsf L$. That shared basis
is why no re-referencing is ever needed between the PCA axes and the crystal axes (§12, preamble).

**Row-vector vs column-vector — this differs between modules, so check the section:**

- **§7 (symmetry finder)** treats fractional coordinates as **column vectors**:
  $\mathbf x' = R\mathbf x + \mathbf t$. The lattice matrix $A$ stores lattice vectors as **rows**,
  so a Cartesian position is $\mathbf r = A^{\mathsf T}\mathbf x$ and
  $|\mathbf r|^2 = \mathbf x^{\mathsf T}G\mathbf x$ with $G = AA^{\mathsf T}$.
- **§12 (crystal-frame algebra, `pcaCrystalFrame.js`)** uses the same row-storage convention:
  $\mathbf r = M^{\mathsf T}\mathbf f$, inverse $\mathbf f = M^{-\mathsf T}\mathbf r$; and since the
  rows of $P$ are orthonormal principal axes, $\mathbf q = P\mathbf r$, $\mathbf r = P^{\mathsf T}\mathbf q$.
  Hence `fracToPca` $= PM^{\mathsf T}$ (units: Å per unit-cell fraction) and
  `pcaToFrac` $= M^{-\mathsf T}P^{\mathsf T}$.
  The crystallographic direction of a principal axis is $[u\,v\,w]_i \propto M^{-\mathsf T}\hat{\mathbf e}_i$
  — the **rows** of $M^{-\mathsf T}$ are the reciprocal vectors $a^\ast,b^\ast,c^\ast$ (no $2\pi$).
  Using $M^{-1}$ or $M$ instead gives, respectively, nothing meaningful and the Miller indices
  $(hkl)$ of the *perpendicular plane*. These coincide only for special metrics.
- **§11 (assistant context)** states its own contraction explicitly in the **row-vector** form
  $c_k = \sum_m \Delta f_m \mathsf A_{mk}$ — the same map, written without a transpose.
- **§13 (orientation engine)** rotates directions as $\mathbf u' = \mathbf A\mathbf u$ with the
  rows of $\mathbf A$ the PCA axes; to bring a PCA-frame vector back to crystal coordinates,
  $\mathbf v_\mathrm{cart} = \mathbf A^{\mathsf T}\mathbf v_\mathrm{pca}$.

In all of these the storage rule is the same — **matrices hold basis vectors as rows,
`m[row][col]`** — and only the choice of writing vectors as columns (§7, §12) or rows (§11) varies.
That is the difference to watch when moving a formula between sections.

**The cell-edge-fraction contract at the KDE API boundary (§8 Step 4).** `AGENTS.md` says
*"`z`/`dz` are cell-edge fractions at the API/slider boundary, converted to Ångström inside
`kde.py`."* As the code stands, **the first half holds only for the `a`/`b`/`c` presets and the
second half is stale**: `/api/kde/slice` passes `positions.fractional_positions` (not the Å array)
into `oriented_kde_slice()`, so the slab half-width, the KDE covariance, the evaluation grid and
the contour coordinates are all dimensionless. Ångströms enter only at draw time, when
`StructurePage.jsx` maps $\hat{\mathbf u},\hat{\mathbf v}$ through `unitCell.unitVectors`. The
honest statement of the contract is: **$z_c$ and $\Delta z$ are fractions of the projection range
$\Delta_d$ of the unit cube along the chosen normal**, coinciding with "fraction of a cell edge"
only for the presets (where $\Delta_d = 1$). The Å conversion
$\Delta z\,\Delta_d\lVert\mathbf h\rVert_2 d_{hkl}$ is derived in §8 but performed nowhere in the
app.

**The three display frames** (§10 Step 13, §12 Steps 5–9, §13 Step 2, §14 Step 10):

1. **PCA frame** — the orthonormal right-handed eigenframe $\hat{\mathbf e}_1,\hat{\mathbf e}_2,\hat{\mathbf e}_3$
   of the displacement covariance, rows of `axes`, sorted by descending eigenvalue and
   canonicalised (largest-magnitude component positive, then the third axis flipped if the frame
   came out left-handed — so the sign rule holds for PC1 and PC2 only). A canonicalised axis is an
   **undirected line**: a printed sign is a reproducibility convention, never a claim about which
   way atoms moved.
2. **Crystallographic frame** — the true cell edges $\mathbf a,\mathbf b,\mathbf c$, generally
   oblique. Comparisons like "the lobe points along $\mathbf a$" are statements about a *direction*;
   a Cartesian component is not a fractional index, and $[u\,v\,w]$ is the contraction that
   converts one to the other.
3. **Orthonormal display frame** $\mathbf e_0,\mathbf e_1,\mathbf e_2$ — Gram–Schmidt on
   $(\mathbf a,\mathbf b)$, completed by a cross product; $\mathbf c$ is never read. Used only for
   the crystal-mode shadow box, its wall projections, and the crystal-mode default camera. For an
   oblique cell only the $\mathbf e_0$–$\mathbf e_1$ face is a genuine crystallographic plane, and
   the shadow box is sized from the KDE half-widths, so **it is not a unit cell**.

#### 3c. The two runtime modes

The app ships in two runtimes: a **Flask backend** (`web_app/backend/app.py` serving the Python
engines in `rmc_toolkits/`) and a **static browser build** (GitHub Pages, no backend; the heavy
engines run in Web Workers under `web_app/frontend/src/workers/`, with the parsing, symmetry and
frame-algebra modules — `browserData.js` (§5), `symmetry.js` (§7), `pcaCrystalFrame.js` (§12),
`orientationSphere.js` (§14) — on the main thread at `web_app/frontend/src/`).

**Which engine runs is decided by the run *source*, not by the build mode.** This is the single
most-repeated trap in the document, and two sections prove it from the code:

- **§8 Step "Two code paths"** — `StructurePage.jsx` branches on
  `const isLocalStructure = Boolean(localRun)` and nothing else. `isStaticMode()` is used on that
  page only to emit the "Open a run folder" message. The **Demo** button is rendered
  unconditionally and sets `localRun`, so a Flask session that loads the bundled demo gets the
  **browser worker**, not `/api/kde/slice`.
- **§13 "Parity: Python engine vs JavaScript port"** and **§10 "Where the code lives"** — the
  switch lives in `requestPca(kind, params)` in `useSiteCloud.js`: if a local `.rmc6f` *text* is
  loaded (Demo or a picked folder) the shared worker answers **in either runtime**; otherwise
  axios hits `/api/pca/sites` | `/api/pca/orientation` | `/api/pca/kde`. A developer running the
  full Flask stack and opening a folder through the file picker is reading the **JavaScript**
  numbers.

Only a **typed backend directory** goes through the HTTP routes. The PCA statistics panel prints
`· browser` or `· server` so you can always tell which produced the numbers on screen (the JS
result carries `browserPcaKde: true`). §4 Step 0 gives the same split for Auto StoG; §11 notes
that the AI Assistant has **no** Python counterpart at all — the Flask backend is never involved in
an assistant request.

#### 3d. "Reference-grade" vs "visualization-grade"

The phrase is the app's own (§8: the in-app `InfoBadge` and the `local-density-note` both read
*"The Flask app uses SciPy KDE for reference-grade values"*). In this document:

- **Reference-grade** — the Python engine under `rmc_toolkits/`, reached through a Flask route
  with a server-side run **directory**. These are the "source of truth" implementations: §2/§3
  name `transforms.py` and `scaling.py`, §10 names `pca_kde.py`, §13 names `orientation.py`, §8
  names `kde.py` (`scipy.stats.gaussian_kde`). Float64 throughout, SciPy/LAPACK solvers,
  contourpy contours.
- **Visualization-grade** — the browser port, reached whenever a run was loaded as browser files.
  Same algorithm, deliberately different in the details: float32 arithmetic on the WebGPU KDE
  branch, a different pseudo-random subsample (mulberry32 vs PCG64) and therefore a **different
  bandwidth matrix**, a $e<-60$ kernel cutoff, a cruder per-cell contour tracer, a hand-rolled
  Jacobi eigensolver instead of LAPACK `eigh`, and a 12-point table + linear interpolation instead
  of `scipy.stats.chi2.ppf`.

**The practical rule.** Quote a number in a paper only if it came from the Python path — i.e. from
`rmc-autoscale`/the Flask API against a server-side directory, or from `rmc_toolkits` called
directly. When the browser produced it, quote it only where the document records a *measured*
cross-engine bound and your number is inside it. Those bounds are section-specific, not global:

| Path | Cross-engine agreement, as measured in the sections |
|---|---|
| Auto StoG level sweep (§4) | level rel. err. < 1e-9; window edges to 9 dp |
| Auto StoG $(a,b)$ (§4) | rel. err. < 1e-6; `iterations` exactly equal |
| Auto StoG $\rho_0$ self-consistency (§4) | **only ~1e-4 relative** — and §4 argues this is a genuine algorithmic divergence (the JS filter omits `s0Target` inside the loop), not pure float noise |
| PCA ellipsoids (§10) | round-off equality asserted at `rtol=1e-9, atol=1e-12` |
| Orientation frame (§13) | $\max\lvert\mathbf A_\mathrm{orient}-\mathbf A_\mathrm{ellipsoid}\rvert = 5\times10^{-16}$ — **but** the two suites are *parallel*, not shared-golden |
| Structure KDE (§8) | **no cross-runtime numerical test exists.** Bandwidth matrices genuinely differ (different subsample); GPU branch is float32; contour tracers differ; element-label case handling differs |

So: §8's browser density is a picture, not a measurement. §10/§13 browser numbers are
reproducible to round-off *for the quantities the parity tables list*. §4's browser $\rho_0$
carries a ~1e-4 honest bound. Anything not in a parity table has not been checked.

Two further limits that apply to *both* grades, so they are not fixed by using Python:
§5's "Rwp" is an **unweighted** R-factor normalized by the RMC-*calculated* column, not the
crystallographic $R_\mathrm{wp}$; and the R-value series is *whatever the last whitespace column of
the `.log` happens to hold* — the code never reads the header, and in the bundled demo that column
is `X_ray_(R)1`.

#### 3e. How code is cited

- **A behaviour** is cited as *file + function*: `sine_transform()` in
  [rmc_toolkits/transforms.py](../../rmc_toolkits/transforms.py). Function names are given in both
  languages where both exist (`orientation_histogram` / `orientationHistogram`), Python first.
- **A constant, a threshold, or a specific branch** is cited with a **line number**:
  `NEGLIGIBLE_AMPLITUDE = 1e-9` at [orientation.py:91](../../rmc_toolkits/orientation.py),
  `SMOOTHING_ALPHA = 0.5` at [orientation.py:95](../../rmc_toolkits/orientation.py). Line numbers
  appear **only** for those cases, because they are the citations that go stale first — treat a
  line number as a pointer into the source at the time of writing, not as a verified permanent
  address. Function names do not carry line numbers.
- **Tests** are cited by node id: `tests/test_parsers.py::test_rwp_uses_observed_series_as_denominator`,
  or by test name for JS suites.
- Paths are repo-relative, written as links from this file's own directory (`docs/algorithms/`), so
  code links begin `../../`, e.g. `../../rmc_toolkits/orientation.py`,
  `../../web_app/frontend/src/workers/orientation.js`.
