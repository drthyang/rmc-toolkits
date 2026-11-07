### Copilot: check_ripples, Fourier termination ripples (the dock only)

- The user's concern (2026-11-01): RMC is a minimiser and fits whatever the data carry, Fourier
  termination ripples included, so ripples must not be read as physical properties. A new
  assessment tool, `check_ripples`, says per fitted dataset how large that risk is and what the
  ripples are, and reads a small peak, a shoulder or a split against them.
- How RMCProfile makes the model it compares, established on its own files, not assumed: an X-ray
  F(Q) fitted in real space is compared through the same rectangle-rule transform of the model F(Q)
  (the Demo, to 1.6e-5 of max |G|): matched. A neutron G(r) without `CONVOLVE` is compared with the
  configuration's untruncated histogram (three runs, to 6e-15): the box is asked to make the data's
  ripples, risk high. With `CONVOLVE :: Q_c` the model is the histogram convolved with the sinc of
  Q_c in D(r) = r·G(r) space, the baseline below `R_CUTOFF` (RMCProfile 6.9.0-rc.0, to 1.8e-14):
  matched only against data of that Q_max with no window, `unknown` when the run's files do not say.
  An X-ray G(r) is the histogram transformed on its QMIN–QMAX grid and back. The evidence is
  committed: `tests/fixtures/rmcprofile_model_gr.json`, RMCProfile's output on a synthetic box.
- The engine, `rmc_toolkits/ripples.py` (the source of truth), its port `src/ripples.js` with
  `src/ripplesRun.js` (a browser-loaded run's files) and Flask's `/api/ripples` (a server folder):
  the data's transform (the stog input or Auto StoG provenance that made the file, `CONVOLVE`, or
  the caller's `qmax` / `window`), the r-resolution (Nyquist step π/Q_max, the kernel's FWHM on the
  run's own grid and window, 3.79/Q_max for a rectangle), the strong peaks as Gaussians fitted
  through the transform and their predicted ripples Δ = T[d] − T_[Q_min,∞)[d], the low-r level
  against them (a low r set to the baseline or Fourier-filtered is detected), and each feature:
  ripple-consistent, larger or smaller than the ripple there, a resolved peak, below the
  resolution, across datasets of different Q_max (a ripple moves, a feature stays), or at a
  predicted ripple in a partial g(r) of the box. A verdict, result first, and the next step: fit in
  Q or give `CONVOLVE`; extend Q_max, use a Lorch-windowed transform or the other radiation's data;
  a run from a new start only for a specific suspect feature.
- One golden, `ripples_fixture.json` (the Demo and seven synthetic runs), holds the browser and
  Flask to the Python engine: every class, risk and state exactly, numbers to 1e-7 (the strong
  peaks are an iterative fit). The golden comparison helpers take an `rtol` option for it.
- Known answers (`tests/test_ripples.py`): Gaussian peaks transformed exactly, cut at Q_max and
  back with RMCProfile's transform — the predicted ripples at the actual positions with amplitudes
  within 2 %; a planted shoulder above the resolution not called a ripple; every ripple
  ripple-consistent; two Q_max separating a ripple from a feature; Lorch suppressing the ripples
  as predicted; the low-r level the planted noise and ripple. On the Demo the model, transformed
  alike, shows 1.2× the predicted low-r ripples: an independent check on real data.
- Model-visible: the dock gains one tool (`dock/tools.json`, `first-request.txt`, the tool-results
  snapshots `check_ripples*.json`), +612 tokens on the dock's worst turn (headroom at 48k
  3,258 → 2,647). The validation's live arms leave it out by name (`LIVE_ARM_EXCLUDED`), so
  protocol 2026-11-01.1, arm B's and C's requests and the recorded 0.8.0 results are unchanged. No
  prompt, other tool description or skill changes.
