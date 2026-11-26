### Copilot: check_ripples, Fourier termination ripples (the dock only)

- The user's concern (2026-11-01): RMC is a minimiser and fits whatever the data carry, Fourier
  termination ripples included, so ripples must not be read as physical properties. A new
  assessment tool, `check_ripples`, says per fitted dataset, from the run's files and first
  principles, whether the box is asked to fit ripples, the data's Q range and window, and the
  r-resolution. It does not judge whether a given peak or shoulder is a ripple.
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
  the caller's `qmax` / `window`) and the r-resolution (Nyquist step π/Q_max, the kernel's FWHM on
  the run's own grid and window, 3.79/Q_max for a rectangle, the ripple period 2π/Q_max). The tool
  adds one line per real-space dataset: the ripple period, a peak's ripples near
  r₀ ± (n + ½)π/Q_max, features closer than the FWHM not resolved. A verdict (the risk) and the
  next step: fit in Q or give `CONVOLVE`; give `qmax`; to test a suspect feature, a higher Q_max, a
  Lorch-windowed transform or the other radiation's data (a ripple moves with Q_max); a run from a
  new start only for a specific suspect feature.
- Not shipped: the reading of a given feature (strong peaks fitted through the transform, their
  predicted ripples, the low-r level, feature classes, datasets of different Q_max compared, a
  partial g(r)) was not reliable enough on a known-truth test (synthetic fcc Ni) and is kept, with
  its tests, on the branch `wip/joint-fit`.
- One golden, `ripples_fixture.json` (the Demo and seven synthetic runs), holds the browser and
  Flask to the Python engine.
- Known answers: RMCProfile's model functions (`tests/test_ripples.py`); `resolution()` gives
  2π/Q_max and 3.79/Q_max within 1 %, and a Gaussian peak's ripple extrema lie within 0.05 Å of
  r₀ ± (n + ½)π/Q_max, n = 1–4, for σ 0.03–0.12 Å and Q_max 15–30 Å⁻¹, against an independent numpy
  transform (`tests/test_ripples_truth.py`).
- Model-visible: the dock gains one tool (`dock/tools.json`, the tool-result snapshot
  `check_ripples.json`), +374 tokens on the dock's worst turn (headroom at 48k 3,258 → 2,885). The
  validation's live arms leave it out by name (`LIVE_ARM_EXCLUDED`), so protocol 2026-11-01.1, arm
  B's and C's requests and the recorded 0.8.0 results are unchanged. No prompt, other tool
  description or skill changes.
