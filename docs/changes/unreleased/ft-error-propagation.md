### Dashboard: the G(r) error ribbon, F(Q)'s σ propagated through RMCProfile's transform (protocol 2026-10-24.4)

- New module `rmc_toolkits/ft_uncertainty.py` (exported; the browser port is
  `src/ftUncertainty.js` with its run glue `src/ftUncertaintyRun.js`). It rebuilds the transform
  RMCProfile writes as `<stem>_FT_XFQ<n>.csv`: `X_ray_exp_renorm` = Σ ΔQ·Q·(F + c)·sin(Qr) over the
  run's F(Q) grid, a rectangle rule with no 2/π and no ρ₀ (Keen's D(r) × π/2), no Lorch, c the run's
  fitted offset. The run's `F(Q)_Expt` is the `.dat`'s input file renormalised by the fitted scale
  s and, with `CONVOLVE`, box-convolved at the PDFpartials' last r, minus c. Two reproduction tests
  establish the chain on the run's own files before anything is drawn, each to 10⁻⁶ of the
  function's maximum. The demo's G(r) reproduces to 1.5×10⁻¹⁴, and two local runs' F(Q) from their
  input files to 5×10⁻¹⁵.
- σ comes only from the input file's own third column, clean by `usable_sigma`'s rule, times s. It
  is propagated exactly: σ(r) and the full covariance ½[P(r − r′) − P(r + r′)], with the correlation
  between r points and the number of independent points. The box convolution is left out of the
  propagation and checked exactly to change σ(r) by under 1 % (0.19 % at r ≤ 20 Å on the demo's
  grids).
- The Dashboard's xPDF card draws the ribbon (InteractivePlot's opt-in `ribbon`: a band under the
  experiment, in the legend, ± σ in the readout). When the G(r) file has its own σ column, that one
  is drawn and the propagated one is a toggle. The caption says the r points are correlated, and
  its ? help gives the transform, both reproduction residuals and the tolerance. With no σ it says
  why: on the demo, "F(Q) carries no uncertainties: scale_ft_rmc.fq (named in GTS_250K.dat) is not
  in the run folder". A failure is also a Problems-section warning, and Live Data re-reads the
  ribbon with its files.
- Flask's `/api/plot/data` sends the same payload as `uncertainty` on an xPDF plot (cached on the
  transform and its companions; the input read only inside the data roots). A column naming an
  uncertainty is the experimental series' `sigma`. One golden holds both runtimes to 10⁻⁹:
  `ft_uncertainty_fixture.json` (`tests/generate_ft_uncertainty_fixture.py`). Known answers cover a
  brute-force Jacobian, Monte Carlo, white noise, the Nyquist grid, MATERIA's pystog golden adapted
  to RMCProfile's convention, and StoG (data-gated).
- `compare_fit` calls an X-ray G(r)'s residual feature above noise only against σ(r) (≥ 3 σ). Its
  `noise` block gives the σ's source, the residual's rms in σ, the reproduction residual and
  tolerance, and the correlation. A G(r) without σ(r) says so and calls nothing above noise.
  Model-visible: its description (+56 tokens on the worst turns; 48k headroom 5,200 dock and
  7,999 arm B) and the Demo's G(r) result (new snapshot `compare_fit.xpdf.json`). Arms A and A+ on
  dev are unchanged. Approved by the manager on the user's delegation, 2026-10-24
  (`validation/findings/decisions/2026-10-24-ft-error-propagation.json`).
