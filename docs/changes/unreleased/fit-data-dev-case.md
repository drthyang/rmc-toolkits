### Copilot validation: a fit-data dev case (plan 4.6)

- **dev-fit-unconverged**, a dev-only case (`devOnly` in `validation/cases.js`; no test twin, no
  reported number), is the first case with fit files, so `assess_convergence`, `compare_fit` and
  the rmc-fit-assessment skill are exercised. Rock-salt MgO, 6³ cells, isotropic thermal motion
  only (no order, Fm-3m), seeded. Besides the `.rmc6f` (now with RMCProfile's run counters: 9.4
  accepted moves per atom) and the partials, its run folder holds RMCProfile's X-ray file family
  in the Demo's formats: `rmc_run_FQ1.csv`, its Fourier transform `rmc_run_FT_XFQ1.csv` and the
  χ² log `rmc_run-00.log`. `runFolder.js` passes them as run-folder entries (`caseFiles`), so the
  run module detects their kinds as for a user's run.
- The planted fit truth (`truth.fit`): the X_ray_(R)1 χ² is still falling after a fast drop (the
  plot's end looks flat; the watchdog and `assess_convergence` say improving, ln χ² −0.034 over
  the last tenth), and the measured G(r) holds a peak at 2.55 Å the configuration does not
  reproduce (`compare_fit`: 97 % of the G(r) misfit in 2.17–3.33 Å; in F(Q), the same
  measurement, a misfit at low Q).
- Scored on its structure like every case; the scoring and the verdict schema are unchanged. The
  fit truth is checked by the context, the fit tools and the answer check's convergence rule
  (`validation/__tests__/fitCase.test.js`, scripted model). Arms A and A+ score 1 on it (appended
  to `validation/baseline/arm-a-dev.json`); the other 27 dev records and the test set are
  unchanged.
- Pins updated: `arms.test.js`, `split.test.js` (a dev-only kind), `caseBlind.test.js` (the dev
  cases now offer every tool its sweep calls), `toolDataShapes.test.js` (the fit tools on the
  case), and `contextBudget.test.js`: arm B's largest dev case is now this one, 35,868 tokens
  needed, 3,100 over 32,768.
