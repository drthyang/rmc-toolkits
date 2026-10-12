# Measured numbers from the 2026-10-10 review

The Workbench's two reference runs, as the Copilot's tools measured them during the referee-style
review (llm README, "A scientific review"). Use them to recognise a familiar pattern, never as
values for another run.

## Demo: GaTa₄Se₈, 250 K (`public/demo`, published data)

- **Data:** one x-ray measurement fitted in both spaces (F(Q) and its Fourier transform G(r)),
  so one dataset, not two.
- **Convergence:** still improving at the end of the history (ln χ² −0.11 over the last tenth;
  the spot check's slope −0.00052 per row). Not converged.
- **Sampling:** about 3.5 accepted moves per atom (1,560,869 generated, 183,014 accepted moves
  for 52,000 atoms).
- **Weights:** Ga carries about 5 % of the x-ray signal: trust only Ga's average site.
- **Constraint:** Ga–Se closest approach 2.2 Å; g(r) just above it is 8.7 % of the first peak, so
  the low-r side of the Ga–Se peak is constraint-shaped.
- **Fit:** G(r) Rw 0.024 in the range that holds 62 % of the misfit (2.34–5.28 Å); largest
  residuals at 2.53, 3.66, 3.82 and 4.9 Å.
- **Box-scale waves:** the strongest displacement modes are the box's longest waves; there is no
  zone-boundary (cell-doubling) order.

## 5 K GaNb₄Se₈ (`data/5K_try1`, unpublished, not committed)

- Still improving; its G(r) misfit sits at 10.8–15 Å. F(Q) was not checked in the review, so its
  "achieved" verdict was not earned.
- No cell doubling in the configuration and no cluster distortion beyond independent atoms,
  although diffraction reports P2₁2₁2₁ with a doubled cell below 29 K (Phys. Rev. B 109, 144101
  (2024)): the RMC box (10×10×10 of the cubic cell) did not develop it.
