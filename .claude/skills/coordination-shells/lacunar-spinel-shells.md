# Lacunar spinel shells (GaTa₄Se₈, GaNb₄Se₈)

Reference values for the Workbench's two reference compounds. Measured values come from the
Copilot's tools on the runs (2026-10-10 review); published values from the user's papers. Use
them to recognise a familiar case, never as another run's numbers.

## Structure in brief

AM₄X₈ lacunar spinels (F-43m at high temperature): M₄ tetrahedral clusters (Ta₄, Nb₄) and GaX₄
tetrahedra on an fcc arrangement. Each M sits in a distorted MX₆ octahedron whose six M–X bonds
split into two lengths, so the first M–X shell has two peaks.

## Measured on the Demo (GaTa₄Se₈, 250 K)

| Polyhedron | Pair | Shell window | Peaks in shell | CN | Notes |
|---|---|---|---|---|---|
| TaSe₆ | Ta–Se | ≈ 2.21–3.23 Å (find_gr_peaks, 2026-10-11); 2.33–2.96 Å in the review | 2.52, 2.74 Å | 6 (5.89 in 2.33–2.96 Å) | first-peak window 2.33–2.71 Å holds only 3.9: it cuts the shell |
| GaSe₄ | Ga–Se | 2.18–2.79 Å | — | 3.99 | mean 2.434 Å, sd 0.095 Å; angle 109.4° (sd 4.1°); low-r side shaped by the 2.2 Å constraint |
| Ta₄ | Ta–Ta | shell mean 3.00 Å (r²g) | g peak 2.93 Å | 3 | clusters 23 % less distorted than independent atoms (z = −64), breathing +27 % |

With the 2.3–2.7 Å window only 47 % of the Ta atoms share the most common Se count (3), spread
from 1 to 6; whole-shell windows give 83–96 %.

## Published (the user's papers)

- **GaTa₄Se₈:** Phys. Rev. Research 4, 033123 (2022). F-43m above T* = 50 K, P-42₁m below (c
  doubled; a = 10.3437, c = 20.6878 Å at 10 K). A local non-cubic Ta₄ distortion at all
  temperatures to 300 K: monoclinic (one mirror), elongated along <110>, 3 long (≈ 3.04–3.06 Å)
  and 3 short (≈ 2.95–2.97 Å) Ta–Ta bonds; intercluster Ta–Ta ≈ 4.28–4.30 Å; a double PDF peak
  near 3 Å. Small-box PDF fits, not RMC.
- **GaNb₄Se₈:** Phys. Rev. B 109, 144101 (2024). F-43m → P2₁3 at T_Q = 50 K (Nb₄ elongated along
  <111>, a Jahn–Teller X5 mode) → P2₁2₁2₁ at T_M = 29 K (cell doubled; Nb₄ elongated along <011>,
  2 long, 2 short and 2 intermediate bonds; intercluster dimers along <011>). X-ray only.

The Workbench's RMC runs show neither the 5 K cell doubling nor a cluster distortion beyond
independent atoms: say so when a question expects the published distortion.
