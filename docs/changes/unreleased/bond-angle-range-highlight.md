### Bond Geometry: pick an angle range and see its triplets in 3D

- Drag across the angle plot (sin-corrected or density) to pick a range of angles, drawn as the
  Slab In Cell band. Drag the band or its edges to move it, by mouse or finger; double-click it or
  press Clear to remove it. The card header reads the range in one line ("85–95°").
- The folded cell draws every A–B–C triplet in the range as thicker tubes in its bond colours and
  dims the other bonds; the legend gives the count ("2,310 angles in 85–95°"). With no range,
  pointing at a bin (tapping it on a phone) previews that bin.
- The triplets are the ones the histogram counts: the same windows, periodic images and central
  atoms, walked by the engine's own loop (tested against an all-pairs search, bin by bin). The
  range moves in whole bins, so its count is the histogram's. Each triplet is drawn between the
  average sites, so a drawn angle is the sites' average one.
- A large cell keeps to the bonds it draws (those nearest its centre) and to at most 750
  triplets, and says how many it shows ("40 of 2,310 angles in 85–95° shown").
- The triplets are found in a worker of their own after Compute, so Compute is not slowed. Not
  yet for a server folder (a run opened in the browser only). Nothing the AI Copilot is sent
  changes.
