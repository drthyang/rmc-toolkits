### Bond Geometry: the partial g(r) is a tab of the angle card

- The partial g(r) window helper is now the third tab of the angle card (sin-corrected, the
  default, then density, then g(r)) instead of a card of its own; same curves, guides and window
  chip. The folded cell takes the right column at full height. The KPI rail stays in place on
  every tab.
- The AI Copilot's card ids are unchanged: a highlight on `bond-windows` (now the g(r) tab panel)
  brings that tab forward, and one on `bond-angles` brings the angle view back. Nothing the model
  is sent changes.
- The folded cell draws the bonds as shaded tubes in their bond-window colours (A–B blue, B–C
  orange, as in the legend and the g(r) guides) instead of thin translucent lines, and fades the
  clouds of elements outside the bonds, so the selected triplet stands out.
- The bond search is a linked-cell one: on a synthetic 8000-site P1 cell (50 Å) the bonds now
  appear in ~65 ms instead of a ~17 s freeze, with the same bonds (tested against the old search).
  A large cell draws at most 1,500 bonds, those nearest its centre, and the legend says how many of
  how many it shows. Draw calls are unchanged (one per window).
