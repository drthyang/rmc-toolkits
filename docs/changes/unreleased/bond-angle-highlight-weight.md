### Bond Geometry: the angle highlight weighted by share, and rounder bonds

- Each triplet in the picked angle range (or the hovered bin) is drawn as strongly as its share:
  the part of its angles that fall in the range. A triplet almost wholly inside shows bright and
  thick; one that only reaches the range with the tail of its spread shows faint; under 2% it is
  not drawn (still counted). A bond used by several triplets shows the brightest. On the Demo
  (Se–Ta–Se), 70–80° and 100–110° now pick out 48 of the 96 Ta–Se bonds each instead of all 96;
  80–100° and 150–180° still hold every bond (each takes part in a cis angle near 88° and in one
  trans angle). The legend counts are unchanged; the "?" help says what the brightness means.
- A large cell keeps the 750 triplets with the largest shares (equal shares nearest the centre, as
  before).
- The folded cell's bonds are drawn like ball-and-stick sticks: lit with a soft highlight, round
  ends, edges that darken so crossing bonds stay apart, a little thinner, and each end stopping
  at the edge of its site's cloud instead of running through it. Same draw calls; on a software
  renderer the tubes take fewer segments and skip multisampling, so the view stays smooth.
