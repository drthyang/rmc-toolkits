### Demo QA for 0.8.0: Safari and phone fixes, hidden 3D panels idle

- The public Demo was run end to end in Chromium and WebKit at desktop, iPhone and iPad sizes, with
  the Copilot on a local model and a scripted one; the checklist, timings and known issues are in
  [docs/demo-qa.md](demo-qa.md).
- Displacement Directions: Safari no longer reports "ResizeObserver loop completed with undelivered
  notifications" (the sphere resizes on the next frame). On phones and tablets the Axis views sit
  side by side in their wide row, instead of three thin strips.
- Atomic Density: in Safari the Normal menu closes on a pick (it opened again: the menu sat inside a
  `<label>`, which forwarded the tap to its first button).
- Bond Geometry: on a phone, with the Copilot's sheet open, the angle plot keeps a readable height
  (the stacked card's floor is 30rem, was 24rem).
- PCA Ellipsoid: the 3D panel's legend names what it draws — the ellipsoid only while it is drawn,
  and the isosurface while it is on.
- 3D panels draw nothing while their page is hidden: after visiting every page, the Dashboard kept
  the main thread busy ~570 ms in every 5 s; ~100 ms now.
- Copilot (display only; nothing the model reads or the verification judges changes): an answer's
  Outcome sits above the checks its note calls "below"; a number the model types before its own
  reference ("2.64 Å {{r1.…}}") shows once, the engine's value, instead of "2.64 Å 2.64".
