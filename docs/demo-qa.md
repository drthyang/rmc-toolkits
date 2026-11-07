# Demo QA for 0.8.0

The public Demo (the bundled GaTa₄Se₈ 250 K run, the **Demo** button) was run end to end in real
browsers before presenting 0.8.0: every page, every export, and the AI Copilot with a free local
model and a scripted one. The console and the network were captured throughout; main-thread stalls
were measured with long-task entries. This page records what was checked, what was fixed, and what
is known.

## Environment

| | |
|---|---|
| Sites | the deployed site (GitHub Pages, 0.8.0) and a local production build of this branch (`VITE_STATIC_MODE=true`, `vite preview`, the same `/rmc-toolkits/` base) |
| Browsers | Chromium 149 and WebKit 26.5 (Safari's engine), driven by Playwright; Chromium on the Apple GPU (ANGLE Metal) |
| Sizes | desktop 1440 × 900; iPhone Pro 402 × 874 (touch, ×3); iPad Pro 1024 × 1366 (touch, ×2, WebKit) |
| Copilot | gemma4:26b on a local Ollama 0.40 (free, nothing sent to a paid API), one question at a time with a 6-minute limit and 30 minutes for the whole Copilot QA (it took about 21); a scripted LM-Studio-like server (fixed replies, no model) for the context meter, the change card on every size and the refusal |
| Data | the bundled public Demo run only |

Timings are from an Apple-silicon laptop. Under a software renderer (headless Chromium without a
GPU) every 3D frame is a 70–130 ms task; with a GPU none of the pages comes close.

## Checklist

Pass = works with no console error or warning; fixed = broken before this QA, fixed (commit
below); known = documented under *Known issues*.

### Load and Dashboard

| Check | Result |
|---|---|
| Page loads, no 404, no console message (both sites, every size) | pass |
| Demo loads the run; Model information, Detected SG, six plots, χ² history | pass |
| Detected SG ladder: P1, Cm, F-43m picked; the CIF button follows the pick | pass |
| ⤓ CIF downloads `GTS_250K_F-43m.cif` | pass |
| Save all figures (zip), a plot's Save as SVG | pass |
| Phone and iPad layout: cards stack, nav scrolls sideways | pass |

### Atomic Density

| Check | Result |
|---|---|
| KDE slice, Slab In Cell and the folded 3D cell draw | pass |
| Element (Ta, All), Grid (220), Slice slider by keyboard | pass |
| Normal menu: a, c, Plane (hkl) with h k l boxes | fixed in Safari: the menu opened again after every pick |
| KDE Save as high-res PNG, 3D Save as PNG | pass |

### Bond Geometry

| Check | Result |
|---|---|
| Compute Se–Ta–Se 2.00–3.00 Å: 236,431 angles, CN 5.95, Ta–Se 2.637 Å | pass |
| Tabs: sin-corrected, density, g(r) | pass |
| Angle range drag 80–100°: readout, band, Clear (mouse; touch swipe on phones) | pass |
| 3D bonds; the range's angle arcs in the folded cell | pass |
| Phone with the Copilot sheet open: the angle plot's height | fixed: the plot was ~60 px tall, its y label cut |

### PCA Ellipsoid

| Check | Result |
|---|---|
| Site picker (52 sites), statistics table, 3D density and site ellipsoids | pass |
| Isosurface on (clears wireframe and shell), Crystal frame | pass |
| 3D panel legend | fixed: it read "50% ellipsoid" with no ellipsoid drawn and never named the isosurface |
| Save as PNG | pass |

### Displacement Directions

| Check | Result |
|---|---|
| Sphere, readouts, site ellipsoids | pass |
| Site change, Resolution Auto (ν = 2 for a Se site) | pass |
| Per-cell readout: hover on desktop, tap on touch | pass |
| Safari console | fixed: "ResizeObserver loop completed with undelivered notifications" as a page error |
| Axis views on phones and iPad | fixed: three thin strips, each sphere ~40 px across; now side by side |
| Save as PNG | pass |

### AI Copilot

| Check | Result |
|---|---|
| Connects to local Ollama from the local build; model list; context meter neutral (9–10 %) | pass |
| `set_bond_geometry` → approval card → **Apply**: the page takes the inputs and computes | pass |
| **Undo**: the earlier inputs are back ("Undone: the earlier inputs are back.") | pass |
| **Decline**: nothing changes ("Declined: nothing changed.") | pass |
| Stale: the bin width edited by hand while the card waits, then Apply → "Not applied: the page changed while this waited for you." | pass |
| The Workbench's own question (which Ta–Se window), answered with the default | pass |
| Numbers by reference: no raw `{{…}}`, no "number missing" | pass |
| A number the model typed before its reference ("2.64 Å 2.64", "CN 6 6") | fixed (display only) |
| Short answers: the conclusion first, the rest behind Show details | pass |
| Outcome note "a check below…" with the checks drawn above it | fixed (display only) |
| An answer check misreads "Ta–Se … coordination number of 6" | known |
| Context meter amber (82 %, a 25,000-token loaded context) and the refusal past 100 % (18,000) with its plain message | pass (scripted server) |
| Change card, Apply, Undo on phone and iPad (bottom sheet) | pass (scripted server) |
| Deployed site with a local Ollama | known: needs the origin allowed; the settings drawer's Test says how |

## Timings

Chromium, local build unless noted; the longest main-thread task in each step.

| Step | Time |
|---|---|
| Page load (network idle) | 0.8 s local, 1.1 s deployed |
| Demo button → first plot drawn | 0.27–0.42 s desktop, 0.9–0.95 s phone; longest task 120–180 ms |
| First visit: Atomic Density / PCA Ellipsoid / Displacement Directions | longest task ≤ 240 / 380 / 460 ms (the deployed site's first visits the slowest) |
| Bond Geometry Compute (236,431 angles, in its worker) | 0.14–1.0 s wall, no long task |
| Angle range drag | one task ≤ 200 ms |
| KDE element, normal, grid and slice changes | no task over 250 ms |
| Main thread kept busy on the Dashboard after visiting every page (per 5 s) | ~570 ms before, ~100 ms after the fix |
| Copilot, gemma4:26b: change request → approval card | 7–85 s (the model's thinking) |
| Copilot, gemma4:26b: a whole turn (clarifier, analysis, answer) | 60–100 s |

No step stalled the main thread for anything near 1 s. WebKit has no long-task API; its runs were
checked by wall clock and console only.

## Fixes made

| Commit | Fix |
|---|---|
| 1c58f87 | Displacement Directions: the sphere resizes on the next frame (Safari's ResizeObserver loop error) |
| 67014d1 | PCA Ellipsoid: the legend names what the panel draws |
| 9ff978a | Atomic Density: the Normal menu closes on a pick in Safari (it sat inside a `<label>`) |
| 4d7bc27 | Copilot: the Outcome sits above the checks its note points to |
| 0387c42 | Copilot: a number typed before its own reference shows once |
| 2e0d9e9 | Displacement Directions: the axis views side by side in a wide panel |
| b40fefd | Bond Geometry: the angle plot keeps a readable height on phones |
| 05be220 | 3D panels draw nothing while their page is hidden |

Nothing the model reads changed: the two Copilot fixes are in the dock's display (the verification,
the judged answer, the history and the validation protocol 2026-11-01.1 are as they were).

## Known issues and workarounds

- **The Copilot on the deployed site with a local Ollama.** Ollama answers only the origins it
  allows, and Chrome asks the user to let a public site reach the local network. Until both, the
  Copilot shows "Connect a model to start" and the console logs the blocked request. Workaround:
  `launchctl setenv OLLAMA_ORIGINS "https://drthyang.github.io"`, restart Ollama, and allow
  Chrome's prompt (QuickStart, *4. AI Copilot*); or present the Copilot from a local build on
  `http://localhost`, which Ollama allows by default.
- **An answer check misreads a pair's coordination.** For "the Ta–Se first shell … coordination
  number of 6", the check takes the element nearest the number (Se) as the centre and reports "The
  answer gives Se a coordination of 6, but … 3 Ta per Se", so a right answer reads *Not verified*.
  The checks are part of the validation protocol and stay as they are for 0.8.0. Asking for the
  count as "Se per Ta" does not avoid it: in a later check the answer said "6 Se atoms per Ta", but
  the Outcome's own summary line still read "coordination is 6" and drew the same Check. For a
  demonstration, ask an interpretation question that does not quote a pair's coordination (the
  [demo script](demo-script.md) uses "Has this run converged, and how good is the fit?").
- **A number shown twice after its unit.** Fixed in 92b81f1 (display only). Fix 0387c42 covered a
  number typed just before its own reference, or with one of a few units between; where a label sat
  between them it still showed twice ("3.25 Se per Ta 3.25", "127 deg² 127"; gemma4:26b, local
  build). The dock now shows it once when a short label (up to three words, no digits or clause
  punctuation) sits between, and only when the typed number matches the engine's value.
- **gemma4:26b's own habits.** It can think for over a minute before a change card, may propose
  the change again after a stale refusal (Decline the second card), and sometimes runs one analysis
  more than asked. Asked whether the TaSe₆ octahedron is distorted, it once passed
  `window_confirmed` to the change tool, which refused the unknown argument (a *Not applied* card),
  and then answered from the analyses. None is an app error.
- **One 404 when connecting to LM Studio.** Server detection asks `/api/version` (Ollama's) first;
  LM Studio answers 404, which the console shows. Expected.
- **Idle 3D pages still draw every frame.** The page in view keeps its render loop (PCA Ellipsoid
  ~10 % of the main thread at rest). No stall; drawing on demand is for after the freeze.
- **Build notices.** Vite warns that the main chunk is over 500 kB (2.2 MB, 0.67 MB gzipped), and
  ESLint and Vite note that `baseline-browser-mapping`'s data is old. Neither shows in the browser;
  splitting the bundle or updating the dependency is for after the freeze.
- **iPhone Safari with a picked folder** (AGENTS.md, *Current known issues*) is outside the Demo and
  was not retested.
