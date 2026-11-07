# Demo script: the Workbench in 5–10 minutes

A click path through the public Demo run (GaTa₄Se₈, 250 K, bundled with the app) for a live
presentation of 0.8.0. Each step has its time, what to click, and one line to say. The prompts are
the ones tested in the Demo QA ([demo-qa.md](demo-qa.md)) and in a check of this script on a local
build with gemma4:26b on Ollama.

Total: about 8 minutes, with 2 minutes of slack for questions or a slow model.

## Prep checklist (the day before, and 15 minutes before)

- [ ] **The model.** Ollama running with `gemma4:26b` pulled (`ollama list`).
- [ ] **The connection**, one of:
  - **Hosted site:** `launchctl setenv OLLAMA_ORIGINS "https://drthyang.github.io"`, then quit and
    reopen Ollama; open the site in Chrome, open the Copilot, and allow Chrome's prompt to let the
    site access your local network. The gear's **Test** must say it is connected.
  - **Local build** (no origin setting or prompt needed; also the offline fallback):

    ```bash
    cd web_app/frontend
    npm ci
    VITE_STATIC_MODE=true npm run build
    npm run preview        # http://localhost:4173/
    ```

- [ ] **Warm the model.** Ask one throwaway question (e.g. "Has this run converged?") and wait for
  the answer, then press **New conversation**. A cold gemma4:26b can take 7–85 s to reach a change
  card; a warm one answers a convergence question in about 40 s.
- [ ] **Settings drawer:** provider Ollama, model `gemma4:26b`, *Let the Copilot run analyses* on,
  *Changes: Ask first*.
- [ ] **Demo pre-loaded** on the Dashboard (press **Demo**), Copilot collapsed to its rail.
- [ ] **Browser zoom** so the room can read it: about 110 % on a 1440-wide screen. Check that the
  Copilot still docks beside the page, not as a bottom sheet (that needs a window at least 1080 px
  wide at the zoomed size).
- [ ] **Fallbacks open in a second tab:** this page and the screenshots in `assets/` (below).
- [ ] Notifications off, other tabs closed, laptop on power.

## The click path

| # | Time | Step |
|---|---|---|
| 0 | 0:30 | Open the site, press Demo |
| 1 | 1:30 | Dashboard: fits, χ², Detected SG ladder |
| 2 | 2:00 | Bond Geometry: Compute, g(r) tab, angle ranges, bonds and arcs |
| 3 | 1:00 | PCA Ellipsoid and Displacement Directions |
| 4 | 2:30 | Copilot: a change with Apply and Undo; an interpretation question |
| 5 | 0:30 | Validation, one sentence |
| 6 | 0:15 | Close |

### 0. Open (0:30)

- Open [drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/) (or the local
  build). Press **Demo**.
- **Say:** "This is an interpretation workbench for RMCProfile runs. It runs entirely in the
  browser; a real run folder is read locally and never uploaded. This is a public example run,
  GaTa₄Se₈ at 250 K, 52,000 atoms."

### 1. Dashboard (1:30)

- Point at the F(Q) and xPDF fits and their Rw badges (0.0457 and 0.0119). Press **Show** on the
  χ² history.
- Point at **Detected SG**: F-43m, and the *Space group vs. tolerance* ladder P1 → Cm → Cmm2 →
  P-42₁m → F-43m. Click **Cm**, then back on **F-43m**; the CIF button follows the pick.
- **Say:** "First question about any run: is the fit good, and has it converged? Then: what
  symmetry does the configuration actually keep? The average is F-43m; the ladder shows which
  symmetry survives as the tolerance tightens, and exports a CIF for any rung."

### 2. Bond Geometry (2:00)

- Open **Bond Geometry**. The inputs read Se–Ta–Se, 2.00–3.00 Å, 1°. Press **Compute**:
  236,431 angles, coordination 5.95 Se per Ta, mean Ta–Se 2.637 Å.
- Click the **g(r)** tab: the Ta–Se partial g(r) with the bond window marked. The first shell is
  split (peaks near 2.52 and 2.74 Å), and the 2.00–3.00 Å window holds all of it. Back to
  **sin-corrected**.
- Drag across **70–80°** on the distribution, then drag the band (or its edges) to **80–100°**.
  The readout counts the angles in the range; in the 3D cell the bonds brighten by the share of
  their angles in the range, and an arc at the central Ta marks each angle.
- **Say:** "An ideal octahedron has angles of 90° and 180°. Here the 90° peak splits into three
  (near 76°, 89° and 104°) and the 180° peak sits near 160°: the TaSe₆ octahedron is distorted,
  with a split Ta–Se shell. Dragging a range shows exactly which bonds make those angles, in the
  configuration itself, not in a model of it."

### 3. PCA Ellipsoid and Displacement Directions (1:00)

- Open **PCA Ellipsoid** and pick a Se site from *Site* (the screenshots use #6 Se).
  **Say:** "Per-site thermal ellipsoids from the PCA of each site's displacement cloud, with the
  3D density and each principal axis in the crystal frame."
- Open **Displacement Directions**, hover the sphere.
  **Say:** "Where atoms move rather than how far: directions binned on a sphere, each preferred
  direction with a calibrated significance, so noise isn't read as a hop direction."

### 4. Copilot (2:30)

Open the **AI Copilot** rail. If a conversation is left over, press **New conversation**.

**4a. A change it must ask for.** Type exactly:

> Set Bond Geometry to the Se–Ta–Se triplet with the Ta–Se window 2.4 to 2.9 Å.

- A **Bond Geometry change** card appears ("Needs your approval"): one line, the details behind
  "?", **Apply** / **Decline**. Nothing has changed yet.
- **Say:** "It has one tool that changes anything, and it asks first. Auto is opt-in and resets
  on reload."
- Press **Apply**: the page takes 2.40–2.90 Å and computes (215,252 angles, coordination 5.66).
  Press **Undo**: "Undone: the earlier inputs are back." (The page is back at 2.00–3.00 Å; press
  **Compute** to redraw.)
- **Say:** "Coordination dropped to 5.66 because 2.4–2.9 Å cuts both ends of the split Ta–Se
  shell (2.21–3.23 Å). Undo is one click, and if I had edited the page while the card waited, Apply would refuse."

**4b. An interpretation question.** Type exactly:

> Has this run converged, and how good is the fit?

- It runs the **Convergence check** (a step card; expand it for the exact numbers), then answers:
  not converged, the total χ² is still falling beyond its noise over the last tenth of the run;
  Rw 0.0456 for F(Q) and 0.0119 for xPDF. **Outcome: Achieved**, "1 analysis ran and passed its
  checks, and no check contradicts the answer".
- **Say:** "The numbers come from the analyses, not from the model's memory, and every answer ends
  with an Outcome the app computes: if the answer contradicted a check, it would say Not
  verified, whatever the model claimed."

### 5. Validation (0:30)

- **Say:** "On 14 held-out synthetic cases with planted truth, the Copilot with Claude Sonnet 5.5
  scored 0.920, against 0.800 for the default pages read by fixed rules. An exhaustive sweep of
  every engine with no model scores 0.948, so it beats the default pages but not yet the sweep."

### 6. Close (0:15)

- **Say:** "Judge the fit, read the local structure, decide the next run, with an assistant that
  computes, asks before it changes anything, and is scored against known answers."

## Fallbacks

**The model is slow or absent.** Keep talking through the step, and switch to the screenshots
(all from the Demo, © 2026 Tsung-Han Yang):

| Moment | Screenshot |
|---|---|
| The change card | [assets/workbench-copilot-approval.webp](../assets/workbench-copilot-approval.webp) |
| Applied, with Undo | [assets/workbench-copilot-applied.webp](../assets/workbench-copilot-applied.webp) |
| Undone | [assets/workbench-copilot-undo.webp](../assets/workbench-copilot-undo.webp) |
| The convergence answer and its Outcome | [assets/workbench-copilot-answer.webp](../assets/workbench-copilot-answer.webp) |
| Bond Geometry, 80–100° range, bonds and arcs | [assets/workbench-bond-geometry.webp](../assets/workbench-bond-geometry.webp) |

A cold model is the usual cause: gemma4:26b took 7–85 s to reach a change card in the QA, and
60–100 s for a turn that asks a question and runs an analysis. Press **Stop** past about 2
minutes.

**The network is down, or the hosted site cannot reach Ollama.** Use the local build (prep
checklist), which needs no network once built and which Ollama allows by default. If `Test` fails,
the settings drawer says why (server not running, or the origin not allowed).

**"Not verified" on a right coordination answer.** A question such as "What is the mean Ta–Se
bond length and the Ta coordination number?" gets a right answer (2.64 Å, 6 Se per Ta) marked
**Not verified**: the check reads "coordination of 6" as Se's, which has 3 Ta. Rephrasing as "How
many Se per Ta…" did not avoid it in the check of this script (the Outcome's own summary line
still says "coordination is 6"). Avoid that question live. If it comes up, explain it honestly:
"The checks are deliberately strict and frozen with the validation protocol; this is a known false
alarm on a pair's coordination, listed in the QA, and it errs towards *Not verified*, never towards
a false *Verified*."

**Other things a live model may do** (none is an app error): propose the change again after a
refusal (press **Decline** on the second card); run one analysis more than asked; take a minute
to think before a change card.
