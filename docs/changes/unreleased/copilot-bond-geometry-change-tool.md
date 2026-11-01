### Copilot: its first change tool, the Bond Geometry inputs, asked first and undoable

- The Copilot can now set the Bond Geometry inputs, the one thing the Workbench lets a user change
  (the user's decision, 2026-11-01, MATERIA-style). `set_bond_geometry` takes the page's controls —
  `end1`, `apex`, `end2` (the triplet, B central), `r12_min` / `r12_max` (the A–B window),
  `distinct_bc` with `r23_min` / `r23_max` (the B–C window), `bin_width` — each optional, the
  others kept as the page has them. Its schema is strict (an unknown argument is refused), with the
  page's own ranges (windows 0–15 Å, bins 0.1–45°) and the run's elements. It sets the controls
  through the page's own handlers, as the user's typing would, and the page then computes as its
  Compute button does. Its method advice points the model to the coordination-shells skill for the
  windows (the whole first shell, `find_gr_peaks`' `first_shell.window_A`, unless the user asked
  for another), and to `compute_bond_angles` for an answer.
- It asks first. The change waits on an approval card in the dock: one plain line, the page after
  the change (each control's before → after behind "?"), with Apply and Decline. A declined change
  returns `{"declined": true}` to the model. The settings drawer's new *Changes* field switches to
  Auto, which applies changes without the card. Auto is never saved, so every reload starts at Ask
  first. When a change waits while the dock is collapsed, the rail says so.
- Stale-state check: a card shown on one state of the page is never applied to another. If the
  inputs, the result shown, the run or its configuration moved while the card waited, the change is
  refused ("Not applied: …") and the model gets the inputs now on the page.
- It is refused in plain words when no run is open, when the page cannot read the run, or when the
  values contradict each other or would leave the page unable to compute (the refusal names the box
  and the argument that sets it). When the page already shows exactly that, nothing is asked.
- Undo: each applied change is the Copilot's entry in a small change history (`actor: 'agent'`),
  and its card has Undo. Undo puts back the inputs from before the change while the page still holds
  the Copilot's inputs on the same run and configuration. If the page showed a result before, it is
  computed again; if not, the Copilot's result is taken off. A refused Undo says why.
- Every Copilot tool now has its effect — `compute`, `view` or `change` — in one list
  (`src/llm/copilot/effects.js`). The dock reads its tools through it with every effect. The
  validation's live arms read theirs compute-only (`liveArmTools`), and the turn offers a change
  tool only with an approve port, which the harness never passes. A change step is not an analysis
  (the footer, the Outcome, the evidence gate).
- Model-visible: one new tool in the dock's request (its definition, `dock/tools.json`, and the tool
  list in `dock/first-request.txt`). The system prompt, every existing tool description and skill,
  and the arms' requests are unchanged; the harness snapshots pass as they were, and the protocol
  version stays. The dock's worst turn takes 721 more tokens: it fits the live runs' 48k with 3,520
  to spare (was 4,241), and arm B's is unchanged at 7,044.
- Tests: the effects list against every builder's tool and both adapters
  (`toolEffects.test.js`, `validation/__tests__/computeOnly.test.js`); the change-tool contract in
  the turn and the loop (`changeTools.test.js`); the strict schema (`toolCalls.test.js`); the page's
  controller (`BondGeometryControl.test.jsx`) and inputs model (`bondGeometryInputs.test.js`); the
  tool on the page (`copilotChangeTools.test.jsx`: schema, refusals, Ask first, Auto, a decline,
  the stale check, the page's controls set and recomputed, Undo); the approval card through the
  dock (`changeCard.test.jsx`); and the dock beside the page end to end
  (`copilotChangeDock.test.jsx`).
