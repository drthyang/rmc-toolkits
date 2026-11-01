### Copilot: numbers by reference always filled or sent back, and short answers in the dock

- The user saw `σ²θ = {{r5.data.average_structure.angle_variance_deg2.mean}} deg²` raw in the dock
  (2026-11-01), and gemma4:26b's first live runs lost otherwise right answers the same way. The model
  reads each result as `{"ref": "r5", "checks": […], "data": {…}}` and wrote its paths from that root,
  while the app read them below `data`. A path written from the tool message's root
  (`{{r5.data.<path>}}`) now names the value at `<path>`: the user's sentence reads "σ²θ = 111 deg² and
  a trans angle of 159.6°". This is exact, not a guess: no tool's data has a top-level `data` key (a
  test keeps it so).
- An answer whose references still name no value goes back to the model once, before anyone reads
  it, with what each path reaches: the object it names and the values it holds, the keys at the step
  that names nothing, or that no result has that ref. The model writes the answer again with exact
  paths, or without the number. The app never fills a near miss. The repair needs no tool, so it also
  runs on the last request. A reference still unfilled after it is a finding under the answer, and
  the Outcome is not verified.
- The dock never shows a raw `{{…}}`. While the answer streams, references the turn's results fill
  show their numbers, and one still being written is held back. Any reference left over is marked
  *number missing*, and its hover names it.
- Answers in the dock are shorter. With the view tools (the dock only), the prompt now asks for the
  conclusion first in one or two plain sentences, then only the few numbers that carry it (at most
  five, unless the user asks for more), with their units and windows, then the caveats in one
  sentence. It also asks for no field-by-field restating, no full distributions and no tables unless
  the user asks: each analysis's step card holds all its numbers.
- A long finished answer shows its conclusion, and the rest sits behind *Show details*. The rest
  folds when it holds a table, a heading or code, when the answer runs to 700 characters, or when the
  rest holds more than 8 numbers. An answer in the prompt's shape shows whole.
- Model-visible:
  - the dock's system prompt gains `ANSWER_RULES` (six system-prompt variants with the view tools,
    and `dock/first-request.txt`). The dock's worst turn takes 262 more tokens, leaving 3,258 to
    spare in the live runs' 48k (was 3,520);
  - the reference repair's message is new (`harness/reference-repair.txt`). Arms B and C read the same
    words, and only when a reference fails.
  - The arms' prompts, the tools, the skills and the scoring are unchanged. Arm B's worst turn stays
    at 7,044 tokens to spare.
- Protocol 2026-11-01.1: a verdict that held `{{rN.data.…}}` now parses, and a run can make one
  repair request. Runs under 2026-10-24.6 are not comparable. The record keeps the repair in `gates`
  (id `reference`) and its withdrawn answer in `retracted`. `evidence.gates` keeps the evidence
  gate's only, and the report counts the answers sent back beside the references filled and unknown.
- Tests:
  - the user's sentence in a real turn on the Demo, and an object path repaired
    (`answerReferences.test.js`);
  - the store and `explain` (`references.test.js`);
  - the repair in the verification, the turn policy and the loop;
  - a scripted harness run repaired and scored, its record pinned;
  - the dock's display (`referenceDisplay.test.jsx`) and the fold (`answerParts.test.jsx`);
  - the dock-only prompt (`copilotPrompt.test.js`).
