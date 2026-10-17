# AI Copilot (experimental LLM module)

The **AI Copilot** connects the RMCProfile Workbench to an **OpenAI-compatible** LLM: a **local**
server the user runs (Ollama, LM Studio) or a **cloud** provider (OpenAI, Gemini, Anthropic) with an
API key. It replaced the chat-only AI Assistant, which answered from one fixed
summary of the run; the Copilot also **runs the Workbench analyses itself**, reads what came back,
and judges whether that answers the question.

It is **docked beside the pages**, not over them, and it **points at the page** while it explains:

- **The dock.** On a window at least 1080 px wide (and taller than 540 px: laptops, monitors,
  tablets in landscape) the Copilot is a side
  panel in its own grid column: the pages shrink and reflow beside it (their width tiers are
  container queries), so it never covers a card. Drag its left edge (or focus the edge and use
  ←/→, Home/End; double-click or Enter resets) to resize it; the pages always keep at least
  760 px. On tablets in portrait, phones and short windows it is a **bottom sheet** in its own grid row,
  resized from its top edge; when the agent brings a card on screen, the sheet comes down until the
  card fits whole above it (to at most 50 % of the window while a card is in focus, never below
  28 %). Collapsed (the chevron, or Esc inside it) it is a slim **rail** at the
  right edge, or a bar at the bottom, with a red dot when an answer arrived while it was
  collapsed. The open state, the width and the sheet height are remembered
  (`localStorage` `rmc-llm-dock-v1`).
- **Pointing at cards.** The agent's `highlight_cards` tool gives one or more cards a breathing
  blue edge (2.4 s, the same language as MATERIA and NEBULA3D; a steady edge under
  `prefers-reduced-motion`), opens the card's page when needed and scrolls it into view; on a
  plot it can also mark x positions (the g(r) peaks it names), and the plot zooms to them. The
  card of an analysis that is running breathes softly. A highlight ends when the agent points at
  something else, a new question starts, the Copilot collapses, or the user clicks anywhere on
  the page.
- **Card chips.** The answer cites cards as `[label](card:<id>)` (or `[[id]]`); each renders as a
  chip that brings the card back on screen.
- **3D views.** The host's `orient_3d_view` turns a 3D view to look down a crystal direction
  `[u v w]` (or a, b, c, PC1–PC3), centres or rings sites in blue, selects a site as a click would,
  or resets it; `current_view` says what is on screen. So the agent can say "looking down [111],
  the displacements…" and show it.
- A step's **Show in \<page\>** button switches the workbench page and loads the same analysis
  there. Bond Geometry gets the triplet and windows and computes; Displacement Directions gets the
  site and resolution. The page shows its own result beside the conversation.
- **It asks instead of guessing** (`ask_user`). When a question leaves open a choice that changes
  the science — which coordination shell or r-window, which site or pair, which tolerance, which
  run — the agent asks one short question with 2–4 options and its default, then uses exactly the
  option the user picks. Where the run data show that such a choice changes the result, the
  Workbench asks itself, before the model guesses. See *Asking the user* below.
- **It follows skills** (`read_skill`): written methods for the questions it gets, shared with
  Claude Code (`.claude/skills/`), drafts awaiting the user's interview. See *Skills* below.
- **Its answers are checked.** Where an answer contradicts what the analyses or the run context
  established, a **Check** line says so under it and the Outcome chip drops to **Not verified**,
  whatever the model said of its own answer. See *Checks on the answer* below.

How one question runs:

```
        question
           │
   ┌───────▼────────┐   tools offered    ┌──────────────┐
   │  agent loop    │ ─────────────────▶ │    model     │
   │ (agentLoop.js) │ ◀───────────────── │ (client.js)  │
   └───────┬────────┘  tool calls / text └──────────────┘
           │ validate → run → fit to budget
   ┌───────▼─────────────────────────────┐
   │ tools: run data (runTools.js) and   │   result + deterministic checks
   │ engines from the host (tools prop)  │ ──▶ `tool` message ──▶ next round
   └─────────────────────────────────────┘
           │ … until the model answers, ending with
           ▼   Outcome: achieved | partial | not achieved — <why>
   transcript: text, step cards (ToolStep), Outcome chip
```

The loop has four stages: **plan → act → observe → assess**. The assessment has two halves:

- every tool result carries **deterministic checks** computed from the full engine output. For
  example, "0 angles in these windows — widen them", "no significant peak (1.2σ < 3σ)", or "the
  first g(r) peak sits on the closest-approach constraint".
- the model ends its answer with a verdict line, which the chat shows as a chip.

A model that cannot call tools still works: the Copilot falls back to plain answers from the run
summary. The old Assistant behaved the same way.

The **convergence watchdog** badge on the Dashboard's χ² chart is part of the same module and is
unchanged.

## Privacy

Local providers keep the promise that run data never leaves the device. The browser talks
**directly** to `http://localhost:*`, and the only party that sees the run-derived summaries is the
server the user runs.

**Cloud** providers receive the run context and the **results of the analyses the Copilot runs**.
Both are summaries of the same kind: site statistics, angle modes, significances, never the
configuration file. Cloud use is opt-in, and the settings drawer warns that data leaves the device.
The *Let the Copilot run analyses* switch turns the tools off. Nothing is sent until the user asks a
question. API keys are stored in `localStorage` only.

**The analyses only compute; the view tools only change what is on screen.** No tool moves a
slider, changes an analysis input or the Detected-SG tolerance, or writes or downloads a file.
The view tools (`highlight_cards`, `orient_3d_view`, `current_view`) open a page, scroll, give a
card its edge, mark a plot, turn a 3D camera, ring sites, and select a site as the user's click
would — and nothing else. A step card's **Show in \<page\>** button, which loads an analysis's
inputs on its page, stays the *user's* action.

## How each seam works (reading order)

1. **`context/runContext.js` + `context/pairCorrelations.js`** build the run context.
   `buildRunContext()` distills the dashboard's parsed state into a compact JSON object:
   - cell and angles;
   - the symmetry block: space group, ladder, and Wyckoff orbits ranked by rms displacement;
   - the symmetry noise floor: below it, the ladder's group changes are sampling noise of the site
     means (rms displacement / √copies);
   - pair correlations, with each first shell's mean distance (from r²g, not the g peak) and
     coordination number;
   - move counters;
   - run settings from the `.dat` file, saying when one measurement is fitted in both spaces;
   - per-dataset Rwp, a Fourier-transform file marked as the same measurement;
   - the Faber–Ziman weight of every partial in each fitted radiation (`scatteringWeights.js`);
   - the ln(χ²) history, downsampled to ≤ 48 points with statistics computed on the full series, the
     watchdog's verdict and the change over the last tenth of the run.

   `contextToJson()` serializes it compactly (arrays and small objects on one line), enforces a
   ~6,000-character budget and records each trim with an `*_omitted` count. The context travels with
   every question, so simple questions need no tool.
2. **`prompts/system.js` + `prompts/copilot.js`**
   - `SYSTEM_PROMPT` holds the domain rules: what Rwp means, that the history is ln(χ²) of one fit
     term, and "quote the numbers, don't guess".
   - `COPILOT_PROMPT` adds how to work with tools: say why before each call, read the `checks`,
     change the arguments on a `fail`, never repeat an identical call. It also asks for the verdict
     line. `copilotSystemPrompt()` adds `ASK_RULES` when `ask_user` is offered.
   - `buildCopilotMessages()` puts the context, inside a fenced JSON block, ahead of the
     conversation. Earlier turns travel as question + final answer only (the newest 8 messages: 4 exchanges).
   - `prompts/templates.js` holds the watchdog's rigid `STATUS: …` prompt and its parser.
3. **`provider/client.js`** is a hand-rolled OpenAI-compatible client (`/models`,
   `/chat/completions`) with an SSE stream parser over `fetch` + `ReadableStream`.
   - `streamChat()` yields `{ content }`, `{ reasoning }`, `{ toolCall }` fragments and
     `{ finishReason }`.
   - `tools` / `tool_choice` join the request body only when given, so watchdog requests are
     unchanged.
   - Each fragment's provider data (`extra_content`, where Gemini puts its `thought_signature`)
     is kept, because it must be echoed back.
   - An HTTP failure carries its `status`.
   - **`provider/contextGuard.js` + `provider/ollama.js`**: every Copilot request goes through
     `streamCopilotChat()`, which never sends a request the model would cut (*The context window*
     below). Ollama's requests use its native `/api/chat` (NDJSON) with a fitting `num_ctx`.
   - `checkConnection()` turns failures into actionable hints: a bare `TypeError` means "server
     down **or** CORS blocked", and 401 / 429 name the key and quota problems.
4. **`copilot/toolCalls.js`** turns streamed fragments into calls.
   - `accumulateToolCall()` keys fragments on `index`, then `id`, then arrival order, and
     tolerates servers that send whole calls or reuse an index.
   - `parseToolArguments()` accepts `''` and fenced JSON.
   - `validateToolArguments()` checks the JSON-Schema subset the tools use: types, `enum`,
     `minimum` / `maximum`, `required`. Unknown keys are dropped with a note. Numeric strings
     are accepted for number fields.
5. **`copilot/tools.js`** holds the tool contract (below) and `executeToolCall()`: parse, validate,
   run with a timeout, and give the tool its own `AbortSignal`, which is aborted on Stop or timeout.
   A thrown error becomes a failed step whose message the model reads. `fitToBudget()` halves the
   largest array, recording `*_omitted`, until a result fits. It never trims silently.
6. **`copilot/agentLoop.js`**: `runCopilotTurn()` is a pure async generator of UI events. The model
   call is injected, so the tests script the model's replies. Its rules:
   - **Budgets.** ≤ 6 rounds of tool calls; ≤ 4 calls per round (extra calls are answered
     "skipped"); 120 s per call; no new round after 10 min. Each result is fitted to 3,000
     characters, and a turn's tool messages to 12,000 in total. Older results are compacted to
     their checks first, so a small local context window keeps the run context and the rules.
   - **Duplicate calls.** An identical call (same validated arguments) returns the cached result
     with a note.
   - **Final round.** When the budget is used, a last request goes out with `tool_choice: 'none'`.
     If it still answers with tool calls only (Ollama ignores `tool_choice`), a fixed
     "budget used" answer with `Outcome: not achieved` closes the turn.
   - **No tool support.** A 4xx saying the model "does not support tools" retries the turn once
     without tools, under the tool-free prompt. Other 400s, such as a schema error, are shown.
   - **Message rules.** Every listed call gets exactly one `tool` message. Arguments are echoed as
     valid JSON (`'{}'` when unparseable). Whitespace-only assistant content is sent as `null`
     (Anthropic rejects empty text blocks).
7. **`copilot/runTools.js`** holds the tools over data the module already receives:
   - `read_run_data(section)`: one context section as built, before the budget trim;
   - `assess_convergence`: the watchdog heuristics plus tail slopes;
   - `find_gr_peaks(pair)`: partial-g(r) peaks and the first coordination shell — its window,
     mean distance, coordination numbers and their sensitivity to the cut, whether it is separated
     from the next shell — and whether the pair's minimum-distance constraint cuts the distribution;
   - `compare_fit(dataset)`: where the model misses the data, as Rw and the share of the misfit
     per range of r or Q, and the largest residual features.
8. **`copilot/askUser.js`**: the `ask_user` tool (*Asking the user* below).
   **`copilot/outcome.js`**: `parseOutcome()` splits the verdict line from the answer. It is
   lenient about Markdown and dashes, strict about the three statuses, and takes the last line.
   **`copilot/turnState.js`** folds loop events into a turn: ordered text segments (with
   reasoning) and tool steps.
9. **`watchdog/`**: slope heuristics over the ln(χ²) history are the **source of truth** for the
   badge. The LLM only writes the tooltip note, and only when the status flips or enough new data
   arrived.
10. **`useAssistant.js`, `useCopilot.js` and `components/`**
    - `useAssistant()` holds the settings, the connection probe (with auto-connect) and the run
      context.
    - `useCopilot()` holds the conversation: Stop, which models cannot call tools, and a new
      conversation when `runKey` changes. It also builds the memory (`conversationMemory`): a
      question without an answer is dropped.
    - `CopilotChat` is the dock (`components/dockLayout.js` holds its pure size rules): the rail,
      then the panel, and the resize handle. The panel mounts on first open, so nothing is fetched
      before that, and stays mounted while collapsed, so a running answer continues. It composes
      the connection bar, the settings drawer (`ConnectionSettings`, with the analyses switch) and
      `CopilotView`. Focus goes to the composer on open (once connected) and back to the rail on
      collapse. It owns the card-focus controller.
    - `CopilotView` renders the transcript (`ToolStep` cards, `MarkdownMessage` with its card chips,
      the Outcome chip), the composer, and a polite live region that announces step and verdict
      changes.
11. **`focus/focus.js` + `focus/cardRefs.js` + `copilot/viewTools.js`** — pointing at the page.
    - `createFocusController()` sets `data-copilot-focus` (primary | secondary | running |
      preview) on the host's `[data-copilot-card]` elements; `Copilot.css` draws the breathing edge.
      `show(ids, { marks })` opens the main card's page through `onOpenPage(key, null)` when it is
      not on screen (a card under `[aria-hidden="true"]` is off screen), waits up to 3 s for it to
      mount, lowers a tall sheet, and scrolls it into view (not at all when it is already in view;
      smooth unless reduced motion). A MutationObserver re-applies the attribute to a card the host
      re-renders. `clear()` runs the tools' `onClear` undo callbacks (a 3D site ring).
    - `cardRefs.js` reads a model's spelling of a card (`card:` scheme, brackets, case, the title, a
      bare plot file name) and rewrites `[[id]]` as a `card:` link; `MarkdownMessage` lets the
      `card:` scheme through react-markdown's URL filter and the sanitizer and renders known cards as
      chip buttons (unknown ones as their label).
    - `viewTools.js` builds `highlight_cards` from the host's `cards`: the description lists every
      card for the model, and run() resolves the ids, keeps marks for plot cards (`marks: true`),
      and reports what was shown and what was not on screen.

## Asking the user (`ask_user`)

The spot check's worst answer was a guess: asked for "Se–Ta–Se angles against an ideal
octahedron", a model picked the first-peak window (2.33–2.71 Å) and cut the TaSe₆ shell in two,
then read the cut as "incomplete octahedra". Which window to use is the user's call when the data
offer two, so the agent asks:

- **The tool** (`copilot/askUser.js`, offered with every run's tools):
  `ask_user({ question, options: [{ label, description? }] (2–4), allow_other = true, default? })`.
  `normalizeAskRequest()` trims the question, reads plain strings as options (and splits
  `"label — description"`), drops repeated labels, and matches `default` to a label or its number;
  it refuses an empty question or fewer than two distinct options with a message the model can act
  on.
- **The card** (`components/QuestionCard.jsx`) shows in the conversation, in the docked panel and
  in the bottom sheet: the question, the options as buttons with their one-line meaning (the
  default marked), an *Other…* field and *Skip*. Keyboard: Tab, ↑/↓, 1–4 to pick, Enter in the
  field. It takes focus when it appears unless the user is typing in a field, and the log never
  scrolls its question out of view. A collapsed dock marks the rail ("a question for you").
  Answered, it keeps the question and the reply.
- **The loop** waits. `useCopilot` passes `askUser(request, { signal })` in the tools' context; the
  card settles it, and the answer returns as the tool's result (`answerResult()`: the chosen label
  and its meaning with "use exactly this choice", a typed answer verbatim, or a skip as "use your
  default and say you assumed it"). Stop, a new run or an unmount rejects it with an AbortError:
  the turn closes as *Stopped*, the card reads *Not answered*, and the next question works.
- **The contract.** `interactive: true` marks a tool that waits for the user: never answered from
  the duplicate-call cache, no time limit, its wait not counted toward the turn's 10 min, at most
  `maxQuestions` (2) questions per request and one per step (`agentLoop.js`). It is not a view
  tool and changes no data; the turn footer does not count it as an analysis.
- **The rules** (`ASK_RULES` in `prompts/copilot.js`, added when the tool is offered): ask only
  when the answer changes which analysis runs or its parameters and neither the context nor an
  analysis can settle it; settle what the data can first (e.g. `find_gr_peaks`, then offer its
  shells as the options); never on a clear question or for permission to run an analysis; one
  question at a time, at most two per request; offer a default; then use the chosen option and
  say which.

### With a real model (2026-10-11)

gemma4:26b (Ollama, through the context guard at `num_ctx` 24,576) on the Demo, docked at
1440 × 900, each question in a new conversation, a question card answered with its default:

| Question | Asked? | What it did |
|---|---|---|
| Has the run converged? · first Ta–Se peak, constraint-limited? · most anharmonic site and largest U_iso · Ga preferred direction? · space group at 0.2 and 0.02 Å | no | answered from the context or one analysis, correctly (improving; 2.52 Å; ref 23; no direction at −0.89σ; the P-42₁m rung below the noise floor) · 6–37 s |
| Se–Ta–Se angles against an ideal octahedron · "What are the Se–Ta–Se bond angles?" | no | used the whole first shell, 2.2–3.23 Å (CN 6), from the context or `find_gr_peaks` · 18–28 s |
| "Compute the bond angles around Ta" · coordination number of Ta · Ta–Ta–Ta angles · Ta–Se bond lengths | no | chose the main polyhedron (Se–Ta–Se) and said so, or covered every option (Se, Ta and Ga neighbours; both Ta–Se peaks and the shell mean) · 9–28 s |

It never called `ask_user` in 16 questions, even after the rules named "the bond angles around
Ta" as a case to ask about: where the context settles the choice (the whole shell) that is right,
and elsewhere it stated its assumption. The card and the loop are exercised with a scripted model
(*Tests*). One answer was wrong in the same way twice: it read the raw angle counts (modes 74.5,
88.5, 104.5°; mean 103 ± 29.6°) as "significantly distorted" octahedra, missing the trans angles
the sin θ factor depletes near 180°.

### The Workbench asks itself (2026-10-11)

gemma4 never asked; so where the computed signals show a real ambiguity that changes the result,
the Workbench shows the question card itself, with the options it computed and a default
(`copilot/clarify.js`, `buildClarifier`). The card is an `ask_user` step marked `app: true`
("The Workbench asks before the analysis"); it counts toward the two questions per request, and
the model's own `ask_user` still works beside it.

- **Before the first request** (`beforeTurn`), on the user's question (`readQuestion`: the
  intent, and the run's elements, pairs and triplets it names):
  - *a bond-angle or coordination request naming only an element:* which neighbours — one
    option per bonded partner (a first shell starting within 1.4× of the nearest one), nearest
    first and the default, each with its window and CN; "Every neighbour type" for coordination.
    The Demo's "bond angles around Ta": Se–Ta–Se (2.21–3.23 Å, CN 6), Ta–Ta–Ta (2.65–3.33 Å, CN
    2.8); not Ga at 3.8 Å;
  - *a named pair or triplet whose first shell is split into peaks, or runs into the next:*
    which window — the whole first shell (default), the first peak only, or a cut 0.2 Å either
    side of the minimum — only when the coordination number changes by ≥ max(0.5, 15 %) between
    them. The Demo's "Se–Ta–Se bond angles": whole shell 2.21–3.23 Å (CN 6) or first peak
    2.33–2.71 Å (CN 3.9). The Ta–Ta shell runs into the next (6 % of the peak at the minimum) but
    a cut 0.2 Å either side changes its CN by < 0.5: not asked;
  - *displacement directions of an element with several sites whose site symmetry is not
    cubic:* one site of each orbit (default; its sites are equivalent by symmetry) or all pooled
    — pooling mixes the sites' local frames and can cancel a per-site preference, the review's
    pooled-Ta error. Ga (site symmetry -43m) is not asked about;
  - *the space group with no tolerance named:* which tolerance, when more than one group holds
    clearly above the sampling noise floor (to ≥ 1.25× it). On the Demo only F-43m does: not asked.

  The answer joins the user's question for the model: `(The Workbench asked: "…" — the user
  chose Se–Ta–Se: Se around Ta in the whole first Ta–Se shell, 2.21–3.23 Å (CN 6). Use exactly
  this.)`. The conversation still shows the question as typed.
- **Before an analysis runs** (`beforeCall`, for tools that declare `shellWindows(args)`:
  `compute_bond_angles`, `analyze_clusters`): a window that overlaps a first shell but holds less
  of it by the same margin is put to the user — the whole first shell (default) or the window as
  asked — and the analysis runs in the one chosen, its tool message saying so. The spot check's
  2.3–2.7 Å Se–Ta–Se window: "it holds 3.8 Se per Ta of the shell's 6". Once per question: the
  choice is remembered by (centre, partner, window), and a later call with that window runs in the
  chosen one unasked (a model that repeats its cut window every round is asked once); a call
  identical to one already answered from the cache asks nothing. A typed "Other" window goes to the
  window the question was about (r12 or r23).
- **Not asked:** a question that settles the choice itself — an ideal polyhedron it names that
  exactly one option's CN matches ("…compare them with an ideal octahedron": only the whole Ta–Se
  shell holds 6), the shell it names ("only the first Ta–Se peak" keeps a partial window as the
  user's, confirmed; "the whole first shell", "the entire shell", "both peaks" run a cut window in
  the whole shell), an element pair or triplet with a clear shell, a site number, a window or a
  tolerance in Å, a window the user chose earlier in the question. The six spot-check questions
  ask nothing.

## Checks on the answer

The reviews found answers that contradicted what the tools had just computed. The tools now
return `facts` beside their data (not sent to the model): `assess_convergence` its status,
`detect_symmetry` its group with the point group's facts, `analyze_displacement_directions`
whether a direction is significant (peak or map ≥ 3σ); the run context adds the watchdog's status
and the tolerance ladder's groups (`contextAuditFacts`). When the model answers,
`auditAnswer()` (`copilot/answerAudit.js`) reads the answer, sentence by sentence and clause by
clause, for a clear positive claim they contradict:

| Check | Fires on |
|---|---|
| convergence | "converged", "reached convergence", "a plateau" while the history is still improving or diverging |
| chirality | a group called chiral whose point group is not a Sohncke group (it has mirrors, glides or rotoinversions) |
| centrosymmetry | a group called centrosymmetric whose point group has no inversion centre |
| inversion lost | an inversion centre "lost" or "broken" when no group it could be lost from has one: neither the groups the answer names nor the parent (the context's average group, or the ladder's top) |
| polar | a group called polar whose point group is not polar |
| direction | a "preferred direction" when no directions analysis of the question found a peak at ≥ 3σ (the whole-map test alone is anisotropy, which thermal ellipsoids show too: review L3) |
| superstructure | a superstructure, cell doubling or zone-boundary order when no `analyze_modulations` of the question found an ordered point |

A negated claim ("has not converged", "achiral", "no preferred direction"), a question, a
hypothetical or conditional mention ("to check whether it has converged", "until χ² plateaus", "I
tested whether site 23 has a preferred direction", "a converged run would…"), or a sentence that
already says the run still moves raises nothing. A negation counts only in the stretch of the
clause just before the claim, after its last comma or "and"/"so": "With no remaining drift in χ²,
the run has converged" is a claim. A chirality claim in a clause that names a group the checks do
not know ("its subgroup P2₁2₁2 is chiral") is not pinned on the sentence's other groups. Symbols
are matched however written
(`P-42₁m`, `P−4 2_1 m`, `$P\bar{4}2_1m$`). Each finding shows as a **Check** line under the
answer, and the Outcome chip is **lowered** (`lowerOutcome`): an answer the model called
*achieved* or *partial* (or gave no verdict) becomes **Not verified**, with "the model said
'achieved'" in its note — the model's own verdict never overrides a computed contradiction;
*not achieved* stays. gemma4:26b's correct answers of 2026-10-11 raise no check.

**Facts the model reads.** `detect_symmetry` and the context's ladder state each point group's
facts (`context/pointGroups.js`: crystal system, inversion centre, Sohncke class, polarity), so the
model need not recall them (the spot check called P-42₁m chiral). `compute_bond_angles` returns the
10° histogram sin θ-corrected as well (1 for random directions: raw counts deplete 180°; the
engine's bins are shared out by overlap, so 15° or 30° bins leave no comb), and for a whole shell
of about six or four neighbours, in an A–B–A triplet with one window (a mixed or split triplet's
angles are a subset of the polyhedron's), the split against an ideal octahedron (80 % cis at
90°, 20 % trans at 180°) or tetrahedron: the Demo's TaSe₆ shell is 80 % cis at a mean 89° and 20 %
trans at 158.8°, which gemma4 had read from the raw modes as "highly distorted, 75–105°". A cut
window gets no polyhedron reading (its CN 3.9 is not a tetrahedron).

## Skills

**Skills** are Claude Code skills at the repository root, `.claude/skills/<name>/SKILL.md` (front
matter `name` and `description`), so the in-app agent and Claude Code read the same text. Four are
drafted, each opening with **"Draft: awaiting the user's interview"**, with its order of steps,
judgment calls marked "(to confirm with the user)", an acceptance checklist and known failure modes
from the spot check and the review:

| Skill | For | Its firm rules in code |
|---|---|---|
| `rmc-fit-assessment` | convergence, where the fit misses (`compare_fit`), constraint-limited peaks, sampling, overfitting, the report | the convergence check on the answer (*Checks on the answer*) |
| `coordination-shells` | bond windows from the partial g(r), coordination against the ideal polyhedron, reading angle distributions | `compute_bond_angles` **refuses** a window that cuts the shell (fewer than 60 % of the central atoms share the most common count, or most have none; the 60 % is an interview item) and returns the counts that show the cut, not the angles, unless `window_confirmed` — which the Workbench's window question sets when the user keeps their window; the polyhedron reading |
| `local-symmetry` | the tolerance ladder, what breaks between rungs, the noise floor, point-group facts | the point-group facts and the chirality / centrosymmetry / polar / inversion checks |
| `displacement-anisotropy` | ellipsoids, kurtosis, preferred directions and their σ, effect sizes | the preferred-direction check; the Workbench's one-site-or-pooled question |

- **Progressive disclosure** (`copilot/skills.js`). The system prompt lists only each skill's name
  and description (`skillIndex`, 2,150 characters for the four); the agent reads a skill
  (4,700–5,500 characters each) with `read_skill(name)` and a reference file beside it with
  `read_skill(name, reference)`. `read_skill` has its own result budget (`resultBudget` 9,000) and
  is never compacted later in the question (`keep: true` in the tool contract), so the agent can
  follow the skill throughout; its step is not counted as an analysis.
- **The host bundles them** (`src/copilotSkills.js`): the SKILL.md files with the app, the reference
  files as lazy chunks loaded when read; it passes `{ files, references, load }` as the `skills`
  prop, and the module parses them (`parseSkillFiles`: a misnamed skill or an unbundled reference
  is an error, shown in the panel). Vite reads them from the repository root
  (`server.fs.allow`); the Docker build copies them in.
- **Replays** (`src/__tests__/copilotSkillReplays.test.js`): per skill, a scripted model on the
  Demo run with the real engines; the skill's right behaviour passes, and the old wrong one
  (a cut window, "converged — a plateau", "P-42₁m is chiral and centrosymmetry is lost", a
  preferred Ga direction) is refused or flagged.
- **Interviews** turn the drafts into the user's method: `docs/SKILL_INTERVIEWS.md`.

## The context window: never cut silently

A local server runs a model with a fixed context window, and a prompt longer than the window is not
refused: it is shortened to fit, without an error, and the model answers without part of what it
was given. The Copilot's prompt is long. On the Demo the first request is ~28,400 characters (the
rules 8,500, the run context 5,000, fourteen tool schemas 14,800; measured before `read_skill` and
the skills index were added, which bring 15 tools and ~2,150 characters more), which gemma4's
tokenizer reads as ~8,600 tokens (the Ollama log recorded 8,025–8,421 for the same request before
`ask_user`);
each analysis adds up to 3,000 characters more.

- **Ollama's OpenAI-compatible endpoint cannot carry a context length.** It runs the model at the
  server's default: `OLLAMA_CONTEXT_LENGTH`, else one Ollama picks from the GPU memory it finds.
  On the Apple M5 Pro this was written on, 0.40.2 picks 262,144 (51.8 GiB of GPU memory) and the
  models' own maxima then apply (gemma4:26b 262,144, qwen3:32b 40,960), so nothing was cut there:
  Ollama logs `truncating input prompt` when it cuts, and the server logs since 2026-08-26 hold no
  such line from the Copilot. But the same Mac's log of 2026-09-27 shows Ollama starting without
  the GPU (`total_vram="0 B"`) and picking 4,096, half the Copilot's first request; older releases
  defaulted to 2,048–4,096 everywhere.
- **What a cut does (verified 2026-10-11, gemma4:26b, the Demo's first request: 28,382
  characters, 14 tools).** The OpenAI endpoint reports 8,437 prompt tokens at the 262,144 default.
  At `num_ctx` 4,096 the same request still returns 200 with no error, but Ollama first drops whole
  messages (the run context among them) and then cuts what is left to half the window, keeping
  its first 5 tokens (`truncating input prompt limit=2051 prompt=5694 keep=5`): the model reads
  2,051 tokens and, without its rules and the run, calls `current_view` instead of `find_gr_peaks`.
  At 10,240, 16,384 and 24,576 (the guard's choice) it reads all 8,437 and calls `find_gr_peaks`.
  A request whose `num_ctx` equals the context the model is loaded with reuses the runner and its
  prompt cache (no reload; 1 token evaluated).
- **So the Copilot sizes the context itself.** For an Ollama server (`GET /api/version` answers a
  version) it sends its requests to the native `/api/chat`, with `options.num_ctx` = the context
  the model is loaded with when that holds the request and 8,192 tokens of reply
  (`REPLY_RESERVE_TOKENS`), else the next multiple of 8,192 that does (at least 16,384, kept for
  the rest of the conversation, at most the model's maximum from `/api/show`), and
  `options.num_predict` = that context minus the prompt, so a long reply stops at the limit
  instead of pushing the prompt's start out of the window. The size comes from a deliberately
  high estimate (`estimatePromptTokens`: 2.5 characters per token plus 8 per message; ~32 % above
  gemma4's real count on the Demo). Tool results and call arguments, JSON dense with numbers,
  count at 1.4 characters per token (`JSON_CHARS_PER_TOKEN`; gemma4 reads them at 1.64), so the
  estimate errs high there too and `num_predict` never reaches into the prompt.
- **It refuses rather than cuts.** When the request cannot fit even the model's maximum, the turn
  stops before anything is sent, with the sizes and what to do (a new conversation, analyses off,
  a model with a longer context). LM Studio fixes the context when a model is loaded: when the
  loaded model's window (`GET /api/v0/models`) cannot hold the prompt and 4,096 tokens of reply,
  the turn stops with the context length to load it with. A model not loaded yet is loaded by the
  request itself (JIT), with LM Studio's default context: its window is read at the reply's first
  fragment, and one too small stops the reply before any of it is shown. Cloud providers and other
  servers are sent the request unchanged (only `localhost` or port 11434 is ever probed); an API key
  set for an Ollama behind a proxy is sent on the native path too.
- **A backstop.** Ollama's `prompt_eval_count` is what the model read: the whole prompt, cached
  or not, and after a cut what was left. Below 0.4 of the estimate (`READ_SHARE_LIMIT`; an uncut
  request reads 0.74 of it, the cut one above 0.18) the turn stops with an error that says the
  reply did not see all of its instructions.
- **A cut reply says so.** A reply that stops at the length limit (`done_reason` /
  `finish_reason` `length`) ends with a warning in the turn's footer instead of passing as a
  whole answer.

## Robust on local models

A local model can stall, loop or reply with nothing: qwen3.8-27b, a dense reasoning model, took 1–5
minutes a question in the 2026-10-10 spot check and once over 15. `copilot/robustness.js` holds the
pieces the agent loop uses so a turn always ends with an answer or a reason:

| Failure | What the loop does | Limit |
| --- | --- | --- |
| No token for minutes (server busy, stopped) | Cancels the request; the error says the server may be busy | `stallMs` 180 s |
| A reply that never ends | Cancels it; the error suggests a narrower question or a faster model | `requestMaxMs` 600 s |
| Reasoning (or text) repeating a passage | Stops it when caught (`repeatsItself`: its last 200 characters 4 times in the last 8,000), asks again once with `antiRepeat` (temperature ≥ 0.6) | one reminder |
| Reasoning only, no answer, no call | One reminder, then `done` with `stoppedBy: 'empty'` and a `reason` the footer shows | one reminder |
| A call written as JSON text | Run as a call (`textToolCalls`: the whole reply must be the call, naming a known tool); the text is retracted | |
| Calls that repeat earlier ones or fail | Two such rounds in a row end the analyses: the next request asks for the answer (`stoppedBy: 'no_progress'`) | `maxNoProgressRounds` 2 |
| A long turn | Earlier results are compacted to their checks and top-level numbers (`keyNumbers`), once | `turnCharBudget` |
| A long conversation | Earlier answers above 2,000 characters are re-sent as their sentences with numbers and their Outcome line (`compactAnswer`) | `MEMORY_ANSWER_CHARS` |
| Out of memory, context overflow, no model loaded, connection lost | The error gets a line on what to do (`explainModelError`) | |

Notes to the model (a reminder, the last-round note) join the last user message when that is the
last message: Gemma's chat template refuses two user messages in a row.

## The tool contract (the `tools` prop)

The host passes engine tools as plain objects; the module validates and runs them:

```js
{
  name: 'compute_bond_angles',           // /^[a-z][a-z0-9_]{0,63}$/, unique
  description: '…',                      // what it computes and when to use it, for the model
  parameters: { type: 'object', properties: {…}, required: [...] },
  label: (args) => 'Bond angles Se–Ta–Se',   // step-card title (optional)
  page: 'geometry',                      // host page key for the step's button (optional)
  pageParams: (args) => ({ … }),         // what that page loads to show the same analysis (optional)
  cards: (args) => ['bond-angles'],      // host cards this call reads: they breathe while it runs (optional)
  view: true,                            // a view tool: changes only the screen; never cached (optional)
  interactive: true,                     // waits for the user (ask_user): never cached, no time limit (optional)
  shellWindows: (args) => [{ centre, partner, window: [lo, hi], apply: ([lo, hi], { confirmed }) => args }],
                                         // its bond windows, checked against the shells before it runs (optional)
  resultBudget: 9000,                    // its own result budget in characters (read_skill) (optional)
  keep: true,                            // its result is never compacted later in the question (optional)
  run: async (args, { signal, focus, askUser }) => ({ data, checks, facts }),
}
```

- `data` is a compact JSON summary built for the model, never a raw engine payload. Put the most
  important rows first: budget trimming halves arrays from the end.
- `checks` is a list of `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`, computed from the
  full result. `fail` means "this does not answer as asked".
- `signal` aborts on Stop or timeout. Stop the work if you can.
- `focus` is the dock's card-focus controller, for view tools: `await focus.show(ids, { marks })`,
  `focus.onClear(undo)` (runs when the highlight ends), `focus.state()`.
- `askUser(request, { signal })` shows a question card and resolves with the user's reply
  (`{ answer, source: 'option' | 'other' | 'skipped', index }`), for interactive tools.
- `facts` (optional, not sent to the model) are findings the final answer must not contradict
  (*Checks on the answer*): `{ convergence }`, `{ groups: [{ symbol, point_group, centrosymmetric,
  sohncke, polar }] }`, `{ direction: { subject, significant, anisotropic, peak_sigma, map_sigma } }`.
- `shellWindows(args)` names the bond windows a call will use, so the Workbench can ask before a
  window that cuts a coordination shell runs; `apply(window, { confirmed })` returns the
  arguments for the window the user chose (`confirmed`: they kept their own).
- Parameters may be arrays (`items`, `minItems`, `maxItems`) and objects (their own `properties` /
  `required`); a list sent as JSON text or a bare comma list is accepted.

`pages` (`{ key: label }`) and `onOpenPage(key, params)` let a step offer **Show in \<page\>**.

- `params` is the tool's `pageParams(args)`, computed from the validated arguments. Without
  `pageParams`, the button reads **Open \<page\>** and `params` is null.
- `runKey` names the open run; when it changes, the conversation restarts.
- `onOpenChange(open)` reports the dock opening and collapsing.
- `cards` lists what the agent may point at: `{ id, page, title, about, kind?, file?, marks?, view3d?, sites? }`.
  Each card element carries `data-copilot-card="<id>"`. `kind` / `file` tie a card to the data the
  run tools read (a plot kind, `r_value`, `symmetry`, `model`), so `assess_convergence`,
  `find_gr_peaks`, `compare_fit` and `read_run_data` breathe the right card while they run.
- `onMarks({ [cardId]: [{ x, label }] })` hands the host the marks to draw inside a plot card.
- `onLayoutChange({ mode, open, width, sheet })` reports where the dock sits; `className` and
  `storageKey` are the host's.

In this repository the host tools live in `src/copilotTools.js`:

| Tool | What it runs |
|---|---|
| `list_sites` | the per-site PCA table, through the worker or `/api/pca/sites` |
| `analyze_displacement_directions` | the orientation engine with calibrated significances and `[u v w]` (Show in Displacement Directions selects the site and resolution) |
| `compute_bond_angles` | the triplets engine (Show in Bond Geometry fills the triplet and windows, then computes); the sin θ-corrected histogram, and a whole shell of ~6 or ~4 neighbours read against an ideal octahedron or tetrahedron |
| `analyze_clusters` | the cluster engine (`workers/clusters.js`, `/api/clusters`): the shape of an element's finite clusters (M₄ tetrahedra, dimers, trimers) against a shuffle null of independent atoms, and the average cluster's own distortion |
| `analyze_modulations` | the modulation engine (`workers/modulations.js`, `/api/modulations`): each element's displacement power at the zone-boundary points and its strongest modes against the configuration's diffuse background — order that enlarges the cell, which the folded average cannot show |
| `detect_symmetry` | the symmetry finder at any tolerance, without changing the card (browser-loaded runs only), with the point group's facts (inversion, Sohncke, polar) |
| `orient_3d_view` (view) | turns a 3D view (`src/copilotViewTools.js`): look down `[u v w]` / a, b, c / PC1–3 (c up, or b looking down c), centre or ring sites, select a site, reset |
| `current_view` (view) | the page on screen, its cards, the highlight, and each 3D view's direction and site |

The 3D views join a registry while mounted (`src/viewer3d.js`, `src/useCopilotViewer.js`,
`src/cameraPose.js`), under their card ids: `density-cell` (Atomic Density), `bond-cell` (Bond
Geometry), `pca-density` and `pca-sites` (PCA Ellipsoid), `direction-sphere` and `direction-sites`
(Displacement Directions). A turn is a 650 ms slerp of the camera (at once under reduced motion)
that lands exactly on the requested pose; in the sphere's PCA frame a crystal direction is
projected onto the site's principal axes. The host's card catalog is `src/copilotCards.js`.

A browser-loaded run is analysed in a worker of the Copilot's own, which Stop or a timeout
terminates. A server folder goes to Flask with the request cancelled on abort.

## Provider setup

- **Ollama** (local): <https://ollama.com>. Pull a model that can call tools
  (`ollama pull qwen3`, `llama3.1`, …). Start the server with this site allowed as an origin:

  ```bash
  OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve
  ```

  The Copilot sets the context length of its own requests through Ollama's native API (*The
  context window* above), so `OLLAMA_CONTEXT_LENGTH` is not needed for it.
- **LM Studio** (local): <https://lmstudio.ai>. Load a tool-capable model with a large context (the
  Copilot refuses a turn the loaded context cannot hold, and says how large to load it), enable
  **CORS** in the Developer/server settings, and start the local server
  (default `http://localhost:1234/v1`).
- **OpenAI / Gemini / Anthropic** (cloud): pick the provider and paste an API key. There is no server
  to run, but your run context and analysis results are sent to that provider. Some providers or
  networks block direct browser calls (CORS).

In the app: click the **AI Copilot** rail at the right edge (a bar at the bottom on portrait tablets and phones), open the settings (gear),
pick a provider, press **Test**, then choose a model. Settings
persist in `localStorage` (`rmc-llm-settings-v1`; the analyses switch is `copilotTools`).

HTTPS pages may call `http://localhost` in Chrome, Edge and Firefox (the trustworthy-origin
exemption); Safari may block it. Running the app locally (`npm run dev`, or the Flask/Docker mode)
always works with local providers.

## Tested with real models

The Copilot first ran against real local models on 2026-10-10: the static dev build
(`VITE_STATIC_MODE=true npm run dev`) with the **Demo** run (GaTa₄Se₈, 250 K), default settings
(temperature 0.2, analyses on), on an Apple M5 Pro with 64 GB. A Playwright script asked six
questions, each in a new conversation, and read the step cards, the answer, the Outcome chip and the
raw stream. Other jobs shared the machine and the models, so the times are upper bounds.

✓ the answer and its verdict hold; ~ the numbers hold but the answer adds a wrong claim;
✗ wrong. Times are per question, from send to the end of the answer.

| Question | qwen/qwen3.8-27b (LM Studio) | gemma4:26b (Ollama) |
|---|---|---|
| Has the run converged? | `assess_convergence` · achieved ✓ ("improving", slope −0.00052 per row) · 91 s | `assess_convergence` · achieved ✓ · 22 s |
| First Ta–Se g(r) peak; constraint-limited? | `find_gr_peaks` · achieved ✓ (2.52 Å, constraint 2.2 Å) · 140 s | no tool: the context holds it · achieved ✓ · 8 s |
| Most anharmonic site, largest U_iso | `list_sites` · achieved ✓ (ref 23, Ta: kurtosis 5.62, 0.0175 Å²) · 55 s | `list_sites` · achieved ✓ · 8 s |
| Do the Ga atoms have a preferred direction? | `analyze_displacement_directions` · achieved ✓ (no: −0.89σ, map 0.56σ) · 63 s | same call · achieved ✓ · 11 s |
| Se–Ta–Se angles against an ideal octahedron | `find_gr_peaks`, then `compute_bond_angles` at 2.33–2.71 Å; saw coordination 3.86 and widened to 2.33–2.98 Å itself (coordination 5.91). No answer within 15 min: one reasoning step took 10.5 min | before the fixes below: one call at 2.3–2.7 Å, coordination 3.75 read as "incomplete octahedra", achieved ✗ · 48 s. After: `find_gr_peaks`, then its two-peak window 2.33–2.96 Å (coordination 5.89), achieved ~ (says there are no near-180° angles; 20 % lie at 150–170°) · 95 s |
| Space group at 0.2 Å and 0.02 Å; what breaks? | `detect_symmetry` at 0.02 Å (0.2 Å from the context) · achieved ~ (F-43m and P-42₁m are right; calls P-42₁m chiral and says centrosymmetry is lost, but F-43m has none) · 321 s | no tool: the context's tolerance ladder · achieved ✓ · 45 s |

- **Tool calling works.** All 16 calls had valid JSON arguments that passed the schemas: no
  malformed call, no unknown tool, no repeated call, no budget overrun. Both models used the
  context when it already held the answer, as the prompt asks.
- **The verdict is a self-assessment.** Neither model ever wrote anything but *achieved*, even
  over a cut coordination shell. gemma4 often ends with a bare `Outcome: achieved`, without the
  reason. The step cards' checks are the part to trust.
- **Speed.** gemma4:26b (a mixture of experts) answered in 8–95 s. qwen3.8-27b, a dense reasoning
  model at LM Studio's default reasoning effort, took 1–5 min, and over 15 min on the bond-angle
  question.
- **qwen3:32b (Ollama)** was not run: other jobs held the memory it needs.

The runs found four bugs, now fixed:

1. The verdict stayed in the answer text without a chip when the model copied the prompt's
   backticks around the line (`parseOutcome` now accepts inline code and a heading marker).
2. `find_gr_peaks`'s first-peak window (2.33–2.71 Å) cut the TaSe₆ shell, whose second Ta–Se peak
   is at 2.74 Å, and no check said so. It now adds `first_two_peaks_window_A` and a note when the
   second peak overlaps the window.
3. `compute_bond_angles` said "ok" for such a window. It now warns when fewer than 60 % of the
   central atoms share the most common bond count.
4. A server error sent mid-stream (LM Studio ran out of GPU memory) was skipped, the turn ended
   empty, and the chat crashed rendering it (`turn.items is not iterable`). The error now shows,
   and an empty reply shows its note.

### A scientific review (2026-10-10, later the same day)

The same Demo and the 5 K GaNb₄Se₈ run (`data/5K_try1`, loaded as a folder) were then put through
referee-style questions with gemma4:26b, and the answers were judged against the published
structures of both compounds (GaTa₄Se₈: Phys. Rev. Research 4, 033123 (2022); GaNb₄Se₈: Phys. Rev.
B 109, 144101 (2024)). The first round was fluent and wrong in specific, fixable ways; the fixes
(new context fields, a first-shell analysis, three new analyses — `compare_fit`, `analyze_clusters`,
`analyze_modulations` — and prompt rules) are in the changelog. The third round, after the fixes:

| Question | Before | After (tools called · answer) |
|---|---|---|
| Referee review: datasets, Rwp, convergence, sampling, constraints | "two datasets", "converged — a plateau", "no peak constraint-limited" | none (the context holds it) · one measurement fitted in both spaces; still improving (−0.11 over the last tenth); 3.5 moves/atom; Ga–Se tail shaped by the 2.2 Å constraint (g at 8.7 % of the peak); Ga weakly constrained (~5 % of the signal) ✓ |
| Is the Ta₄ cluster distorted? | "peak shifted to 2.93 Å: distorted" (no tool) | `find_gr_peaks` → `analyze_clusters` with its shell window · 23 % *less* distorted than independent atoms (z = −64), breathing +27 %, weak <110> elongation and <111> compression ✓ (it read the mean edge ladder as the range of bonds — the tool now says what the ladder is) |
| Preferred Ta displacement direction; trigonal or tetragonal? | pooled map isotropic → "no distortion" | `list_sites`, pooled directions (with the new note), `analyze_clusters` · <111> compression preferred (19 % vs 14 %, z = 11), no <001> signal ✓ |
| Is GaSe₄ regular? | window 2.21–2.63 Å (CN 3.87), "± 0.02 Å" invented | `find_gr_peaks` → `compute_bond_angles` at the shell window 2.18–2.79 Å · CN 3.99, 2.434 Å (sd 0.095 measured), 109.4° (sd 4.1°), constraint caveat ✓ |
| Where does the model miss G(r)? | "the context lacks it" (partial) | `compare_fit` · 62 % of the misfit at 2.34–5.28 Å (Rw 0.024 there), largest residuals at 2.53, 3.66, 3.82, 4.9 Å ✓ |
| What can the data say about Ga? | — | none · weak weights, light sampling, constraint-shaped Ga–Se tail; trust only the average site ✓ |
| Is the average structure still cubic? | lower rungs read as symmetry breaking | none · F-43m; the P-42₁m … P1 rungs lie below the 0.031 Å noise floor ✓ |
| 5 K: does the configuration show the reported P2₁2₁2₁ Nb₄ distortion? | — | `find_gr_peaks` → `analyze_clusters` · no: clusters 7 % less distorted than independent atoms, average cluster regular (3.025–3.035 Å), inconsistent with the report ✓ |
| 5 K: converged, and where in r and Q is the misfit? | — | `assess_convergence`, `compare_fit` (G(r) only) · improving; G(r) misfit at 10.8–15 Å ~ (did not check F(Q), yet called it achieved) |
| 5 K: cubic? A hidden superstructure? | no tool could say | `read_run_data`, `analyze_modulations` · F-43m; no zone-boundary order; the strongest modes are box-scale waves ✓ |

What a referee should take from the runs themselves: neither RMC configuration carries a
cluster distortion beyond independent atoms, nor the cell doubling diffraction finds at 5 K; both
carry box-scale displacement waves; the Ta–Ta and Ga–Se distributions show the closest-approach
and shell-gap effects above. The checks on the step cards remain the part to trust.

### The view tools with a real model (2026-10-11)

The docked Copilot's view tools were checked with **gemma4:26b** (Ollama, OpenAI-compatible
endpoint) on the static build with the Demo run, in Chrome at 1440 × 900 (docked), each question
in a new conversation. Other jobs shared the model server, so the times are upper bounds.

| Question | Tools called | Result |
|---|---|---|
| First Ta–Se g(r) peak; highlight the card and mark it | `highlight_cards` (the partials card, a mark at 2.52 Å) | ✓ the card breathes, fully in view beside the panel; the plot zooms to the marked peak; the answer cites the card as a chip · 223 s |
| Largest U_iso; cite the PCA Ellipsoid card | `list_sites` | ✓ ref 23, 0.0175 Å²; from the Dashboard, its `pca-stats` chip opened PCA Ellipsoid and brought the statistics card into view · 342 s |
| The direction sphere of Ta site 1 down [111] | `list_sites`, `analyze_displacement_directions`, `orient_3d_view` (select site 1, along `[1 1 1]`) | ✓ the page opened, the sphere turned to [1 1 1] with [-1 -1 2] up (a, b, c at 120°) · 847 s |
| Centre site 1 in the site view and ring it | `orient_3d_view` (`direction-sites`, centre and ring site 1) | ✓ site 1 at the canvas centre with the blue halo; the model wrote only the Outcome · 434 s |
| What am I looking at? | `current_view` | ✓ the page, the sphere down [1 1 1] (site 1) and the site view down [3 2 3], as on screen · 553 s |

No page error and no failed step. Two display bugs were found and fixed: the Outcome note showed
LaTeX raw (`$U_{iso}$`; it now renders as Markdown inline), and chips written as the bare
reference (`<card:direction-sphere>`) read "card:direction-sphere" (they now show the card's
title; a re-run of the last question confirmed it). The model's numbers for the sphere differ
from the drawn map because the tool lets the engine choose the resolution while the page keeps
its own (ν = 10, smoothed), as noted under *Show in* in
[ai-assistant.md Step 17](../../../../docs/algorithms/ai-assistant.md#step-17--the-copilot-tool-loop-analyses-checks-budgets-verdict).

## Import boundary and extraction recipe

This module is designed to be lifted into its own repository:

- Run data and engine tools enter **only as props** to `CopilotChat` (`runName`, `plotFiles`,
  `rValueFile`, `structure`, `symmetry`, `runSettings`, `pcaSites`, `liveData`, `tools`, `pages`,
  `onOpenPage`, `runKey`, `onOpenChange`, `cards`, `onMarks`, `onLayoutChange`, `className`,
  `storageKey`, `skills`) and `WatchdogBadge` (`rValueFile`). The host renders `CopilotChat` once, as the
  last child of its app grid, and places it by its `data-mode`: a column right of the pages
  (`docked`) or a row under them (`sheet`). The page cards it may point at carry
  `data-copilot-card`. Nothing in `src/llm/` imports
  `browserData.js`, `symmetryModel.js`, `useSiteCloud.js`, or any host component. Cell math is
  deliberately duplicated from `ModelSummary.jsx`.
- The host imports only `src/llm/index.js`.
- Styling: the module's own CSS lives in `components/Copilot.css` (the conversation, settings,
  watchdog badge, the breathing `[data-copilot-focus]` edge and the card chips) and
  `components/CopilotChat.css` (the dock, rail, sheet and resize handle). A kit token it reads
  keeps its old literal as a fallback (`var(--warning, #b45309)`). It uses exactly one host UI-kit
  class, **by class name only** and never kit JS (see `src/ui/README.md`): `ui-pill`, on a step's
  Show in / Open button. An extracted copy needs a host stylesheet that defines it. The sheet's
  height is a share of `--viewport-h` (100vh without it).

To extract: copy `src/llm/`, give your app grid a column and a row for the dock, render
`CopilotChat` once with props from your own data source, tag your cards with
`data-copilot-card`, and pass your own tools (and view tools, if you have 3D views).

## Tests

`__tests__/` covers everything that doesn't need a live model:

- **Context and prompts:** context assembly and the character budget; the watchdog heuristics and
  reply parsing; the prompt builders and conversation memory.
- **Streaming:** SSE parsing, including tool-call fragments split mid-arguments and Gemini-style
  whole calls with `extra_content`; tool-call assembly and argument validation.
- **The context window** (`contextGuard.test.js`): the estimate errs high; the plan keeps a loaded
  context that fits, asks for a fitting one over a 4,096 default, keeps it across the
  conversation, caps it at the model's maximum and refuses beyond it; server detection never
  probes a cloud URL; the native NDJSON stream (thinking, content, whole tool calls, usage,
  `length`, errors, a tools-unsupported 400) and the message conversion; LM Studio's refusal; a
  reply whose server read only part of the prompt stops the turn.
  `copilotChat.test.jsx` runs a turn on a mocked Ollama loaded at 4,096 (native requests with
  `num_ctx` ≥ 16,384 and no `/v1` call), a refused turn, and a cut reply.
- **The agent loop, against a scripted model:**
  - the call → result → answer cycle; argument errors fed back and corrected;
  - `ask_user` replays on the Demo's Ta–Se shell: an ambiguous Se–Ta–Se question finds the shells,
    asks which window, and the analysis runs in the one chosen (an option, the default, a typed
    window, or Skip); a clear question runs without asking; Stop while the question waits; no time
    limit and the wait off the turn clock; one question per step, two per request, never cached;
  - unknown tools, thrown errors, duplicate calls, the per-step limit;
  - the final `tool_choice: 'none'` round and an empty last round;
  - the tools-unsupported fallback, and a schema 400 that must surface;
  - abort and timeout, which also abort the tool's signal;
  - per-result and per-turn budgets.
- **Skills:** skill files (front matter, references, a misnamed skill), `read_skill` (its budget,
  references on demand, never compacted) (`skills.test.js`); the repository's own skills parse,
  are drafts with every section, and reach the agent whole (`src/__tests__/copilotSkills.test.js`);
  one replay per skill on the Demo (`copilotSkillReplays.test.js`); the cut-shell refusal and its
  confirmed escape (`copilotTools.test.js`).
- **Checks and the Workbench's questions:** the answer checks with their negated and qualified
  phrasings and `lowerOutcome` (`answerAudit.test.js`); the clarifier on synthetic shells — split,
  merged, clean, word-like element symbols, the tolerance ladder, typed windows
  (`clarify.test.js`); and on the Demo with the real engines and context
  (`src/__tests__/copilotChecksReplays.test.js`): a contradicted answer is lowered, gemma4's
  correct answers keep *achieved*, each Workbench question (an element alone, a split shell, a
  cut window kept or replaced, one site or pooled), the clear questions that never ask, the
  model's own `ask_user` beside them within two per request, and cancel. In the dock
  (`copilotChat.test.jsx`): the Workbench's card before any model request, the lowered chip.
- **Rendering:** Outcome parsing, the Markdown pipeline with its card chips, and
  `copilotChat.test.jsx`, which drives the dock end to end against a mocked server: nothing is
  fetched before it opens, Show in hands over the page parameters, Esc keeps the conversation, a
  run change restarts it, an answer that lands while collapsed marks the rail, the open state and
  width are remembered and resize from the keyboard, and the agent's highlight_cards makes a host
  card breathe (with marks), its chip brings the card back from another page, and a click on the
  page clears it. The `ask_user` card: it waits on one request, takes focus, moves with the arrow
  keys, sends the option picked by its number key back to the model, which computes in that window;
  *Other…*; Stop while waiting; a run change drops the question; a collapsed dock marks the rail.
- **Pointing:** the focus controller (`focus.test.js`: levels, page switch and wait, clearing,
  undo callbacks, re-rendered cards), card references (`cardRefs.test.js`), `highlight_cards`
  (`viewTools.test.js`), the dock's size rules (`dockLayout.test.js`), and list/object arguments
  (`toolCalls.test.js`).

The host tools are tested on the committed demo run with the real engines in
`src/__tests__/copilotTools.test.js`; the 3D view tools and their math in `viewer3d.test.js`,
`cameraPose.test.js`, `useCopilotViewer.test.jsx`, `copilotViewTools.test.js` and
`copilotCards.test.js` (which also checks that every card in the catalog is tagged in its page),
and the plot marks in `interactivePlotMarks.test.jsx`. Run everything with `npm test` (vitest).
