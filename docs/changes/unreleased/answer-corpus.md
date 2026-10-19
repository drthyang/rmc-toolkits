### Copilot: a corpus of real answers for the gate's and the checks' claim readings (plan 3.5)

- `src/__tests__/answerCorpus.test.js` runs real model answers through the evidence gate and the
  checks on the answer as they are today. The answers sit one file each in
  `src/__tests__/answerCorpus/entries/`, each with its text, source (file, lines, commit), model,
  question and run (the Demo or a dev case, never a test case). Where the record shows them, an
  entry also lists the steps that ran, and those steps must reproduce the transcript's numbers
  with the real engines.
- One adapter (`answerCorpus/readAnswer.js` `readAnswer`, public functions only) gives the reading
  each entry pins: the claims the gate reads, its send-back given the steps (or "steps
  unknown"), and the checks that fire. Each reading carries a person's judgement, so false alarms
  are counted.
- Nine answers today: gemma4:26b's three correct ones of 2026-10-11 (verbatim), Muse Glimmer's
  GaSe₄ answer from the README GIF (transcribed), and five wrong ones of the 2026-10-10 spot check
  and review, as the docs and skills quote them. Four readings fire, all right: chirality,
  inversion lost, convergence, and the gate's cluster send-back. There is no false alarm and one
  miss ("no cluster distortion" from a pooled map needs no analysis). No regex, check, prompt or
  snapshot changed.
- [copilot-validation.md](algorithms/copilot-validation.md#the-answer-corpus-plan-item-35)
  says how to add an answer from the live runs' records. A regex change must keep the pins or
  update them with a judgement.
