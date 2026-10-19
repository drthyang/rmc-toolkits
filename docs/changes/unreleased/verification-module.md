### Copilot: one verification module answers "is this answer verified?" (architecture candidate #4)

- `src/llm/copilot/verification/` holds the evidence gate, the checks on the answer, numbers by
  reference and the computed Outcome behind one factory and two operations:
  `createVerification({ context, catalog, offered, phrasings, maxSendBacks, method })`,
  `record(step, show)` once per tool step (refs r1, r2, … in record order, none for a step with an
  error, whose facts are dropped; `show(ref)` builds the tool message and the ref resolves to
  exactly its data) and `judge(text, { canSendBack, rounds, stoppedBy, cut })` → a send-back, the
  turn's `done` event, or null for a blank reply with no answer standing. judge's order is fixed
  inside: fill the references, decide a send-back on the filled text, run the checks, add the
  unfilled-reference and missing-evidence findings, compute the Outcome. One verification serves
  one turn (record or judge after `done` throws); a phrasing for an unknown claim, or an asserted
  pattern that is not a global RegExp, is a TypeError.
- The internals are two declarative catalogues: CLAIMS (`claims.js`: one row per claim in
  send-back order — what, scope, phrasings, `backedBy` through `anyOf` / `firstOffered` /
  `allOf`, the prompt's phrase) and FACTS (`facts.js`: one row per facts key in check order — the
  producing tools, the shape, the phrasings it reads, the check; the `lowering` row has none yet).
  `evidence.js` and `answerAudit.js` are gone; `verdict.js` (with `isAnalysis`, the one
  definition of an analysis), `outcome.js`, `contextAuditFacts` and `shownContext` moved inside;
  `references.js` stays (without its unused `add`); the shell catalogue gives the finite clusters
  (`shells.js` `clusterElementsOf`, `CLUSTER_CN_MAX`). The turn builds one verification where it
  built four loop parameters, the method riding on it (both send-back wordings byte for byte);
  the agent loop records every tool step and judges every answer with it, and lost its reference
  store, step log, facts, gate state and `finish`.
- `EVIDENCE_RULES` stays written by hand; a binding test renders it from CLAIMS byte for byte
  (`evidenceSentence`) and pins what it does not say (`promptGaps`: coordination also takes
  `compute_bond_angles`, a lowering falls back to `detect_symmetry`, "no order" needs a per-site
  map and each finite cluster but neither `find_gr_peaks` nor `measure_symmetry_lowering`).
- Three behaviour changes, each pinned: the dock footer counts `tool_end.analysis` (views,
  unknown or skipped calls and a step stopped mid-run no longer count), as the computed Outcome
  does; the app's claim catalogue no longer reads the validation verdict's words — the harness
  passes them (`validation/verdict.js` `VERDICT_PHRASINGS`, the turn's `config.phrasings`), so
  its runs read as before; `DEFAULT_LIMITS` loses `maxEvidenceGates`, now the verification's
  `DEFAULT_MAX_SEND_BACKS` (2) — `config.limits.maxEvidenceGates` still sets it and the run policy
  still states it (`runPolicy.js` `TURN_LIMITS`, records unchanged). Also fixed: the dock fills an
  earlier segment's `{{ctx.…}}` from the context the model read.
- Unchanged: every send-back text, finding, computed verdict and `done` field; what the model
  sees (every model-visible snapshot); arms A and A+, the harness's expectations, the scripted
  replays and the answer corpus pins (its adapter moved onto the interface). No protocol bump.
  Tests at the interface (`verification.test.js`, `answerChecks.test.js`), the binding test
  (`verificationCatalogue.test.js`), the real tools' facts against FACTS on the Demo
  (`verificationFacts.test.js`) and the verdict's words (`validation/__tests__/verdictPhrasings.test.js`).
  Details: [the llm README](../web_app/frontend/src/llm/README.md), *Evidence, numbers by
  reference, the computed verdict*.
