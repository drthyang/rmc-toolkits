### Validation: `run.js live --limits` for a slow local model

- `live` takes `--limits key=value,…`, the batch's overrides of the turn's budgets (the keys of
  `TURN_LIMITS`, positive integers, milliseconds for the `*Ms` ones), e.g.
  `--limits maxTurnMs=2400000,requestMaxMs=1200000` for a local model that thinks at ~12 tokens a
  second, which the default 10-minute turn would cut short. The harness already took per-run
  overrides and kept them in each record's `policy.limits`; only the command line did not pass them.
  The same values go to every run of a batch; an unknown key or a value that is not a positive integer
  is refused before anything runs. Without the flag nothing changes: a record is byte-identical, and
  the turn, the tools, the prompt and the protocol are untouched (no protocol bump).
