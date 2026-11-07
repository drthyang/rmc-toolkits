### Docs: the README, QuickStart and a demo script, current with 0.8.0

- README: what the Workbench is for (judge the fit, read the local structure, decide the next
  analysis or RMC run), a 3-step start, current screenshots from the public Demo (WebP, in
  `assets/workbench-*.webp`, replacing the older PNGs and the Copilot GIF on the page; the files
  stay for older links), *How the Copilot stays safe* and *Validation* (the 14-case held-out test,
  protocol 2026-11-01.1: Claude Sonnet 5.5 0.920 against 0.800 for the default pages and 0.948 for
  the exhaustive sweep).
- QuickStart: every section brought to 0.8.0 (the Detected SG card on the Dashboard and for
  server folders, the ladder and CIF, Save on every view) and the Copilot's connection steps: on the
  hosted site, `launchctl setenv OLLAMA_ORIGINS "https://drthyang.github.io"`, a restart of Ollama
  and Chrome's local-network prompt; or a local build on `http://localhost`, which Ollama allows by
  default; the settings drawer's Test says what failed.
- New [demo-script.md](demo-script.md): a timed 5–10 minute click path through the Demo, with a
  prep checklist, the exact prompts, one line to say per step and fallbacks.
- [demo-qa.md](demo-qa.md): asking "Se per Ta" does not avoid the coordination check's false
  *Not verified* (checked again); a number shown twice when a unit sits between it and its
  reference is listed as a known display issue.
- The Copilot guide: its five skills are interviewed (no longer drafts), and the Ollama setup names
  `launchctl` and Chrome's prompt. REFERENCE: the repository layout names the Copilot module, the
  validation harness and the skills.
