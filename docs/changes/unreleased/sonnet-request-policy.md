### Validation: Claude Sonnet 5.5's request policy, no temperature and a reply cap (protocol 2026-10-24.6)

- The user's decisions of 2026-10-24 ("Fix both first", "Raise to $20"). Claude Sonnet 5.5 takes
  only its default sampling, and a request carrying `temperature` is refused with HTTP 400. Its
  `validation/models.js` entry now sets `temperature: null` (NEBULA's setting) and
  `maxTokens: MAX_REPLY_TOKENS` (8,192). Its spend cap is $20 (RMC's share of the cycle's $40 total
  may be up to $28; $12 before), and its cache-read price is $0.20 per million tokens.
- The model port (`providerModel`) takes the two settings, only when given, so the dock's requests
  are unchanged. With `temperature: null` no request carries a temperature key, the warmer retry
  included: the retry is still asked, but unchanged. `client.js` and `ollama.js` leave a null
  temperature out of the body. `maxTokens` goes out as `max_tokens` in an OpenAI-compatible body,
  with thinking counted within it; on Ollama it caps the plan's `num_predict`. `max_tokens` is never
  sent by default, and a cap that is not a positive integer is refused when the port is made.
- The spend tracker already priced each request at 8,192 tokens of reply. That figure is now the
  cap each Sonnet request is sent with, so the reservation is a bound the API enforces. Its
  arithmetic is unchanged. Its comment now says what the charge cannot see: thinking billed as
  output but returned with no text.
- The run policy records what is sent. Sonnet's temperature is
  `{ base: null, warmer_retry_min: null, sent: false }`, and its reply cap is
  `{ sent: 'max_tokens', max_tokens: 8192 }` (`replyCapOf`). The transcript says "no temperature
  sent" and "max_tokens 8192 in every request", and the report's temperature column says "not
  sent". The LM Studio and Ollama rows, the local models' requests (byte for byte), the prompt, the
  tools and the scoring are unchanged.
- Tests (all scripted, fake servers only): Sonnet's request as sent and its body, for both the
  first request and the warmer retry: no temperature key, `max_tokens` 8192. LM Studio's and
  Ollama's bodies are pinned as before. The policy, the transcript and the report are checked from
  a metered live run, the $20 boundary of the cap check is tested, and the port, `client.js` and
  `ollama.js` have unit tests.
