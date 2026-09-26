import { afterEach } from 'vitest'

// Vitest's worker handles its reporter RPC replies only when the event loop turns, and
// between synchronous tests the runner yields microtasks only. A file of long synchronous
// tests (the symmetry suites take 70+ s on a 2-core CI runner) therefore never lets the
// worker read the reply to its onTaskUpdate call, which fails after 60 s with
// "[vitest-worker]: Timeout calling onTaskUpdate" although every test passed.
// Yield one macrotask after each test.
afterEach(() => new Promise((resolve) => setTimeout(resolve, 0)))
