### Tests: the χ² clamp compared to 1e-12 relative (oldest-deps CI)

- `test_log_data_carries_every_chi_term_and_their_total` compared the zero curvature term's clamp,
  `np.log(1e-12)`, with `math.log(1e-12)` bit for bit. Under the oldest pinned numpy its SIMD log
  can differ in the last ulp, depending on the runner's CPU, which turned main's Python 3.9 job red.
  The test now checks the same things exactly — one value on every row, the right row count — and
  the clamp's value to 1e-12 relative. The code is unchanged.
