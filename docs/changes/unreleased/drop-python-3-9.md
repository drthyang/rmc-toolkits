### Python 3.10 or newer

- `rmc_toolkits` now needs Python 3.10 or newer (`requires-python = ">=3.10"`). Python 3.9 reached
  end of life, and GitHub's Ubuntu 26.04 runners no longer provide it, so CI could not test it.
- CI's oldest-dependencies job runs on Python 3.10 with the same floors pinned exactly (numpy
  1.22.4, SciPy 1.8.1, matplotlib 3.6.3, contourpy 1.0.7). The floors are unchanged.
