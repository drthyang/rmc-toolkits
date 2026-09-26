# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""/api/scaling/*: a malformed ``.inp`` reports its own parse error.

With the default ``kind: 'auto'`` a ``.inp`` that failed to parse was re-read
as S(Q) data: inspect answered 200 ``{'kind': 'data'}`` and preview 400
"data mode requires qmin and qmax" -- never the real problem (a truncated file,
a D-exponent, n_files = 2, a short line 22, ...). A ``.inp`` is never an S(Q)
file, so its error now surfaces; other names that merely contain "input"
(which may be data) keep the fallback.
"""

from pathlib import Path
import os
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402

RUN = ROOT / "results" / "scaling_inp_errors"
REL = "results/scaling_inp_errors"
GOOD = (
    "1\nsynth.dat\n0.60 30.0\n-9 0.1\n0\nscale.fq\nscale.gr\n25\n1000\nN\n0.05\n0\nN\nY\n1.0\n"
    "scale_ft.sq\nscale_ft.gr\n0.02\nscale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n2.48 2.65 3.1\n"
)


class MalformedInpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        shutil.rmtree(RUN, ignore_errors=True)
        RUN.mkdir(parents=True)
        (RUN / "truncated.inp").write_text("\n".join(GOOD.splitlines()[:6]) + "\n")
        (RUN / "two_files.inp").write_text(GOOD.replace("1\nsynth.dat", "2\nsynth.dat", 1))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(RUN, ignore_errors=True)

    def test_inspect_and_preview_name_the_parse_error(self):
        for name in ("truncated.inp", "two_files.inp"):
            for payload in ({"inspect": True}, {}):
                with self.subTest(name=name, **payload):
                    response = self.client.post(
                        "/api/scaling/preview", json={"path": f"{REL}/{name}", **payload}
                    )
                    body = response.get_json()
                    self.assertEqual(response.status_code, 400, body)
                    self.assertNotIn("data mode requires", body["error"])
                    self.assertNotEqual(body.get("kind"), "data")


if __name__ == "__main__":
    unittest.main()
