# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""How ``python web_app/backend/app.py`` binds: loopback, debug off, unless opted in.

The Werkzeug interactive debugger runs arbitrary Python for anyone who can
reach it, so the development server must never listen beyond this machine or
enable the debugger unless the user asks for it through the environment. These
tests exercise the resolution logic only; no server is started.
"""

from pathlib import Path
import os
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


class ServerSettingsTests(unittest.TestCase):
    def settings(self, **environ):
        return backend_app.server_settings(environ)

    def test_defaults_are_loopback_port_5000_debug_off(self):
        self.assertEqual(self.settings(), ("127.0.0.1", 5000, False))

    def test_host_is_opt_in(self):
        self.assertEqual(self.settings(RMC_TOOLKITS_HOST="0.0.0.0")[0], "0.0.0.0")
        self.assertEqual(self.settings(RMC_TOOLKITS_HOST="  ::1 ")[0], "::1")
        # A blank value is the default, never an empty bind address.
        self.assertEqual(self.settings(RMC_TOOLKITS_HOST="  ")[0], "127.0.0.1")

    def test_debug_is_opt_in(self):
        for raw in ("1", "true", "TRUE", "yes", "on"):
            self.assertTrue(self.settings(RMC_TOOLKITS_DEBUG=raw)[2], raw)
        for raw in ("", "0", "false", "no", "off"):
            self.assertFalse(self.settings(RMC_TOOLKITS_DEBUG=raw)[2], raw)
        with self.assertRaisesRegex(ValueError, "RMC_TOOLKITS_DEBUG"):
            self.settings(RMC_TOOLKITS_DEBUG="maybe")

    def test_port_prefers_port_then_rmc_toolkits_port(self):
        self.assertEqual(self.settings(RMC_TOOLKITS_PORT="5050")[1], 5050)
        self.assertEqual(self.settings(PORT="8080", RMC_TOOLKITS_PORT="5050")[1], 8080)
        for bad in ("abc", "0", "70000", "5.5"):
            with self.assertRaisesRegex(ValueError, "port"):
                self.settings(RMC_TOOLKITS_PORT=bad)

    def test_startup_line_names_the_bind_address(self):
        line = backend_app.startup_message("127.0.0.1", 5050, False)
        self.assertIn("127.0.0.1:5050", line)
        self.assertIn("debug off", line)
        exposed = backend_app.startup_message("0.0.0.0", 5000, True)
        self.assertIn("0.0.0.0:5000", exposed)
        self.assertIn("every network interface", exposed)
        self.assertIn("debug ON", exposed)

    def test_main_block_uses_the_resolved_settings(self):
        source = (ROOT / "web_app" / "backend" / "app.py").read_text(encoding="utf-8")
        main_block = source.split('if __name__ == "__main__":', 1)[1]
        self.assertNotIn("0.0.0.0", main_block)
        self.assertNotIn("debug=True", main_block)
        self.assertIn("server_settings(", main_block)


if __name__ == "__main__":
    unittest.main()
