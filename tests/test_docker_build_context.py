# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The Docker image builds from a clean clone and never bakes in private runs.

``data/`` is gitignored: a fresh clone has no such folder, so a ``COPY data``
fails the documented ``docker build``, and on a machine that has one it copies
private RMCProfile runs into the image. Runs are mounted at run time instead.
Static checks only; docker is not invoked.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def _instructions(dockerfile: str) -> list[str]:
    lines = []
    for raw in dockerfile.splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def _ignore_patterns() -> set[str]:
    path = ROOT / ".dockerignore"
    if not path.exists():
        return set()
    return {
        line.strip().rstrip("/")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    }


class DockerBuildContextTests(unittest.TestCase):
    def test_no_gitignored_folder_is_copied(self):
        instructions = _instructions((ROOT / "Dockerfile").read_text(encoding="utf-8"))
        copies = [line for line in instructions if line.upper().startswith("COPY")]
        for line in copies:
            if "--from=" in line:
                continue
            sources = line.split()[1:-1]
            for source in sources:
                self.assertFalse(
                    re.match(r"^(\./)?data(/|$)", source),
                    f"Dockerfile copies the gitignored data/ folder: {line}",
                )

    def test_image_has_an_empty_mount_point_for_runs(self):
        instructions = _instructions((ROOT / "Dockerfile").read_text(encoding="utf-8"))
        self.assertTrue(
            any(line.startswith("RUN") and "mkdir" in line and "/app/data" in line for line in instructions),
            "the image should create an empty /app/data for mounted runs",
        )

    def test_dockerignore_keeps_private_and_bulky_folders_out(self):
        patterns = _ignore_patterns()
        for required in ("data", ".git", ".venv", "results"):
            self.assertIn(required, patterns)
        self.assertTrue(
            {"node_modules", "**/node_modules", "web_app/frontend/node_modules"} & patterns
        )


if __name__ == "__main__":
    unittest.main()
