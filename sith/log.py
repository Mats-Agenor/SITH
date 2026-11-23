from __future__ import annotations

import datetime as _dt
import shlex
import subprocess
from pathlib import Path


class Logger:
    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("", encoding="utf-8")

    def write(self, message: str = "") -> None:
        stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(message, flush=True)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(f"[{stamp}] {message}\n")

    def command(self, cmd: list[str], cwd: Path | None = None) -> None:
        self.write("$ " + " ".join(shlex.quote(str(x)) for x in cmd))
        proc = subprocess.run(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if proc.stdout:
            print(proc.stdout, end="", flush=True)
            with self.path.open("a", encoding="utf-8") as fh:
                fh.write(proc.stdout)
                if not proc.stdout.endswith("\n"):
                    fh.write("\n")
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
