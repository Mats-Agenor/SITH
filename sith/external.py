from __future__ import annotations

import os
import shutil
from pathlib import Path


def discover_executable(explicit: str | None, candidates: list[str], label: str) -> str:
    if explicit and explicit.lower() != "auto":
        p = Path(explicit).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p.resolve())
        found = shutil.which(explicit)
        if found:
            return found
        raise FileNotFoundError(f"Configured {label} executable not found: {explicit}")

    for name in candidates:
        found = shutil.which(name)
        if found:
            return found

    raise FileNotFoundError(
        f"{label} executable not found. Tried: {', '.join(candidates)}. "
        "Install it or set its path in JEDIS.conf."
    )


def discover_sextractor_config_dir(explicit: str | None) -> Path:
    candidates: list[Path] = []
    if explicit and explicit.lower() != "auto":
        candidates.append(Path(explicit).expanduser())
    candidates.extend(
        [
            Path("/usr/share/source-extractor"),
            Path("/usr/local/share/source-extractor"),
            Path("/usr/share/sextractor"),
            Path("/usr/local/share/sextractor"),
        ]
    )

    for d in candidates:
        if (
            d.is_dir()
            and (d / "default.nnw").is_file()
            and (d / "tophat_3.0_3x3.conv").is_file()
        ):
            return d.resolve()

    raise FileNotFoundError(
        "SExtractor was found, but the exact SPLUSCALIB resources "
        "default.nnw and tophat_3.0_3x3.conv were not found. "
        "Set sextractor_config_dir in JEDIS.conf."
    )


def package_root() -> Path:
    return Path(__file__).resolve().parent


def stilts_jar() -> Path:
    p = package_root() / "external" / "stilts.jar"
    if not p.is_file():
        raise FileNotFoundError(f"Bundled STILTS jar not found: {p}")
    return p
