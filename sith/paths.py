from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FieldPaths:
    root: Path
    work_images: Path
    detection: Path
    photometry: Path
    sexconf: Path
    catalogs: Path
    aperture: Path
    master: Path
    gaiaxp: Path
    calibration: Path
    diagnostics: Path
    logs: Path

    @classmethod
    def create(cls, output_root: Path, field: str) -> "FieldPaths":
        root = output_root / field
        obj = cls(
            root=root,
            work_images=root / "work" / "images",
            detection=root / "detection",
            photometry=root / "photometry",
            sexconf=root / "photometry" / "sexconf",
            catalogs=root / "photometry" / "catalogs",
            aperture=root / "aperture_correction",
            master=root / "photometry" / "master",
            gaiaxp=root / "gaiaxp",
            calibration=root / "calibration",
            diagnostics=root / "diagnostics",
            logs=root / "logs",
        )
        for p in obj.__dict__.values():
            Path(p).mkdir(parents=True, exist_ok=True)
        return obj
