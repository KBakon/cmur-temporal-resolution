"""Verify every archived deposit file listed in CHECKSUMS.sha256."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    lines = (ROOT / "CHECKSUMS.sha256").read_text(encoding="utf-8").splitlines()
    seen = set()
    for line in lines:
        expected, relative = line.split("  ", 1)
        if relative in seen or relative.startswith("/") or ".." in Path(relative).parts:
            raise ValueError("Invalid or duplicate checksum path: " + relative)
        seen.add(relative)
        path = ROOT / relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"Checksum mismatch: {relative}")
    actual_files = {
        p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file()
        and p.name != "CHECKSUMS.sha256" and "__pycache__" not in p.parts
    }
    if seen != actual_files:
        raise AssertionError(f"Unlisted files: {sorted(actual_files-seen)}; missing files: {sorted(seen-actual_files)}")
    print(f"PASS: {len(seen)} deposit-file SHA-256 checks")


if __name__ == "__main__":
    main()
