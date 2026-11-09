"""Download the dataset files from Zenodo and verify them against the published checksums."""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

RECORD = "24016356"
API = f"https://zenodo.org/api/records/{RECORD}"


def _digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(target: Path) -> list[Path]:
    """Download every file of the Zenodo record into ``target`` (skips files whose MD5 already matches)."""
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(API, timeout=60) as response:
        record = json.load(response)
    paths = []
    for entry in record["files"]:
        path = target / entry["key"]
        algorithm, expected = entry["checksum"].split(":", 1)
        if not (path.exists() and _digest(path, algorithm) == expected):
            print(f"downloading {entry['key']} ({entry['size']} bytes)", flush=True)
            urllib.request.urlretrieve(entry["links"]["self"], path)
        if _digest(path, algorithm) != expected:
            raise ValueError(f"{entry['key']}: {algorithm} checksum differs from Zenodo")
        paths.append(path)
    return paths


def extract_counties(target: Path) -> int:
    """Unpack counties.zip next to it (``counties/<geoid>_<name>.csv``, as listed in counties_index.csv)."""
    import zipfile

    target = Path(target).resolve()
    with zipfile.ZipFile(target / "counties.zip") as archive:
        names = archive.namelist()
        for name in names:
            if not (target / name).resolve().is_relative_to(target):
                raise ValueError(f"counties.zip: unexpected path {name}")
        archive.extractall(target)
    return len(names)


def verify(target: Path) -> dict[str, bool]:
    """SHA-256 of every package file against manifest.json (``output_sha256``)."""
    target = Path(target)
    manifest = json.loads((target / "manifest.json").read_text(encoding="utf-8"))
    result = {name: _digest(target / name, "sha256") == digest for name, digest in manifest["output_sha256"].items()}
    if not all(result.values()):
        raise ValueError(f"SHA-256 mismatch: {[n for n, ok in result.items() if not ok]}")
    return result
