# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path
from typing import Any, Dict

CACHE_ENV = "MALARIA_POPGEN_CACHE"
DEFAULT_CACHE = Path.home() / ".cache" / "malaria-popgen-toolkit"


def _cache_dir() -> Path:
    return Path(os.environ.get(CACHE_ENV, DEFAULT_CACHE)).expanduser().resolve()


def _load_json(filename: str) -> Dict[str, Any]:
    here = Path(__file__).parent
    with open(here / filename, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    try:
        urllib.request.urlretrieve(url, tmp)  # nosec: URLs are controlled via manifest.json
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()


def resolve_species(species: str) -> Dict[str, str]:
    """
    Resolve genome resources for a species into local cached file paths.

    Returns a dict with keys: fasta, fai, gff3, gene_product.
    Cached under ~/.cache/malaria-popgen-toolkit/<species>/...
    """
    registry = _load_json("species_registry.json")
    manifest = _load_json("manifest.json")

    if species not in registry:
        raise ValueError(
            f"Unknown species '{species}'. Available: {', '.join(sorted(registry))}"
        )

    cache = _cache_dir()
    resolved: Dict[str, str] = {}

    for logical_key, manifest_key in registry[species].items():
        if manifest_key not in manifest:
            raise KeyError(
                f"Manifest key '{manifest_key}' for species '{species}' not found in manifest.json"
            )

        entry = manifest[manifest_key]
        url = entry["url"]
        filename = entry["filename"]

        target = cache / species / filename

        if not target.exists():
            _download(url, target)

        resolved[logical_key] = str(target)

    return resolved
