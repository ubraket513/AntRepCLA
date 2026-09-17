"""Pickle-backed caching of expensive intermediates."""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

from loguru import logger


def load_cache(path: Path) -> Any | None:
    """Return the cached object at ``path``, or None if it cannot be read.

    A corrupt or unreadable cache is not fatal: the caller recomputes.
    """
    if not path.is_file():
        return None
    try:
        with open(path, 'rb') as f:
            obj = pickle.load(f)
        logger.info(f"Loaded cache from {path}")
        return obj
    except Exception as e:
        logger.warning(f"Ignoring unreadable cache {path}: {e}")
        return None


def save_cache(obj: Any, path: Path) -> None:
    """Write ``obj`` to ``path``, warning rather than raising on failure.

    Failing to cache is not a reason to discard a completed computation.
    """
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'wb') as f:
            pickle.dump(obj, f)
        logger.info(f"Saved cache to {path}")
    except Exception as e:
        logger.warning(f"Failed to write cache {path}: {e}")
