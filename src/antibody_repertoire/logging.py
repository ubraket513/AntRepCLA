"""Logging setup.

The modules deliberately do not configure logging at import time. Adding a
file sink as a side effect of an import creates log files wherever the code
happens to be imported from -- including during test collection -- and makes
the sink impossible to redirect. :func:`configure_logging` is instead called
once, by the entry point.
"""

from __future__ import annotations

import sys
from pathlib import Path

from loguru import logger

LOG_FORMAT = "{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}"


def configure_logging(log_dir: Path, filename: str = "pipeline.log", level: str = "INFO") -> None:
    """Send logs to stderr and to a rotating file under ``log_dir``.

    Replaces any existing sinks, so calling this twice does not double up
    output.
    """
    logger.remove()
    logger.add(sys.stderr, format=LOG_FORMAT, level=level)

    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(
        log_dir / filename,
        format=LOG_FORMAT,
        level="DEBUG",
        rotation="10 MB",
        retention="10 days",
    )
