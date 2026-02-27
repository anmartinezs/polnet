"""Centralised logging configuration for the Polnet package.

Provides :func:`setup_logger` and the package-wide ``_LOGGER``
singleton.  Call :func:`setup_logger` once at application start-up
to attach a rotating file handler and a console handler.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

_LOGGER = logging.getLogger("polnet")


def setup_logger(
    log_folder: Path,
    console_level: int = logging.WARNING,
    file_level: int = logging.DEBUG,
):
    """Set up the logger with console and file handlers.

    Args:
        log_folder: Directory for log files.
        console_level: Minimum level for terminal output (default WARNING).
        file_level: Minimum level for file output (default DEBUG).
    """
    if _LOGGER.hasHandlers():
        _LOGGER.handlers.clear()

    log_folder.mkdir(parents=True, exist_ok=True)
    log_file = log_folder / "polnet.log"

    # Root logger level must be the lowest of the two handlers
    _LOGGER.setLevel(min(console_level, file_level))

    # --- Console Handler ---
    console_handler = logging.StreamHandler()
    console_handler.setLevel(console_level)
    console_format = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    console_handler.setFormatter(console_format)

    # --- File Handler (rotating, always verbose) ---
    file_handler = RotatingFileHandler(
        log_file, maxBytes=5_000_000, backupCount=5
    )
    file_handler.setLevel(file_level)
    file_format = logging.Formatter(
        "%(asctime)s,%(msecs)03d | %(levelname)-8s | %(name)s | %(module)s:%(lineno)d | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler.setFormatter(file_format)

    _LOGGER.addHandler(console_handler)
    _LOGGER.addHandler(file_handler)
    _LOGGER.propagate = False
