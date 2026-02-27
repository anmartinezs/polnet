"""Helix filament configuration file loader.

Defines :class:`FlmsFile`, which parses ``.flms`` key-value
configuration files and exposes fiber type and geometry
parameters to the filament factory and network pipeline.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import ast
from pathlib import Path


class FlmsFile:
    """Parser for ``.flms`` helix filament configuration files.

    Reads key-value pairs from a text file (ignoring comments)
    and stores them internally for retrieval by the filament
    factory.

    Attributes:
        __params (dict): Parsed key-value configuration.
    """

    def __init__(self):
        """Initialise an empty filament configuration."""
        self.__params = {}

    @property
    def type(self):
        """Return the filament type string, or None if not loaded."""
        return self.__params.get("FLMS_TYPE", None)

    def load(self, in_file: Path) -> dict:
        """Load helix filament parameters from a ``.flms`` file.

        Args:
            in_file (Path): Path to the ``.flms`` configuration
                file.

        Returns:
            dict: Shallow copy of the parsed parameter dictionary.

        Raises:
            ValueError: If in_file does not have a ``.flms``
                extension.
            FileNotFoundError: If in_file does not exist.
        """
        if not in_file.suffix == ".flms":
            raise ValueError("Input file must have a .flms extension.")
        if not in_file.exists():
            raise FileNotFoundError(f"Flms file {in_file} does not exist.")
        with open(in_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "#" in line:
                    line = line.split("#")[0].strip()
                if "=" in line:
                    key, value = [part.strip() for part in line.split("=", 1)]
                    try:
                        self.__params[key] = ast.literal_eval(value)
                    except (ValueError, SyntaxError):
                        self.__params[key] = value
        return self.__params.copy()
