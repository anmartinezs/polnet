"""Membrane-bound protein configuration file loader.

Defines :class:`PmFile`, which parses ``.pms`` key-value
configuration files and exposes monomer identity, density, and
placement parameters to the membrane-bound protein generator.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import ast
from pathlib import Path


class PmFile:
    """Parser for ``.pms`` membrane-bound protein config files.

    Reads key-value pairs from a text file (ignoring comments)
    and stores them internally for retrieval by the protein
    placement generator.

    Attributes:
        __params (dict): Parsed key-value configuration.
    """

    def __init__(self):
        """Initialise an empty membrane-bound protein configuration."""
        self.__params = {}

    @property
    def type(self):
        """Return the monomer identifier string, or None if not loaded."""
        return self.__params.get("MMER_ID", None)

    def load(self, in_file: Path) -> dict:
        """Load membrane-bound protein parameters from a ``.pms`` file.

        Args:
            in_file (Path): Path to the ``.pms`` configuration
                file.

        Returns:
            dict: Shallow copy of the parsed parameter dictionary.

        Raises:
            ValueError: If in_file does not have a ``.pms``
                extension.
            FileNotFoundError: If in_file does not exist.
        """
        if not in_file.suffix == ".pms":
            raise ValueError("Input file must have a .pms extension.")
        if not in_file.exists():
            raise FileNotFoundError(f"Protein file {in_file} does not exist.")
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
