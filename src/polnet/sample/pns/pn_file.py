"""Cytosolic protein configuration file loader.

Defines :class:`PnFile`, which parses ``.pns`` key-value
configuration files and exposes monomer identity, density, and
placement parameters to the cytosolic protein network generator.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import ast
from pathlib import Path


class PnFile:
    """Parser for ``.pns`` cytosolic protein configuration files.

    Reads key-value pairs from a text file (ignoring comments)
    and stores them internally for retrieval by the protein
    placement generator.

    Attributes:
        __params (dict): Parsed key-value configuration.
    """

    def __init__(self):
        """Initialise an empty cytosolic protein configuration."""
        self.__params = {}

    @property
    def type(self):
        """Return the monomer identifier string, or None if not loaded."""
        return self.__params.get("MMER_ID", None)

    def load(self, in_file: Path) -> dict:
        """Load cytosolic protein parameters from a ``.pns`` file.

        Args:
            in_file (Path): Path to the ``.pns`` configuration
                file.

        Returns:
            dict: Shallow copy of the parsed parameter dictionary.

        Raises:
            ValueError: If in_file does not have a ``.pns``
                extension.
            FileNotFoundError: If in_file does not exist.
        """
        if not in_file.suffix == ".pns":
            raise ValueError("Input file must have a .pns extension.")
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
