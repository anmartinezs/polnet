"""Membrane configuration file loader.

Defines :class:`MbFile`, which parses ``.mbs`` key-value
configuration files and exposes membrane type and geometry
parameters to the factory and generator pipeline.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import ast
from pathlib import Path


class MbFile:
    """Parser for ``.mbs`` membrane configuration files.

    Reads key-value pairs from a text file (ignoring comments)
    and stores them internally for retrieval by the membrane
    factory.

    Attributes:
        __params (dict): Parsed key-value configuration.
    """

    def __init__(self):
        """Initialise an empty membrane configuration."""
        self.__params = {}

    @property
    def type(self):
        """Return the membrane type string, or None if not loaded."""
        return self.__params.get("MB_TYPE", None)

    def load(self, in_file: Path) -> dict:
        """Load membrane parameters from a ``.mbs`` file.

        Args:
            in_file (Path): Path to the ``.mbs`` configuration
                file.

        Returns:
            dict: Shallow copy of the parsed parameter dictionary.

        Raises:
            ValueError: If in_file does not have a ``.mbs``
                extension.
            FileNotFoundError: If in_file does not exist.
        """
        if not in_file.suffix == ".mbs":
            raise ValueError("Input file must have a .mbs extension.")
        if not in_file.exists():
            raise FileNotFoundError(f"Membrane file {in_file} does not exist.")
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
