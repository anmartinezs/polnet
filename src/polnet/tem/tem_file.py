"""TEM simulation configuration file loader.

Defines :class:`TEMFile`, which parses ``.tem`` key-value
configuration files and exposes microscope parameters (tilt
range, defocus, electron dose, etc.) to the TEM simulator.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import ast
from pathlib import Path


class TEMFile:
    """Parser for ``.tem`` TEM simulation configuration files.

    Reads key-value pairs from a text file (ignoring comments)
    and stores them internally for retrieval by the TEM simulator.

    Attributes:
        __params (dict): Parsed key-value configuration.
    """

    def __init__(self):
        """Initialise an empty TEM simulation configuration."""
        self.__params = {}

    @property
    def type(self):
        """Return the TEM configuration type, or None if not loaded."""
        return self.__params.get("TEM_TYPE", None)

    def load(self, in_file: Path) -> dict:
        """Load TEM simulation parameters from a ``.tem`` file.

        Args:
            in_file (Path): Path to the ``.tem`` configuration
                file.

        Returns:
            dict: Shallow copy of the parsed parameter dictionary.

        Raises:
            ValueError: If in_file does not have a ``.tem``
                extension.
            FileNotFoundError: If in_file does not exist.
        """
        if not in_file.suffix == ".tem":
            raise ValueError("Input file must have a .tem extension.")
        if not in_file.exists():
            raise FileNotFoundError(f"TEM file {in_file} does not exist.")
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
