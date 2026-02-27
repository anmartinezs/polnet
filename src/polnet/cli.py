"""Command-line interface for polnet.

Entry point registered as ``polnet`` in pyproject.toml.

Usage::

    polnet config/all_features.yaml                     # default (WARNING on console)
    polnet config/all_features.yaml -v                  # INFO on console
    polnet config/all_features.yaml -vv                 # DEBUG on console
    polnet config/all_features.yaml --log-dir /tmp/logs
    polnet config/all_features.yaml -s 12345            # Use random seed
    polnet config/all_features.yaml -n 5                # Generate 5 tomograms
    polnet config/all_features.yaml -v -s 12345 -n 5    # INFO, seed, and multiple tomograms

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import argparse
import logging
from pathlib import Path
import shutil
import sys

import yaml

from . import __version__
from .logging_conf import (
    _LOGGER as logger,
    setup_logger,
)
from .pipeline import gen_tomos


def app() -> int:
    """Parse arguments, initialize logging, and run the tomogram generation pipeline.

    Returns:
        0 on success, 1 on failure.
    """

    parser = config_parser()
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(
            f"Error: Configuration file {config_path} does not exist.",
            file=sys.stderr,
        )
        return 1

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if "folders" not in config:
        print("Error: Config missing 'folders' section.", file=sys.stderr)
        return 1
    if "global" not in config:
        print("Error: Config missing 'global' section.", file=sys.stderr)
        return 1
    if "sample" not in config:
        print("Error: Config missing 'sample' section.", file=sys.stderr)
        return 1
    if "tem" not in config:
        print("Error: Config missing 'tem' section.", file=sys.stderr)
        return 1

    if args.log_dir is not None:
        log_dir = Path(args.log_dir)
    else:
        root = config["folders"].get("root", None)
        if root is None:
            root = Path.cwd()
            config["folders"]["root"] = str(root)

        if "output" not in config["folders"]:
            print(
                "Error: Config missing 'output' path in 'folders' section.",
                file=sys.stderr,
            )
            return 1

        log_dir = Path(root) / config["folders"]["output"]

    if args.verbose >= 2:
        console_level = logging.DEBUG
    elif args.verbose == 1:
        console_level = logging.INFO
    else:
        console_level = logging.WARNING

    setup_logger(log_folder=log_dir, console_level=console_level)

    logger.debug("Config loaded from %s", config_path)
    logger.debug("Console log level: %s", logging.getLevelName(console_level))

    if args.seed is not None:
        config["global"]["seed"] = args.seed
        logger.debug("Random seed overridden to %s", args.seed)
    if args.ntomos is not None:
        config["global"]["ntomos"] = args.ntomos
        logger.debug("Number of tomograms overridden to %s", args.ntomos)

    # Ensure root is resolved before computing the output directory
    root = config["folders"].get("root", None)
    if root is None:
        root = str(Path.cwd())
        config["folders"]["root"] = root

    out_dir = Path(root) / config["folders"]["output"]
    out_dir.mkdir(parents=True, exist_ok=True)
    _save_config_copy(config_path, out_dir, args)

    try:
        gen_tomos(config)
    except Exception as e:
        logger.error("Pipeline failed: %s", e)
        logger.debug("Traceback:", exc_info=True)
        return 1
    return 0


def _reconstruct_command(args: argparse.Namespace) -> str:
    """Rebuild the CLI invocation string from parsed arguments.

    Args:
        args: Parsed argument namespace.

    Returns:
        The reconstructed command as a single string.
    """
    parts = ["polnet", args.config]
    for _ in range(args.verbose):
        parts.append("-v")
    if args.seed is not None:
        parts.extend(["-s", str(args.seed)])
    if args.ntomos is not None:
        parts.extend(["-n", str(args.ntomos)])
    if args.log_dir is not None:
        parts.extend(["-o", args.log_dir])
    return " ".join(parts)


def _save_config_copy(
    config_path: Path, out_dir: Path, args: argparse.Namespace
) -> None:
    """Copy the YAML config into the output folder with a CLI provenance comment.

    A comment block is appended at the end of the copy recording the
    exact command that was used to launch the run, so that any output
    directory is fully self-documenting.

    Args:
        config_path: Original config file.
        out_dir: Destination output directory (must already exist).
        args: Parsed CLI arguments.
    """
    saved_config = out_dir / config_path.name
    shutil.copy2(config_path, saved_config)

    cli_cmd = _reconstruct_command(args)
    with open(saved_config, "a", encoding="utf-8") as f:
        f.write("\n# ---------------------------------------------------------\n")
        f.write(f"# Command used to generate these tomograms:\n")
        f.write(f"#   {cli_cmd}\n")
        f.write("# ---------------------------------------------------------\n")

    logger.debug("Config saved to %s", saved_config)


def config_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser for the polnet CLI.

    Returns:
        An argparse.ArgumentParser object with the defined arguments.
    """
    parser = argparse.ArgumentParser(
        description="Polnet: a comprehensive tool for generating synthetic cryo-electron tomograms."
    )

    parser.add_argument(
        "config",
        type=str,
        help="Path to the YAML configuration file.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase output verbosity. Use -v for INFO, -vv for DEBUG.",
    )

    parser.add_argument(
        "-o",
        "--log-dir",
        type=str,
        default=None,
        help="Directory for log files (default: output folder from config).",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=None,
        help="Override random seed from config (default: None, use config value).",
    )

    parser.add_argument(
        "-n",
        "--ntomos",
        type=int,
        default=None,
        help="Number of tomograms to generate, overriding config value (default: None, use config value).",
    )

    return parser


if __name__ == "__main__":
    sys.exit(app())
