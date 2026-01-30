#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import sys
from dataclasses import dataclass
from typing import Callable, Optional


# NOTE: the resolver exports resolve_species (not resolve)
try:
    from malaria_popgen_toolkit.resources.resolver import resolve_species  # noqa: F401
except Exception:
    # Don't hard-fail CLI help if resources can't import for some reason
    resolve_species = None  # type: ignore


@dataclass(frozen=True)
class CommandSpec:
    name: str
    module: str
    help: str


# Map CLI subcommands to your src/malaria_popgen_toolkit/commands/*.py files
COMMANDS: list[CommandSpec] = [
    CommandSpec("dataset-stats", "malaria_popgen_toolkit.commands.dataset_stats", "Dataset summary statistics"),
    CommandSpec("fws-dotplot", "malaria_popgen_toolkit.commands.fws_dotplot", "Plot Fws dotplot"),
    CommandSpec("pca-plot", "malaria_popgen_toolkit.commands.pca_plot", "Plot PCA from genotype/metadata inputs"),
    CommandSpec("haplotype-map-region", "malaria_popgen_toolkit.commands.haplotype_map_region", "Plot haplotype map by region"),
    CommandSpec("ibd-plot", "malaria_popgen_toolkit.commands.ibd_plot", "Plot IBD diagnostics"),
    CommandSpec("hmmibd-matrix", "malaria_popgen_toolkit.commands.hmmibd_matrix", "Run/prepare hmmIBD matrix workflow"),
    CommandSpec("hmmibd-summary", "malaria_popgen_toolkit.commands.hmmibd_summary", "Summarize hmmIBD results"),
    CommandSpec("hmmibd-ibdplots", "malaria_popgen_toolkit.commands.hmmibd_ibdplots", "Plot hmmIBD outputs"),
    CommandSpec("ihs-selection", "malaria_popgen_toolkit.commands.ihs_selection", "Run iHS selection scan pipeline"),
    CommandSpec("xpehh-selection", "malaria_popgen_toolkit.commands.xpehh_selection", "Run XP-EHH selection scan pipeline"),
    CommandSpec("missense-drugres-af", "malaria_popgen_toolkit.commands.missense_drugres_af", "Missense drug resistance allele frequencies"),
]


def _load_command_module(module_path: str):
    try:
        return importlib.import_module(module_path)
    except Exception as e:
        raise RuntimeError(f"Failed to import command module '{module_path}': {e}") from e


def _attach_command_parser(subparsers, spec: CommandSpec) -> None:
    """
    Attach a subparser for a command module.

    We support several common patterns inside each command module:
      - add_arguments(parser) / add_args(parser) / add_parser(subparsers)
      - build_parser(parser) / configure_parser(parser)
      - run(args) / main(args) / cli(args) / entry(args)

    This makes main.py tolerant to minor refactors across modules.
    """
    mod = _load_command_module(spec.module)

    parser = subparsers.add_parser(spec.name, help=spec.help, description=spec.help)

    # 1) Let the module define its CLI args
    added = False

    # Pattern A: module provides add_arguments(parser) or similar
    for fn_name in ("add_arguments", "add_args", "configure_parser", "build_parser"):
        fn = getattr(mod, fn_name, None)
        if callable(fn):
            fn(parser)
            added = True
            break

    # Pattern B: module expects to create its own parser via add_parser(subparsers)
    if not added:
        fn = getattr(mod, "add_parser", None)
        if callable(fn):
            # Some codebases do: add_parser(subparsers) -> parser
            p = fn(subparsers)
            if p is not None:
                parser = p  # type: ignore
            added = True

    # If module didn't add args, keep parser minimal but usable
    if not added:
        parser.add_argument(
            "--help-module",
            action="store_true",
            help="Show module docstring / available callables (debug aid).",
        )

    # 2) Choose an entry function to run when this subcommand is invoked
    entry: Optional[Callable] = None
    for fn_name in ("run", "main", "cli", "entry", "execute"):
        fn = getattr(mod, fn_name, None)
        if callable(fn):
            entry = fn
            break

    if entry is None:
        # No runnable entrypoint found
        def _no_entrypoint(_args):
            raise RuntimeError(
                f"Command module '{spec.module}' does not define a runnable entrypoint.\n"
                "Expected one of: run(args), main(args), cli(args), entry(args), execute(args)."
            )
        entry = _no_entrypoint

    parser.set_defaults(_command_entrypoint=entry, _command_module=spec.module)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="malaria-popgen-toolkit: Population genomics CLI toolkit for Plasmodium analyses",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version and exit",
    )

    subparsers = parser.add_subparsers(dest="command", required=False)

    for spec in COMMANDS:
        _attach_command_parser(subparsers, spec)

    return parser


def _print_version() -> None:
    try:
        from importlib.metadata import version
        v = version("malaria-popgen-toolkit")
    except Exception:
        v = "unknown"
    print(v)


def main(argv: Optional[list[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        _print_version()
        return 0

    # No subcommand -> show help
    if not getattr(args, "command", None):
        parser.print_help()
        return 0

    # Debug helper for modules that didn't add args
    if getattr(args, "help_module", False):
        mod = _load_command_module(getattr(args, "_command_module"))
        print(f"Module: {mod.__name__}")
        print(getattr(mod, "__doc__", "") or "(no docstring)")
        print("\nCallables:")
        for name in dir(mod):
            obj = getattr(mod, name)
            if callable(obj) and not name.startswith("_"):
                print(f"  - {name}")
        return 0

    entry = getattr(args, "_command_entrypoint", None)
    if not callable(entry):
        parser.error("Internal error: no command entrypoint configured.")
        return 2

    # Run command
    try:
        ret = entry(args)
        return int(ret) if ret is not None else 0
    except KeyboardInterrupt:
        return 130
    except Exception as e:
        # Make errors readable in CLI usage
        print(f"[malaria-pipeline] ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


