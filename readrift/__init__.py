"""readrift -- visualise long sequence reads mapped onto a reference by BLAST.

A Python port of ``read_print_23.pl``.  See ``code_structure.md`` for the
audit of the original and the design of this package, and ``CHANGES.md`` for
every way the output differs.
"""

from readrift.params import Params

__version__ = "1.0.0"
__all__ = ["Params", "main", "__version__"]


def main(argv: list[str] | None = None) -> int:
    """Entry point.  Kept import-light so ``--help`` does not load matplotlib."""
    import sys

    from readrift.cli import SUBCOMMANDS, parse_args, parse_browse_args

    raw = list(sys.argv[1:] if argv is None else argv)

    if raw and raw[0] in SUBCOMMANDS:
        from readrift.pipeline import browse

        return browse(parse_browse_args(raw[1:]))

    from readrift.pipeline import run

    return run(parse_args(raw))
