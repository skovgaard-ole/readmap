"""Input parsers.

Named ``inputs`` rather than ``io`` so there is no chance of confusion with
the standard library module of that name (decision D6).
"""

from readmap.inputs.reference import load_reference

__all__ = ["load_reference"]
