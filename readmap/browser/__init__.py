"""The interactive genome browser.

``readmap`` produces a poster: correct, complete, and 149 inches wide.  Finding
the handful of reads that carry the structural signal means paging through the
whole thing, and once found there is no way to ask a read what its coordinates
were.  This package answers that by serving the same classified data as a
pannable, zoomable view in a web browser.

Three pieces:

* :mod:`readmap.browser.store` -- a columnar cache of one analysed run.  Reading
  and classifying a multi-gigabyte BTOP file costs minutes; the cache turns
  every session after the first into a sub-second start.
* :mod:`readmap.browser.region` -- one visible window, packed into lanes.  Both
  the JSON API and the image exporter consume it, so they cannot disagree.
* :mod:`readmap.browser.server` -- a loopback-only HTTP server and the canvas
  front end in ``static/``.

Nothing here is imported by the PDF path, so a run that never browses never
pays for it.
"""

from __future__ import annotations

__all__ = ["serve"]


def serve(*args, **kwargs):
    """Start the browser.  See :func:`readmap.browser.server.serve`."""
    from readmap.browser.server import serve as _serve

    return _serve(*args, **kwargs)
