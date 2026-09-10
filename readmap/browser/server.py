"""The local HTTP server behind the interactive browser.

Standard library only -- no Flask, no extra dependency to install on a machine
that already runs the analysis.  The server binds the loopback interface, serves
a fixed set of static files and answers region queries out of an already-loaded
:class:`~readmap.browser.store.BrowserStore`.

Two deliberate restrictions:

* **Loopback only.**  ``--host`` exists, but the default is ``127.0.0.1`` and
  anything else prints a warning.  This is a read-only view of the user's own
  files with no authentication; it has no business on a network interface.
* **No path joining from the request.**  Static files resolve through a fixed
  dictionary, so there is no traversal to get wrong.
"""

from __future__ import annotations

import json
import socket
import threading
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from readmap.browser.region import DEFAULT_LIMIT, build_region
from readmap.browser.store import CLASS_NAMES, EVENT_KINDS, JUNCTION_NAMES, BrowserStore
from readmap.render import theme

_STATIC_DIR = Path(__file__).with_name("static")

#: The only files that can be served, and what to send them as.  A fixed table
#: rather than a directory walk: nothing outside it is reachable by any URL.
_STATIC: dict[str, tuple[str, str]] = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/style.css": ("style.css", "text/css; charset=utf-8"),
}

#: Cap on what one request may ask for, so a hand-edited URL cannot make the
#: server build a gigabyte of JSON.
_MAX_LIMIT = 20_000
_MAX_DEPTH_BINS = 4_000


class _Handler(BaseHTTPRequestHandler):
    server_version = "readmap"
    sys_version = ""

    # Keep-alive: every response here sets Content-Length, and do_GET always
    # answers, so HTTP/1.1 is safe -- and panning fires a request per frame,
    # which would otherwise pay a TCP handshake each time.
    protocol_version = "HTTP/1.1"

    # -- plumbing ----------------------------------------------------------

    @property
    def store(self) -> BrowserStore:
        return self.server.store  # type: ignore[attr-defined]

    def log_message(self, fmt: str, *args) -> None:
        if getattr(self.server, "verbose", False):
            super().log_message(fmt, *args)

    def _send(self, body: bytes, content_type: str, status: int = 200,
              extra: dict[str, str] | None = None) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        # The page never loads anything remote; say so, so a stray URL in the
        # data cannot turn into an outbound request.
        # blob: is for the export download, which the page builds itself from a
        # fetched response.  Everything else is same-origin only, so a URL
        # sitting in a read name or a GenBank product cannot become a request.
        self.send_header(
            "Content-Security-Policy",
            "default-src 'self' blob:; img-src 'self' data:; object-src 'none'",
        )
        self.send_header("X-Content-Type-Options", "nosniff")
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, payload, status: int = 200) -> None:
        body = json.dumps(payload, allow_nan=False, separators=(",", ":")).encode("utf-8")
        self._send(body, "application/json; charset=utf-8", status)

    def _error(self, status: int, message: str) -> None:
        self._json({"error": message}, status)

    # -- routing -----------------------------------------------------------

    def do_HEAD(self) -> None:  # noqa: N802
        self.do_GET()

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        route = parsed.path
        # keep_blank_values matters: `classes=` means "the user unticked every
        # class", which is not the same as omitting the parameter.
        query = parse_qs(parsed.query, keep_blank_values=True)

        try:
            if route in _STATIC:
                return self._serve_static(route)

            handler = {
                "/api/meta": self._api_meta,
                "/api/region": self._api_region,
                "/api/depth": self._api_depth,
                "/api/read": self._api_read,
                "/api/annotation": self._api_annotation,
                "/api/search": self._api_search,
                "/api/events": self._api_events,
                "/api/export": self._api_export,
            }.get(route)

            if handler is None:
                return self._error(HTTPStatus.NOT_FOUND, f"no such route: {route}")
            return handler(query)

        except _BadRequest as exc:
            return self._error(HTTPStatus.BAD_REQUEST, str(exc))
        except (BrokenPipeError, ConnectionResetError):  # pragma: no cover
            return
        except Exception as exc:  # pragma: no cover - defensive
            return self._error(HTTPStatus.INTERNAL_SERVER_ERROR, f"{type(exc).__name__}: {exc}")

    # -- static ------------------------------------------------------------

    def _serve_static(self, route: str) -> None:
        filename, content_type = _STATIC[route]
        path = _STATIC_DIR / filename
        try:
            body = path.read_bytes()
        except OSError:
            return self._error(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                f"missing front-end file {filename}; reinstall readmap",
            )
        self._send(body, content_type)

    # -- API ---------------------------------------------------------------

    def _api_meta(self, query: dict) -> None:
        store = self.store
        meta = store.meta
        self._json(
            {
                "sample": meta.get("sample", ""),
                "created": meta.get("created", ""),
                "readmap_version": meta.get("readmap_version", ""),
                "cache": str(store.path),
                "is_genbank": meta.get("is_genbank", False),
                "has_identity": meta.get("has_identity", False),
                "contigs": meta.get("contigs", []),
                "project": meta.get("project", {}),
                "params": meta.get("params", {}),
                "sources": meta.get("sources", {}),
                "summary": meta.get("summary", {}),
                "notes": meta.get("notes", []),
                "classes": list(CLASS_NAMES),
                "junctions": list(JUNCTION_NAMES),
                "event_kinds": list(EVENT_KINDS),
                "read_count": store.read_count,
                "hit_count": store.hit_count,
                "event_count": store.event_count,
                # Colours come from the same module the PDF uses, so the screen
                # and the printed map cannot drift apart.
                "palette": {
                    "classes": {
                        cls.value: theme.CLASS_COLORS[cls]
                        for cls in theme.CLASS_ORDER
                    },
                    "class_labels": {
                        cls.value: theme.CLASS_LABELS[cls]
                        for cls in theme.CLASS_ORDER
                    },
                    "inversion": theme.INVERSION,
                    "events": theme.EVENT_COLORS,
                    "event_labels": theme.EVENT_LABELS,
                    "event_descriptions": theme.EVENT_DESCRIPTIONS,
                    "annotations": theme.ANNOTATION_COLORS,
                    "annotation_symbols": theme.ANNOTATION_SYMBOLS,
                    "annotation_legend": [
                        {"group": g, "label": label}
                        for g, label in theme.ANNOTATION_LEGEND
                    ],
                    "ink": theme.INK,
                    "ink_secondary": theme.INK_SECONDARY,
                    "ink_muted": theme.INK_MUTED,
                    "grid": theme.GRID,
                    "axis": theme.AXIS,
                    "surface": theme.SURFACE,
                },
            }
        )

    def _api_region(self, query: dict) -> None:
        contig = _one(query, "contig")
        start = _int(query, "start", 0)
        end = _int(query, "end", 0)
        limit = min(_int(query, "limit", DEFAULT_LIMIT), _MAX_LIMIT)
        classes = _classes(query)
        inverted_only = _flag(query, "inverted")
        # What the client will actually show, which is narrower than what it
        # fetches; the lane gap is measured against it.  Absent (or nonsense)
        # falls back to the fetched range.
        view_span = _int(query, "span", 0)

        try:
            view = build_region(
                self.store,
                contig,
                start,
                end,
                classes=classes,
                inverted_only=inverted_only,
                limit=limit,
                view_span=view_span,
            )
        except KeyError:
            raise _BadRequest(f"no such contig: {contig!r}") from None

        self._json(view.as_json())

    def _api_depth(self, query: dict) -> None:
        contig = _one(query, "contig")
        start = _int(query, "start", 0)
        end = _int(query, "end", 0)
        bins = min(max(1, _int(query, "bins", 800)), _MAX_DEPTH_BINS)
        if self.store.contig(contig) is None:
            raise _BadRequest(f"no such contig: {contig!r}")
        self._json(self.store.depth(contig, start, end, bins))

    def _api_read(self, query: dict) -> None:
        read_id = _int(query, "id", -1)
        try:
            self._json(self.store.read_detail(read_id))
        except IndexError:
            raise _BadRequest(f"no such read: {read_id}") from None

    def _api_annotation(self, query: dict) -> None:
        anno_id = _int(query, "id", -1)
        try:
            self._json(self.store.annotation_detail(anno_id))
        except IndexError:
            raise _BadRequest(f"no such annotation: {anno_id}") from None

    def _api_search(self, query: dict) -> None:
        text = _one(query, "q", required=False)
        self._json({"results": self.store.search(text or "")})

    def _api_events(self, query: dict) -> None:
        contig = _one(query, "contig")
        position = _int(query, "from", 0)
        direction = -1 if _one(query, "dir", required=False) == "prev" else 1
        kinds = query.get("kinds", [])
        wanted = [k for item in kinds for k in item.split(",") if k]
        event = self.store.next_event(contig, position, direction, wanted or None)
        self._json({"event": event})

    def _api_export(self, query: dict) -> None:
        contig = _one(query, "contig")
        start = _int(query, "start", 0)
        end = _int(query, "end", 0)
        fmt = (_one(query, "format", required=False) or "pdf").lower()
        if fmt not in ("pdf", "png"):
            raise _BadRequest("format must be pdf or png")

        try:
            view = build_region(
                self.store,
                contig,
                start,
                end,
                classes=_classes(query),
                inverted_only=_flag(query, "inverted"),
                limit=min(_int(query, "limit", DEFAULT_LIMIT), _MAX_LIMIT),
            )
        except KeyError:
            raise _BadRequest(f"no such contig: {contig!r}") from None

        from readmap.browser.export import render_region

        sample = str(self.store.meta.get("sample", "readmap"))
        depth = self.store.depth(view.contig, view.start, view.end, 1_200)
        body = render_region(view, sample, fmt, depth=depth)
        name = f"{sample}_{view.contig}_{view.start}-{view.end}.{fmt}"
        self._send(
            body,
            "application/pdf" if fmt == "pdf" else "image/png",
            extra={"Content-Disposition": f'attachment; filename="{_safe(name)}"'},
        )


# --------------------------------------------------------------------------
# Query helpers
# --------------------------------------------------------------------------


class _BadRequest(Exception):
    """A malformed request, reported as 400 rather than a stack trace."""


def _one(query: dict, key: str, required: bool = True) -> str:
    values = query.get(key)
    if not values:
        if required:
            raise _BadRequest(f"missing parameter: {key}")
        return ""
    return values[0]


def _int(query: dict, key: str, default: int) -> int:
    values = query.get(key)
    if not values:
        return default
    try:
        return int(float(values[0]))
    except ValueError:
        raise _BadRequest(f"{key} must be a number, got {values[0]!r}") from None


def _flag(query: dict, key: str) -> bool:
    values = query.get(key)
    return bool(values) and values[0].lower() in ("1", "true", "yes", "on")


def _classes(query: dict) -> set[str] | None:
    """The requested read classes.

    ``None`` means the parameter was absent, i.e. show everything.  An empty
    set means it was present and empty -- every filter unticked -- which must
    show no reads rather than quietly reverting to all of them.
    """
    values = query.get("classes")
    if values is None:
        return None
    wanted = {c for item in values for c in item.split(",") if c}
    unknown = wanted - set(CLASS_NAMES)
    if unknown:
        raise _BadRequest(f"unknown class(es): {', '.join(sorted(unknown))}")
    return wanted


def _safe(name: str) -> str:
    """A filename safe to put in a Content-Disposition header."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


class _Server(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, address, handler, store: BrowserStore, verbose: bool) -> None:
        self.store = store
        self.verbose = verbose
        super().__init__(address, handler)


def serve(
    store: BrowserStore,
    host: str = "127.0.0.1",
    port: int = 0,
    open_browser: bool = True,
    verbose: bool = False,
) -> int:
    """Run the browser until interrupted.  Returns a process exit code."""
    if not _STATIC_DIR.is_dir():
        print(
            f"\nERROR: the front-end files are missing from {_STATIC_DIR}.\n"
            f"       Reinstall readmap, or run from a checkout that has them.\n"
        )
        return 1

    try:
        httpd = _Server((host, port), _Handler, store, verbose)
    except OSError as exc:
        print(f"\nERROR: cannot listen on {host}:{port} -- {exc}\n")
        return 1

    actual_host, actual_port = httpd.server_address[0], httpd.server_address[1]
    if isinstance(actual_host, bytes):  # pragma: no cover
        actual_host = actual_host.decode()
    shown = "127.0.0.1" if actual_host in ("0.0.0.0", "::") else actual_host
    url = f"http://{shown}:{actual_port}/"

    summary = store.meta.get("summary", {})
    print(f"\n  readmap browser  --  {store.meta.get('sample', '')}")
    print(f"  {store.read_count:,} reads on {len(store.contigs)} contig(s), "
          f"{summary.get('events', 0):,} structural events")
    print(f"  cache: {store.path}")
    print(f"\n  {url}\n")
    if host not in ("127.0.0.1", "localhost", "::1"):
        print(
            f"  WARNING: listening on {host}, not just loopback. This server has no\n"
            f"           authentication and exposes your read data to the network.\n"
        )
    print("  Press Ctrl-C to stop.\n")

    if open_browser:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n  Stopped.\n")
    finally:
        httpd.shutdown()
        httpd.server_close()

    return 0


def free_port(host: str = "127.0.0.1") -> int:
    """An unused port, for callers that need to know it before binding."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((host, 0))
        return int(sock.getsockname()[1])
