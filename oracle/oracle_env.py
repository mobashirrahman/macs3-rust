#!/usr/bin/env python3
"""Where the pinned MACS3 oracle lives, for every script in this directory.

The oracle is not vendored: it is a git checkout of MACS3 3.0.5 plus a Python
environment with the pinned NumPy, SciPy, sklearn and hmmlearn and the compiled
extensions. `ENV.lock` records what the pin *is* -- version, commit, dependency
versions -- and deliberately not where any one machine put it, so a path baked
into a script here is a bug twice over: it does not resolve on CI, and it
silently resolves to the author's own checkout on the author's laptop, which is
how a portability problem stays invisible until a runner reports it.

Resolution order, first hit that is really there wins:

    1. the explicit argument, when a script takes one
    2. `MACS3_SRC` / `MACS3_ORACLE_BIN` / `MACS3_ORACLE_PYTHON` / `MACS3_VENV`
    3. `ENV.provisioned`, written by `provision_oracle.sh` next to this file
    4. `<repo>/.oracle`, where `provision_oracle.sh` puts a checkout by default

Import this, do not re-derive it: a script that grows its own copy of this order
is the next thing to drift.

    from oracle_env import ensure_oracle_python, require_src
    ensure_oracle_python()          # before `import numpy`
    src = require_src()             # SystemExit with a usable message if absent
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_ROOT = os.path.join(REPO, ".oracle")


def record(name: str, key: str) -> str | None:
    """`KEY=value` from a shell-style record, or None if absent or empty."""
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for line in fh:
            if line.startswith(f"{key}="):
                return line.split("=", 1)[1].strip() or None
    return None


def _first(candidates, is_there) -> str | None:
    for candidate in candidates:
        if candidate and is_there(candidate):
            return candidate
    return None


def _is_checkout(path: str) -> bool:
    return os.path.isdir(os.path.join(path, "MACS3"))


def _is_venv(path: str) -> bool:
    return os.path.isfile(os.path.join(path, "bin", "python"))


def _is_file(path: str) -> bool:
    return os.path.isfile(path)


def _checkout_of(init_path: str) -> str:
    """`MACS3_PATH` is `<checkout>/MACS3/__init__.py`, so the checkout is two up."""
    return os.path.dirname(os.path.dirname(init_path))


def oracle_src(explicit: str | None = None) -> str | None:
    """The pinned checkout, or None if this machine has none."""
    provisioned_path = record("ENV.provisioned", "MACS3_PATH")
    return _first(
        [
            explicit,
            os.environ.get("MACS3_SRC"),
            record("ENV.provisioned", "MACS3_SRC"),
            _checkout_of(provisioned_path) if provisioned_path else None,
            os.path.join(DEFAULT_ROOT, "macs3-src"),
        ],
        _is_checkout,
    )


def oracle_venv(explicit: str | None = None) -> str | None:
    """The virtualenv holding the pinned NumPy and the compiled extensions."""
    return _first(
        [
            explicit,
            os.environ.get("MACS3_VENV"),
            record("ENV.provisioned", "MACS3_VENV"),
            os.path.join(DEFAULT_ROOT, "venv"),
        ],
        _is_venv,
    )


def oracle_python(explicit: str | None = None) -> str | None:
    """An interpreter that can `import MACS3`."""
    venv = oracle_venv()
    return _first(
        [
            explicit,
            os.environ.get("MACS3_ORACLE_PYTHON"),
            os.path.join(venv, "bin", "python") if venv else None,
        ],
        _is_file,
    )


def oracle_bin(explicit: str | None = None) -> str | None:
    """The oracle's `macs3` entry point.

    The virtualenv's installed copy comes first because it is executable: upstream's
    checkout leaves `bin/macs3` mode 644, so a caller that execs this path directly
    needs the installed one. The two files differ only in their shebang.
    """
    src = oracle_src()
    venv = oracle_venv()
    return _first(
        [
            explicit,
            os.environ.get("MACS3_ORACLE_BIN"),
            os.path.join(venv, "bin", "macs3") if venv else None,
            os.path.join(src, "bin", "macs3") if src else None,
        ],
        _is_file,
    )


def require(what: str, path: str | None) -> str:
    """`path`, or exit with what to do about it."""
    if path is None:
        sys.exit(
            f"{os.path.basename(sys.argv[0])}: no MACS3 oracle ({what}); run\n"
            f"  bash oracle/provision_oracle.sh\n"
            f"or point MACS3_SRC / MACS3_ORACLE_BIN at an existing one."
        )
    return path


def require_src(explicit: str | None = None) -> str:
    return require("source tree not found", oracle_src(explicit))


def require_python(explicit: str | None = None) -> str:
    return require("no interpreter that can import MACS3", oracle_python(explicit))


def require_bin(explicit: str | None = None) -> str:
    return require("no `macs3` entry point", oracle_bin(explicit))


def ensure_oracle_python() -> None:
    """Re-exec under the provisioned interpreter, when we are not already it.

    Every script that does `import numpy` -- or imports MACS3 itself -- belongs
    on the oracle virtualenv's Python. A bare `python3` is the wrong interpreter
    twice over: on a clean runner it does not have NumPy at all, and where the
    host does have it, the version is whatever the host shipped rather than the
    pinned one, which changes results rather than raising.

    Called before the first NumPy import; re-execs with the same argv, so
    `sys.argv[0]` and every flag survive.
    """
    python = oracle_python()
    if python is None:
        print(
            f"note: no provisioned MACS3 interpreter found; continuing on "
            f"{sys.executable}. Expected NumPy results to be unverified.",
            file=sys.stderr,
        )
        return
    # A virtualenv's `bin/python` is normally a symlink to the interpreter it was
    # built from, so "same file" is the wrong test -- running the venv path is what
    # selects the venv, and comparing inodes would skip the re-exec entirely. What
    # identifies us as already being the oracle interpreter is our own prefix.
    venv = oracle_venv()
    if os.path.abspath(python) == os.path.abspath(sys.executable) or (
        venv and os.path.abspath(sys.prefix) == os.path.abspath(venv)
    ):
        return
    # Anything already printed belongs to this process, so flush before execv
    # replaces its image.
    sys.stdout.flush()
    sys.stderr.flush()
    os.execv(python, [python, *sys.argv])