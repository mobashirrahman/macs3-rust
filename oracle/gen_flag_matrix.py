#!/usr/bin/env python
"""Auto-derive the full CLI flag matrix from upstream MACS3's argparse.

The matrix is read out of `bin/macs3`'s `prepare_argparser()` -- the one function
that builds every subcommand's parser -- rather than scraped from `--help`, so it
carries the information help text drops: `nargs`, the real `dest`, whether a flag
is required, `choices`, and the default *as a parsed object* rather than a
stringified one.

Why this is worth automating: the drop-in-compatibility requirement is checked
against this table. A hand-written list drifts, and a drift in `nargs` or a
`choices` set is exactly the kind of difference that changes accept/reject
behaviour on an invalid invocation.

Why not scrape `--help`? Because MACS3's help strings *prose about flags that do
not exist*. `macs3 callpeak --help` contains `--exsize`, `--extension` and
`--shiftsize`; `macs3 randsample --help` contains `--num` and `--percent`. None of
those are option strings -- the real ones are `-e/--extsize`, `-s/--shift`,
`-n/--number`, `-p/--percentage` -- and all of the invented ones appear only inside
`help=` strings. A scraper would report them as required compatibility surface and
the port would accept flags upstream rejects.

Output: one TSV row per (subcommand, option) with a header line.

Usage:
    gen_flag_matrix.py <path-to-macs3-src> [out.tsv]
"""

import argparse
import importlib.machinery
import importlib.util
import os
import sys


def load_prepare_argparser(src):
    """Import `prepare_argparser` out of `bin/macs3`, which has no .py suffix.

    `spec_from_file_location` infers the loader from the extension and returns
    `None` for an extensionless file, so the loader is named explicitly.
    """
    path = os.path.join(src, "bin", "macs3")
    if not os.path.isfile(path):
        raise SystemExit(f"no such file: {path}")
    loader = importlib.machinery.SourceFileLoader("macs3_bin", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    # The script's module-level code only defines functions and imports, so
    # executing it does not run any command.
    loader.exec_module(mod)
    return mod.prepare_argparser


def render(value):
    """Render a default/choices value as one stable, reversible cell."""
    if isinstance(value, (list, tuple)):
        return ",".join(render(v) for v in value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    return str(value)


def action_type(action):
    """A short, stable name for what the action does."""
    if isinstance(action, argparse._StoreTrueAction):
        return "store_true"
    if isinstance(action, argparse._StoreFalseAction):
        return "store_false"
    if isinstance(action, argparse._CountAction):
        return "count"
    if isinstance(action, argparse._AppendAction):
        return "append"
    if isinstance(action, argparse._VersionAction):
        return "version"
    if isinstance(action, argparse._HelpAction):
        return "help"
    if isinstance(action, argparse._SubParsersAction):
        return "subcommand"
    if isinstance(action, argparse._StoreAction):
        return "store"
    return type(action).__name__


def nargs_of(action):
    """`nargs` as a string: '*', '+', '?', an int, or '' for exactly-one."""
    n = getattr(action, "nargs", None)
    if n is None:
        return ""
    if isinstance(n, int):
        return str(n)
    return n


def walk(parser, rows, subcommand):
    for action in parser._actions:
        if action_type(action) == "subcommand":
            for name, sub in action.choices.items():
                rows["subcommands"].add(name)
                walk(sub, rows, name)
            continue
        if not action.option_strings:
            # positional argument
            rows["out"].append(
                (
                    subcommand,
                    "<" + (action.metavar or action.dest) + ">",
                    action_type(action),
                    nargs_of(action),
                    action.dest,
                    "true" if action.required else "false",
                    render(action.default),
                    render(getattr(action, "choices", None)),
                    first_line(action.help),
                )
            )
            continue
        for flag in action.option_strings:
            rows["out"].append(
                (
                    subcommand,
                    flag,
                    action_type(action),
                    nargs_of(action),
                    action.dest,
                    "true" if action.required else "false",
                    render(action.default),
                    render(getattr(action, "choices", None)),
                    first_line(action.help),
                )
            )


def first_line(helptext):
    """The first line of help, with tabs and newlines removed for TSV safety."""
    if not helptext:
        return ""
    text = helptext.strip().split("\n")[0]
    return " ".join(text.split())


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: gen_flag_matrix.py <macs3-src> [out.tsv]")
    src = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else None

    prepare = load_prepare_argparser(src)
    parser = prepare()

    rows = {"out": [], "subcommands": set()}
    walk(parser, rows, "<top>")

    header = (
        "subcommand",
        "flag",
        "action",
        "nargs",
        "dest",
        "required",
        "default",
        "choices",
        "help",
    )
    lines = ["\t".join(header)]
    for row in rows["out"]:
        lines.append("\t".join(c.replace("\t", " ") for c in row))
    text = "\n".join(lines) + "\n"

    if out_path:
        with open(out_path, "w") as fh:
            fh.write(text)
        per_cmd = {}
        for row in rows["out"]:
            per_cmd[row[0]] = per_cmd.get(row[0], 0) + 1
        print(f"wrote {len(rows['out'])} options across "
              f"{len(rows['subcommands'])} subcommands to {out_path}")
        for name in sorted(per_cmd):
            print(f"  {name:<16} {per_cmd[name]:>3} options")
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
