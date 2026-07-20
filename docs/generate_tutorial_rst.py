#!/usr/bin/env python
"""Generate RST tutorial pages from their notebook sources."""

from __future__ import print_function

import argparse
import difflib
import os
import sys

import nbformat
from nbconvert import RSTExporter


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
TUTORIALS = {
    os.path.join("notebooks", "JWST_Running_Pandexo.ipynb"): {
        "output": os.path.join("docs", "tutorialjwst.rst"),
        "title": "JWST Tutorial",
    },
}


def _strip_generated_toc(notebook):
    """Remove the notebook-only Table of Contents cell before conversion."""
    cells = list(notebook.cells)
    if cells and cells[0].cell_type == "markdown":
        source = cells[0].source
        if "Table of Contents" in source and "tocSkip" in source:
            notebook.cells = cells[1:]
    return notebook


def _demote_rst_headings(body):
    """Demote notebook headings so the docs page can add its own title."""
    demote = {
        "=": "-",
        "-": "~",
        "~": "^",
        "^": '"',
    }
    lines = body.splitlines()
    for index, line in enumerate(lines[1:], start=1):
        stripped = line.strip()
        if stripped and len(set(stripped)) == 1 and stripped[0] in demote:
            if len(stripped) >= len(lines[index - 1].strip()):
                lines[index] = demote[stripped[0]] * len(line)
    return "\n".join(lines) + ("\n" if body.endswith("\n") else "")


def _strip_trailing_whitespace(body):
    return "\n".join(line.rstrip() for line in body.splitlines()) + "\n"


def _normalize_ascii_punctuation(body):
    """Keep generated documentation consistent with the ASCII source style."""
    return body.translate(str.maketrans({
        "\u00a0": " ",
        "\u2018": "'",
        "\u2019": "'",
        "\u201c": '"',
        "\u201d": '"',
        "\u2013": "-",
        "\u2014": "-",
        "\u2026": "...",
    }))


def _render_notebook(notebook_path, title):
    notebook = nbformat.read(notebook_path, as_version=4)
    notebook = _strip_generated_toc(notebook)

    body, _resources = RSTExporter().from_notebook_node(notebook)
    body = _normalize_ascii_punctuation(body)
    body = _strip_trailing_whitespace(_demote_rst_headings(body.lstrip()))

    header = [
        ".. AUTO-GENERATED FILE. DO NOT EDIT.",
        ".. Generated from {} by docs/generate_tutorial_rst.py.".format(
            os.path.relpath(notebook_path, ROOT)
        ),
        "",
        title,
        "=" * len(title),
        "",
        "",
    ]
    return "\n".join(header) + body


def _write_or_check(notebook_relpath, config, check=False):
    notebook_path = os.path.join(ROOT, notebook_relpath)
    output_path = os.path.join(ROOT, config["output"])
    rendered = _render_notebook(notebook_path, config["title"])

    if check:
        try:
            with open(output_path, "r") as handle:
                current = handle.read()
        except IOError:
            current = ""

        if current != rendered:
            diff = difflib.unified_diff(
                current.splitlines(True),
                rendered.splitlines(True),
                fromfile=config["output"],
                tofile="generated {}".format(config["output"]),
            )
            sys.stdout.writelines(diff)
            return False
        return True

    with open(output_path, "w") as handle:
        handle.write(rendered)
    print("Generated {}".format(config["output"]))
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if generated RST differs from the checked-in files",
    )
    args = parser.parse_args(argv)

    ok = True
    for notebook_relpath, config in TUTORIALS.items():
        ok = _write_or_check(notebook_relpath, config, check=args.check) and ok

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
