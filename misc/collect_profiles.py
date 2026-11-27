#!/usr/bin/env python3
"""
collect_profiles.py -- gather the radiative-transfer density and temperature
profiles of the model grid into a single compressed file for transfer.

The profiles are expected in the directory tree

    <root>/cSD_ii/M_jj/kk/mc.denstemps.tab

where ii, jj and kk are the grid indices i, j and k of the file
bes_grid_final3.txt, each written with two digits and a leading zero.

The script writes one concatenated text file in which every profile is preceded
by a marker line

    ### NODE ID=<ID> i=<i> j=<j> k=<k> PATH=<relative path>

so that the blocks can be separated again without knowing the directory tree.
The output is gzip-compressed by default.  Missing profiles are reported and
skipped, they do not stop the run.

Examples
--------
Count what is present, without writing anything:

    python3 collect_profiles.py --grid bes_grid_final3.txt --root . --dry-run

Collect every node of the grid:

    python3 collect_profiles.py --grid bes_grid_final3.txt --root . \\
        --out profiles_final3.txt.gz

Collect a chosen subset:

    python3 collect_profiles.py --grid bes_grid_final3.txt --root . \\
        --ids 440,976,973,1016 --out profiles_subset.txt.gz
"""

import argparse
import gzip
import os
import sys


def read_grid(path):
    """Return a list of dictionaries with the identifier and the indices."""
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    header = lines[0].lstrip("# ").split()
    try:
        ci, cj, ck = header.index("i"), header.index("j"), header.index("k")
        cid = header.index("ID")
    except ValueError:
        sys.exit("ERROR: the grid file must contain the columns ID, i, j and k")
    nodes = []
    for line in lines[1:]:
        if line.lstrip().startswith("#"):
            continue
        f = line.split()
        nodes.append(dict(ID=int(float(f[cid])), i=int(float(f[ci])),
                          j=int(float(f[cj])), k=int(float(f[ck]))))
    return nodes


def profile_path(node, pad, name):
    return os.path.join("cSD_%0*d" % (pad, node["i"]),
                        "M_%0*d" % (pad, node["j"]),
                        "%0*d" % (pad, node["k"]), name)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--grid", default="bes_grid_final3.txt",
                   help="grid file listing ID, i, j and k (default: %(default)s)")
    p.add_argument("--root", default=".",
                   help="directory that contains the cSD_ii tree (default: %(default)s)")
    p.add_argument("--out", default="profiles_final3.txt.gz",
                   help="output file; gzip-compressed if it ends in .gz "
                        "(default: %(default)s)")
    p.add_argument("--name", default="mc.denstemps.tab",
                   help="name of the profile file in each node directory "
                        "(default: %(default)s)")
    p.add_argument("--pad", type=int, default=2,
                   help="number of digits of the zero-padded indices "
                        "(default: %(default)s)")
    p.add_argument("--ids", default="",
                   help="comma-separated list of node identifiers; "
                        "all nodes of the grid file if omitted")
    p.add_argument("--dry-run", action="store_true",
                   help="report how many profiles are present and stop")
    args = p.parse_args()

    nodes = read_grid(args.grid)
    if args.ids:
        wanted = {int(s) for s in args.ids.replace(",", " ").split()}
        nodes = [n for n in nodes if n["ID"] in wanted]
        missing_ids = wanted - {n["ID"] for n in nodes}
        if missing_ids:
            print("WARNING: identifiers not in the grid file: %s"
                  % ", ".join(str(i) for i in sorted(missing_ids)))

    present, missing, total_bytes = [], [], 0
    for n in nodes:
        rel = profile_path(n, args.pad, args.name)
        full = os.path.join(args.root, rel)
        if os.path.isfile(full):
            present.append((n, rel, full))
            total_bytes += os.path.getsize(full)
        else:
            missing.append((n, rel))

    print("nodes requested : %d" % len(nodes))
    print("profiles found  : %d" % len(present))
    print("profiles missing: %d" % len(missing))
    print("uncompressed size of the profiles: %.1f MB" % (total_bytes / 1.0e6))
    if missing:
        show = missing[:10]
        print("first missing paths:")
        for n, rel in show:
            print("   ID=%-6d %s" % (n["ID"], rel))
        if len(missing) > len(show):
            print("   ... and %d more" % (len(missing) - len(show)))
    if args.dry_run or not present:
        return

    opener = gzip.open if args.out.endswith(".gz") else open
    with opener(args.out, "wt") as out:
        out.write("### COLLECTION grid=%s nodes=%d\n"
                  % (os.path.basename(args.grid), len(present)))
        for n, rel, full in present:
            out.write("### NODE ID=%d i=%d j=%d k=%d PATH=%s\n"
                      % (n["ID"], n["i"], n["j"], n["k"], rel))
            with open(full) as fh:
                data = fh.read()
            out.write(data)
            if not data.endswith("\n"):
                out.write("\n")
    print("written: %s (%.1f MB)" % (args.out, os.path.getsize(args.out) / 1.0e6))


if __name__ == "__main__":
    main()
