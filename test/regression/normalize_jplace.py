#!/usr/bin/env python3
"""Canonicalise a jplace file for the regression tests.

The placement order inside one query, and the like-weight-ratio (a ratio of
sums accumulated in hash-map order), are not stable between runs. Everything
else - the fields list, the tree, the edge numbers, the lengths and the
likelihoods - is. This prints the file as sorted JSON with the unstable parts
normalised away, so two runs can be compared byte for byte.

Usage: normalize_jplace.py <file.jplace>
"""

import json
import sys


def canonical(document):
    # Both are metadata rather than behaviour: the invocation contains the
    # output path and the version changes with every release.
    document["metadata"]["invocation"] = "x"
    document["metadata"]["version"] = "x"
    for record in document["placements"]:
        record["n"] = sorted(record["n"])
        record["p"] = sorted(
            [
                [p[0], round(p[1], 4), round(p[2], 4), round(p[3], 4), "LWR", round(p[5], 4)]
                for p in record["p"]
            ]
        )
    document["placements"] = sorted(document["placements"], key=lambda r: r["n"])
    return document


def main():
    with open(sys.argv[1]) as handle:
        document = json.load(handle)
    json.dump(canonical(document), sys.stdout, sort_keys=True, indent=1)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
