"""krepp: k-mer-based search, distance estimation and phylogenetic placement.

In-memory results — no stdout parsing, no intermediate files:

    import krepp

    opts = krepp.Options(num_threads=4, p_value=True)
    s = krepp.Session("path/to/index", opts)

    for qd in s.dist("queries.fasta"):          # like `krepp dist`
        for row in qd.rows:
            print(qd.query, row.ref_name, row.dist, row.p_value)

    for qp in s.place("queries.fasta"):         # like `krepp place`
        for p in qp.placements:
            print(qp.query, p.edge_num, p.distal_node, p.like_weight_ratio)

    for qm in s.matches("queries.fasta"):       # k-mer match summaries
        print(qm.query, qm.onmers, qm.enmers, len(qm.rows))

Build with `make python` from the repository root.

Errors: invalid paths and option values raise ValueError; corrupt indexes,
unreadable/truncated query files and query-time record errors raise
RuntimeError (they are relayed out of the parallel regions, not swallowed).
Session calls are serialised process-wide; OpenMP parallelism within one
call is unaffected. With num_threads > 1 the order of per-query results
across batches is nondeterministic, exactly like the CLI's output order.
"""

from . import _krepp
from ._krepp import (
    DistRow,
    MatchRow,
    Placement,
    QueryDist,
    QueryMatch,
    QueryPlace,
    Session,
    __version__,
)


class Options(_krepp.Options):
    """Query options; field names and defaults match the CLI flags.

    Accepts keyword arguments for convenience:
        krepp.Options(num_threads=4, p_value=True)
    """

    def __init__(self, **kwargs):
        super().__init__()
        for key, value in kwargs.items():
            if not hasattr(type(self), key) and not hasattr(self, key):
                raise TypeError(f"Options has no field '{key}'")
            setattr(self, key, value)


__all__ = [
    "Session",
    "Options",
    "DistRow",
    "QueryDist",
    "Placement",
    "QueryPlace",
    "MatchRow",
    "QueryMatch",
    "__version__",
]
