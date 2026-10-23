"""Data frames from mapping batches (``MappingBatch.to_pandas``)."""

from __future__ import annotations

import weakref

import numpy as np


# Per index: the transcript names as a categorical dtype (or an object array when names repeat),
# built once; an index holds ~250,000 transcripts.
_names = weakref.WeakKeyDictionary()


def _transcripts(index, ids):
    """Transcript names for ``ids`` (-1: none), categorical over every transcript in the index."""
    import pandas as pd

    names = _names.get(index)
    if names is None:
        tx = index.transcript_names()
        try:
            names = pd.CategoricalDtype(categories=tx)
        except ValueError:
            names = np.asarray(tx, dtype=object)
        _names[index] = names
    if isinstance(names, pd.CategoricalDtype):
        return pd.Categorical.from_codes(ids, dtype=names)
    out = names[np.maximum(ids, 0)]
    out[ids < 0] = None
    return out


def batch_to_pandas(batch, unmapped=False):
    import pandas as pd

    hits = batch.hit_columns()
    reads = batch.read_columns()
    row_read = hits.pop("read").astype(np.int64)
    tx = hits.pop("transcript_id").astype(np.int64)
    if unmapped:
        # One row per unmapped read, in read order among the hits.
        miss = np.flatnonzero(reads["num_hits"] == 0)
        if len(miss):
            order = np.argsort(np.concatenate([row_read, miss]), kind="stable")
            n = len(miss)
            row_read = np.concatenate([row_read, miss])[order]
            tx = np.concatenate([tx, np.full(n, -1)])[order]
            for k, v in hits.items():
                fill = np.zeros(n, dtype=v.dtype)
                hits[k] = np.concatenate([v, fill])[order]
    names = batch.names
    read = np.asarray(names, dtype=object)[row_read] if names is not None else row_read
    df = pd.DataFrame({"read": read, "read_length": reads["length"][row_read]})
    df["rank"] = hits.pop("rank")
    df["transcript"] = _transcripts(batch.index, tx)
    df["transcript_id"] = tx
    reverse = hits.pop("reverse")
    df["strand"] = pd.Categorical(np.where(reverse, "-", "+"), categories=["+", "-"])
    for k in ("score", "query_start", "query_end", "target_start", "target_end", "anchors",
              "anchored", "max_discrepancy", "extended"):
        df[k] = hits.pop(k)
    df["low_coverage"] = reads["low_coverage"][row_read]
    df["rescued"] = reads["rescued"][row_read]
    if unmapped:
        none = tx < 0
        df.loc[none, "strand"] = None
        df["mapped"] = ~none
    return df
