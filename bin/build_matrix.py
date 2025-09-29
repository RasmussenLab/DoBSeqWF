#!/usr/bin/env python3
"""
Build a matrix context from pools (ground truth) and an optional decode table.

Inputs
------
1) Pools TSV (required): columns: pool_id<TAB>dimension
   - dimension must be "row" or "column"
   - pool_id must be unique

2) Decode TSV (optional): columns: sample_id<TAB>row_pool_id<TAB>col_pool_id
   - each (row_pool_id, col_pool_id) should appear at most once

Outputs
-------
A) JSON context (for machines):
   - schema_version, layout_hash
   - pools (rows/columns with index + labels)
   - matrix metadata and full cell list (row/col pool, indices, labels, alias, sample_id)

B) TSV (for humans / spreadsheets):
   - one row per intersection with: row/col pool IDs, indices, labels, cell_label (A01),
     alias (same as cell_label by default), sample_id (if any)

Usage
-----
python build_matrix.py \
  --pools pools.tsv \
  --decode decode.tsv \
  --out-json matrix_context.json \
  --out-tsv matrix_cells.tsv \
  --pad-width auto
"""
import argparse, csv, hashlib, json, sys, re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

def read_pools(pools_path: Path) -> Tuple[List[str], List[str]]:
    rows, cols = [], []
    with pools_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        first = next(reader, None)
        if first is None:
            raise ValueError("Pools TSV is empty.")
        records = [first]
        for r in reader:
            records.append(r)
        for rec in records:
            if len(rec) < 2:
                raise ValueError(f"Bad pools row (need 2 columns): {rec}")
            pool_id, dim = rec[0].strip(), rec[1].strip().lower()
            if dim not in {"row", "column"}:
                raise ValueError(f"Invalid dimension '{dim}' for pool '{pool_id}'")
            (rows if dim == "row" else cols).append(pool_id)
    # Check uniqueness
    all_ids = rows + cols
    if len(all_ids) != len(set(all_ids)):
        dup = [p for p in set(all_ids) if all_ids.count(p) > 1]
        raise ValueError(f"Duplicate pool_id(s) in pools TSV: {dup}")
    if not rows or not cols:
        raise ValueError("Need at least one row pool and one column pool.")
    return rows, cols

def read_decode(decode_path: Optional[Path]) -> Dict[Tuple[str,str], str]:
    mapping = {}
    if not decode_path: return mapping
    with decode_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")

        first = next(reader, None)
        if first is None:
            return mapping

        records = [first]
        for r in reader:
            records.append(r)
        for rec in records:
            if len(rec) < 3:
                raise ValueError(f"Bad decode row (need 3 columns): {rec}")
            sample_id, rpid, cpid = rec[0].strip(), rec[1].strip(), rec[2].strip()
            key = (rpid, cpid)
            if key in mapping and mapping[key] != sample_id:
                raise ValueError(f"Decode conflict for {key}: '{mapping[key]}' vs '{sample_id}'")
            mapping[key] = sample_id
    return mapping

def column_label_from_index(idx: int) -> str:
    """Create excel style column label 1->A, 26->Z, 27->AA"""
    label = []
    n = idx
    while n > 0:
        n -= 1
        label.append(chr(65 + (n % 26)))
        n //= 26
    return "".join(reversed(label))

def layout_hash(rows_ordered: List[str], cols_ordered: List[str]) -> str:
    m = hashlib.sha256()
    m.update(b"rows\0"); [m.update((r+"\0").encode()) for r in rows_ordered]
    m.update(b"cols\0"); [m.update((c+"\0").encode()) for c in cols_ordered]
    return m.hexdigest()[:16]

def build_context(rows: List[str],
                  cols: List[str],
                  decode: Dict[Tuple[str,str], str],
                  pad_width: Optional[int]) -> dict:
    # Force deterministic ordering
    rows_ord = sorted(rows, key=str.lower)
    cols_ord = sorted(cols, key=str.lower)

    # Indices and labels
    row_index = {pid: i for i, pid in enumerate(rows_ord)}
    col_index = {pid: i for i, pid in enumerate(cols_ord)}
    col_label = {pid: column_label_from_index(col_index[pid] +1 ) for pid in cols_ord}
    n_rows = len(rows_ord)
    width = pad_width if pad_width and pad_width > 0 else len(str(n_rows))
    row_label = {pid: str(row_index[pid] + 1).zfill(width) for pid in rows_ord}

    # Pools section
    pools_rows = [{"pool_id": pid, "index": row_index[pid], "label": row_label[pid]} for pid in rows_ord]
    pools_cols = [{"pool_id": pid, "index": col_index[pid], "label": col_label[pid]} for pid in cols_ord]

    # Cells (Cartesian product, ie. all combinations)
    cells = []
    for rpid in rows_ord:
        for cpid in cols_ord:
            ri, ci = row_index[rpid], col_index[cpid]
            rl, cl = row_label[rpid], col_label[cpid]
            cell = {
                "row_pool_id": rpid,
                "col_pool_id": cpid,
                "row_index": ri,
                "col_index": ci,
                "row_label": rl,
                "col_label": cl,
                "cell_label": f"{cl}{rl}",
                "alias": f"{cl}{rl}",
                "sample_id": decode.get((rpid, cpid), None),
            }
            cells.append(cell)

    ctx = {
        "schema_version": "1.0",
        "layout_hash": layout_hash(rows_ord, cols_ord),
        "pools": {
            "rows": pools_rows,
            "columns": pools_cols,
        },
        "matrix": {
            "n_rows": len(rows_ord),
            "n_cols": len(cols_ord),
            "row_label_width": width,
            "cells": cells,
        },
    }
    return ctx

def write_json(ctx: dict, path: Path):
    with path.open("w") as fh:
        json.dump(ctx, fh, indent=2, sort_keys=False)

def write_tsv(ctx: dict, path: Path):
    fields = [
        "row_pool_id","col_pool_id",
        "row_index","col_index",
        "row_label","col_label",
        "cell_label","alias","sample_id"
    ]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(fields)
        for cell in ctx["matrix"]["cells"]:
            w.writerow([
                cell["row_pool_id"], cell["col_pool_id"],
                cell["row_index"], cell["col_index"],
                cell["row_label"], cell["col_label"],
                cell["cell_label"], cell["alias"], cell.get("sample_id")
            ])

def write_decode(ctx: dict, path: Path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for cell in ctx["matrix"]["cells"]:
            w.writerow([
                cell["alias"], cell["row_pool_id"], cell["col_pool_id"]
            ])

def main(argv=None):
    p = argparse.ArgumentParser(description="Build matrix context from pools and optional decode.")
    p.add_argument("--pools", required=True, type=Path, help="Pools TSV (pool_id, dimension)")
    p.add_argument("--decode", required=False, type=Path, help="Decode TSV (sample_id, row_pool_id, col_pool_id)")
    p.add_argument("--out-json", required=True, type=Path, help="Output JSON path")
    p.add_argument("--out-tsv", required=True, type=Path, help="Output TSV path")
    p.add_argument("--out-decode", required=True, type=Path, help="Output decode TSV path")
    p.add_argument("--pad-width", default="auto", help='Column label zero-padding width (int or "auto")')
    args = p.parse_args(argv)

    rows, cols = read_pools(args.pools)
    decode = read_decode(args.decode) if args.decode else {}
    pad_width = None if str(args.pad_width).lower() == "auto" else int(args.pad_width)

    # Validate decode keys are known pools
    for (rpid, cpid) in decode.keys():
        if rpid not in rows:
            raise ValueError(f"Decode references unknown row_pool_id '{rpid}'")
        if cpid not in cols:
            raise ValueError(f"Decode references unknown col_pool_id '{cpid}'")

    ctx = build_context(rows, cols, decode, pad_width)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_decode.parent.mkdir(parents=True, exist_ok=True)
    write_json(ctx, args.out_json)
    write_tsv(ctx, args.out_tsv)
    write_decode(ctx, args.out_decode)

if __name__ == "__main__":
    main()
    