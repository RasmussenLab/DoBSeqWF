#!/usr/bin/env python3
import json, sys, argparse
from pathlib import Path
from dataclasses import dataclass, fields
from pysam import VariantFile
from typing import Dict, Set, List, Tuple, Optional, Any, Union

## This is pin_basic.py - script for processing pipeline output for easy downstream analysis.
# mads - 2025-09-05

@dataclass
class VariantEntry:
    VARID: str
    CHROM: str
    POS: int
    REF: str
    ALT: str
    QUAL: float = 0.0
    AC: int = 0
    FS: float = 0.0
    SOR: float = 0.0
    GQ: int = 0
    DP: int = 0
    REF_AD: int = 0
    ALT_AD: int = 0

def get_other_variants(pools, idx):
    return set().union(*(pools[i] for i in range(len(pools)) if i != idx))

def pin(vcf_folder: Path, caller: str, matrix_context: str) -> Dict[str, Dict[str, Set[str]]]:
    """
    Performs the actual pinpointing logic and returns a dictionary of variants per sample.
    
    Parameters:
    vcf_folder (Path): Path to the folder containing pool vcfs.
    caller (str): Variant caller name.
    
    Returns:
    A dictionary mapping sample IDs to their corresponding unique variants and all pinpointable variants.
    """
    sample_variants: Dict[str, Dict[str, Set[str]]] = {}

    row_vcfs = [vcf_folder / f"{pool}.{caller}.vcf.gz" for pool in matrix_context.row_pools]
    column_vcfs = [vcf_folder / f"{pool}.{caller}.vcf.gz" for pool in matrix_context.col_pools]

    row_pools = [set(collect_variant_ids(vcf)) for vcf in row_vcfs]
    column_pools = [set(collect_variant_ids(vcf)) for vcf in column_vcfs]

    # Pinning logic:
    for cell in matrix_context._ctx["matrix"]["cells"]:
        sample_id = cell["alias"]
        row_idx = cell["row_index"]
        col_idx = cell["col_index"]
        # Extract all variants from all other pools in each dimension
        other_row_variants = get_other_variants(row_pools, row_idx)
        other_col_variants = get_other_variants(column_pools, col_idx)
        
        # Extract unique and pool-specific variants
        unique_row_variants = row_pools[row_idx].difference(other_row_variants)
        unique_col_variants = column_pools[col_idx].difference(other_col_variants)
        unique_pins = unique_row_variants.intersection(unique_col_variants)

        # Extract all pinnable variants. Only required to be unique in one dimension
        unique_one_dimension_row = unique_row_variants.intersection(column_pools[col_idx])
        unique_one_dimension_col = unique_col_variants.intersection(row_pools[row_idx])
        all_pins = unique_one_dimension_row.union(unique_one_dimension_col)

        sample_variants[sample_id] = {
            "unique_pins": unique_pins,
            "all_pins": all_pins
        }
    
    return sample_variants


def collect_variant_ids(vcf_path: str, sample_id: Optional[str] = None) -> Set[str]:
    """
    Read VCF and return set of variant IDs for a single sample in format 'chrom:pos:ref:alt'.
    """
    with VariantFile(vcf_path) as vcf:
        if sample_id is not None:
            if sample_id not in vcf.header.samples:
                raise ValueError(f"Sample '{sample_id}' not found in VCF header.")
            chosen_sample = sample_id
        else:
            if len(vcf.header.samples) == 0:
                raise ValueError("No samples found in VCF header.")
            elif len(vcf.header.samples) > 1:
                raise ValueError(
                    f"Multi-sample VCF ({len(vcf.header.samples)} samples) "
                    f"requires specifying sample_id. Available: {list(vcf.header.samples)}"
                )
            chosen_sample = next(iter(vcf.header.samples))
        
        variant_ids: Set[str] = set()
        
        for rec in vcf.fetch():
            if not rec.alts or len(rec.alts) != 1:
                raise ValueError(
                    f"Non-biallelic record at {rec.chrom}:{rec.pos} "
                    f"(REF={rec.ref}, ALTS={rec.alts}). Expected exactly one ALT."
                )
            
            # Only include if sample has a non-missing call
            sample = rec.samples[chosen_sample]
            if sample.alleles and None not in sample.alleles:
                variant_id = f"{rec.chrom}:{rec.pos}:{rec.ref}:{rec.alts[0]}"
                variant_ids.add(variant_id)
                    
        return variant_ids

def _first_scalar(x, default=0):
    if x is None:
        return default
    if isinstance(x, (list, tuple)):
        return x[0] if x else default
    return x

def collect_variant_entries(vcf_path: str, sample_id: Optional[str] = None) -> List[VariantEntry]:
    """
    Read VCF and return VariantEntry per record.

    If sample_id is given, extract that sample's fields; otherwise use the first sample if present.
    For site-only VCFs, sample-level fields remain 0.
    """
    vcf = VariantFile(vcf_path)

    # Choose sample once (or None if site-only)
    if sample_id is not None:
        if sample_id not in vcf.header.samples:
            raise ValueError(f"Requested sample_id '{sample_id}' not found in VCF header.")
        chosen_sample = sample_id
    else:
        chosen_sample = next(iter(vcf.header.samples)) if len(vcf.header.samples) > 0 else None

    out: List[VariantEntry] = []

    for rec in vcf.fetch():
        # Only accept one alt allele - multi-allelic sites is split before input.
        if not rec.alts or len(rec.alts) != 1:
            raise ValueError(
                f"Non-biallelic record at {rec.chrom}:{rec.pos} "
                f"(REF={rec.ref}, ALTS={rec.alts}). Expected exactly one ALT."
            )
        alt = rec.alts[0]

        # Site-level fields
        entry = VariantEntry(
            VARID=f"{rec.chrom}:{rec.pos}:{rec.ref}:{alt}",
            CHROM=rec.chrom,
            POS=rec.pos,
            REF=rec.ref,
            ALT=alt,
            QUAL=rec.qual or 0.0,
            AC=int(_first_scalar(rec.info.get("AC") if "AC" in vcf.header.info else 0, 0)),
            FS=float(_first_scalar(rec.info.get("FS") if "FS" in vcf.header.info else 0.0, 0.0)),
            SOR=float(_first_scalar(rec.info.get("SOR") if "SOR" in vcf.header.info else 0.0, 0.0)),
        )

        # Sample-level fields (if any)
        if chosen_sample is not None:
            s = rec.samples[chosen_sample]
            entry.GQ = int(s.get("GQ", 0) or 0)
            entry.DP = int(s.get("DP", 0) or 0)
            ad = s.get("AD", [])

            entry.REF_AD = int(ad[0]) if isinstance(ad, (list, tuple)) and len(ad) > 0 else 0
            entry.ALT_AD = int(ad[1]) if isinstance(ad, (list, tuple)) and len(ad) > 1 else 0
        out.append(entry)
    return out

def get_variant_entry(vcf_path: str, variant_id: str, sample_id: Optional[str] = None) -> VariantEntry:
    """
    Locate a specific variant, extract and return a VariantEntry.
    """
    try:
        chrom, pos_str, ref, alt = variant_id.split(":")
        pos = int(pos_str)
    except ValueError:
        raise ValueError(f"Invalid variant_id format '{variant_id}'. Expected 'chrom:pos:ref:alt'")

    with VariantFile(vcf_path) as vcf:
        if sample_id is not None:
            if sample_id not in vcf.header.samples:
                raise ValueError(f"Requested sample_id '{sample_id}' not found in VCF header.")
            chosen_sample = sample_id
        else:
            chosen_sample = next(iter(vcf.header.samples)) if len(vcf.header.samples) > 0 else None
        
        # Search for the specific variant
        for rec in vcf.fetch(chrom, pos - 1, pos):
            if (rec.pos == pos and 
                rec.ref == ref and 
                rec.alts == (alt,)):
                
                # Only accept one alt allele - multi-allelic sites should be split before input
                if not rec.alts or len(rec.alts) != 1:
                    raise ValueError(
                        f"Non-biallelic record at {rec.chrom}:{rec.pos} "
                        f"(REF={rec.ref}, ALTS={rec.alts}). Expected exactly one ALT."
                    )
                
                # Site-level fields
                entry = VariantEntry(
                    VARID=variant_id,
                    CHROM=rec.chrom,
                    POS=rec.pos,
                    REF=rec.ref,
                    ALT=rec.alts[0],
                    QUAL=rec.qual or 0.0,
                    AC=int(_first_scalar(rec.info.get("AC") if "AC" in vcf.header.info else 0, 0)),
                    FS=float(_first_scalar(rec.info.get("FS") if "FS" in vcf.header.info else 0.0, 0.0)),
                    SOR=float(_first_scalar(rec.info.get("SOR") if "SOR" in vcf.header.info else 0.0, 0.0)),
                )
                
                # Sample-level fields (if any)
                if chosen_sample is not None:
                    s = rec.samples[chosen_sample]
                    entry.GQ = int(s.get("GQ", 0) or 0)
                    entry.DP = int(s.get("DP", 0) or 0)
                    ad = s.get("AD", [])
                    
                    entry.REF_AD = int(ad[0]) if isinstance(ad, (list, tuple)) and len(ad) > 0 else 0
                    entry.ALT_AD = int(ad[1]) if isinstance(ad, (list, tuple)) and len(ad) > 1 else 0
                return entry
        raise ValueError(f"Variant '{variant_id}' not found in VCF '{vcf_path}'")

class MatrixContext:
    """
    A helper class for accessing matrix context data generated by build_matrix.py.
    
    Provides convenience functions to navigate DoBSeq matrices safely.
    """
    
    def __init__(self, ctx: dict):
        self._ctx = ctx
        
        # Sort pools by index for consistent ordering
        self._rows = sorted(ctx["pools"]["rows"], key=lambda r: r["index"])
        self._cols = sorted(ctx["pools"]["columns"], key=lambda c: c["index"])
        
        # Pool ID lists (ordered)
        self.row_pools: List[str] = [r["pool_id"] for r in self._rows]
        self.col_pools: List[str] = [c["pool_id"] for c in self._cols]
        
        # Index lookups
        self.row_index_by_pool: Dict[str, int] = {r["pool_id"]: r["index"] for r in self._rows}
        self.col_index_by_pool: Dict[str, int] = {c["pool_id"]: c["index"] for c in self._cols}
        
        # Label lookups
        self.row_label_by_pool: Dict[str, str] = {r["pool_id"]: r["label"] for r in self._rows}
        self.col_label_by_pool: Dict[str, str] = {c["pool_id"]: c["label"] for c in self._cols}
        
        # Combined dimension lookup
        self._dim_by_pool: Dict[str, str] = {
            **{p: "row" for p in self.row_pools},
            **{p: "column" for p in self.col_pools}
        }
        
        # Label lists (ordered)
        self.row_labels: List[str] = [r["label"] for r in self._rows]
        self.col_labels: List[str] = [c["label"] for c in self._cols]
        
        # Cell lookups - single source of truth
        cells = ctx["matrix"]["cells"]
        self.cells: Dict[Tuple[str, str], Dict[str, Any]] = {
            (c["row_pool_id"], c["col_pool_id"]): c for c in cells
        }
        
        # Alternative access methods
        self.by_alias: Dict[str, Dict[str, Any]] = {
            c["alias"]: c for c in cells
        }
        self.by_sample: Dict[str, Dict[str, Any]] = {
            c["sample_id"]: c for c in cells if c.get("sample_id")
        }
    
    @classmethod
    def load(cls, path: Union[Path, str]) -> "MatrixContext":
        """Load matrix context from JSON file."""
        with open(path, "r") as fh:
            data = json.load(fh)
        return cls(data)
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Return matrix shape as (n_rows, n_cols)."""
        return (self._ctx["matrix"]["n_rows"], self._ctx["matrix"]["n_cols"])
    
    @property
    def n_rows(self) -> int:
        """Number of rows in the matrix."""
        return self._ctx["matrix"]["n_rows"]
    
    @property
    def n_cols(self) -> int:
        """Number of columns in the matrix."""
        return self._ctx["matrix"]["n_cols"]
    
    @property
    def layout_hash(self) -> str:
        """Get the layout hash for this matrix configuration."""
        return self._ctx["layout_hash"]
    
    def pool_dim(self, pool_id: str) -> str:
        """Get dimension ('row' or 'column') for a pool ID."""
        if pool_id not in self._dim_by_pool:
            raise ValueError(f"Unknown pool_id: {pool_id}")
        return self._dim_by_pool[pool_id]
    
    def pool_index(self, pool_id: str) -> int:
        """Get the 1-based index for a pool ID."""
        if pool_id in self.row_index_by_pool:
            return self.row_index_by_pool[pool_id]
        elif pool_id in self.col_index_by_pool:
            return self.col_index_by_pool[pool_id]
        else:
            raise ValueError(f"Unknown pool_id: {pool_id}")
    
    def pool_label(self, pool_id: str) -> str:
        """Get the formatted label for a pool ID."""
        if pool_id in self.row_label_by_pool:
            return self.row_label_by_pool[pool_id]
        elif pool_id in self.col_label_by_pool:
            return self.col_label_by_pool[pool_id]
        else:
            raise ValueError(f"Unknown pool_id: {pool_id}")
    
    def cell(self, row_pool: str, col_pool: str) -> Optional[Dict[str, Any]]:
        """Get cell data for the intersection of row and column pools."""
        return self.cells.get((row_pool, col_pool))
    
    def cell_from_alias(self, alias: str) -> Optional[Dict[str, Any]]:
        """Get cell data from cell alias (e.g., 'A01')."""
        return self.by_alias.get(alias)
    
    def cell_from_sample(self, sample_id: str) -> Optional[Dict[str, Any]]:
        """Get cell data from sample ID."""
        return self.by_sample.get(sample_id)
    
    def cell_label(self, row_pool_id: str, col_pool_id: str) -> Optional[str]:
        """Get the cell label (e.g., 'A01') for pool intersection."""
        cell = self.cell(row_pool_id, col_pool_id)
        return cell["cell_label"] if cell else None
    
    def get_samples(self) -> List[str]:
        """Get all sample IDs in the matrix."""
        return list(self.by_sample.keys())
    
    def __repr__(self) -> str:
        return (f"MatrixContext({self.n_rows}×{self.n_cols}, "
                f"{len(self.get_samples())} samples)")

def main():
    p = argparse.ArgumentParser(description="Generate pool and pinpoint variant tables from matrix context + pool VCFs")
    p.add_argument("--context", "-c", type=Path, required=True, help="Path to matrix_context.json")
    p.add_argument("--vcf-folder", "-v", type=Path, required=True, help="Folder containing per-pool VCFs")
    p.add_argument("--output", "-o", type=Path, required=True, help="Output folder for result TSVs")
    p.add_argument("--caller", "-C", type=str, default="GATK", help="Variant caller name (used in file naming, default: GATK)")
    p.add_argument("--pin-type", "-p", type=str, choices=["unique_pins", "all_pins"], default="unique_pins", help="Which pinpoint set to use in pinpoint_variants.tsv (default: unique_pins)")
    args = p.parse_args()

    output_path = args.output
    matrix_context = args.context
    pool_vcfs_path = args.vcf_folder
    pin_type = args.pin_type
    caller = args.caller

    output_path.mkdir(parents=True, exist_ok=True)
    mc = MatrixContext.load(matrix_context)
    
    pools = (mc.row_pools,mc.col_pools)

    ######
    ## Create table with all variants:

    variant_entry_fields = [f.name for f in fields(VariantEntry)]
    header = ["dim", "pool_id", "pool_index", "pool_label", *variant_entry_fields]

    with open(output_path / 'all_pool_variants.tsv', 'w') as fout:
        print(*header, sep="\t", file=fout)
        for n, dim in enumerate(['row','column']):
            for pool_id in pools[n]:
                filepath = pool_vcfs_path / (f"{pool_id}.{caller}.vcf.gz")
                variants = collect_variant_entries(filepath)
                for v in variants:
                    values = [getattr(v, col) for col in variant_entry_fields]
                    print(dim, pool_id, mc.pool_index(pool_id), mc.pool_label(pool_id), *values, sep="\t", file=fout)
    
    ######
    # Create table with pinpointables only:
    pins = pin(vcf_folder=pool_vcfs_path, matrix_context=mc, caller=caller)
    
    common_fields = ['VARID','CHROM','POS','REF','ALT']
    variant_entry_header_fields = [f"{dim}.{f.name}" for f in fields(VariantEntry) for dim in ["row","column"] if f not in common_fields]
    header = ["uvarid", "sample_alias", "sample_id", "row_id", "row_index", "row_label", "column_id", "column_index", "column_label", *common_fields, *variant_entry_header_fields]

    with open(output_path / 'pinpoint_variants.tsv', 'w') as fout:
        print(*header, sep="\t", file=fout)
        for cell in mc._ctx["matrix"]["cells"]:
            sample_alias = cell["alias"]
            sample_id = cell["sample_id"]
            row_label = cell["row_label"]
            col_label = cell["col_label"]
            row_idx = cell["row_index"]
            col_idx = cell["col_index"]
            row_id = cell["row_pool_id"]
            col_id = cell["col_pool_id"]
            for pinpoint in pins[sample_alias][pin_type]:
                row_info = get_variant_entry(vcf_path=pool_vcfs_path / (f"{row_id}.{caller}.vcf.gz"), variant_id=pinpoint)
                col_info = get_variant_entry(vcf_path=pool_vcfs_path / (f"{col_id}.{caller}.vcf.gz"), variant_id=pinpoint)
                
                assert row_info.VARID == col_info.VARID, "Row and column variant IDs do not match"
                
                uvarid = f"{sample_alias}:{row_info.VARID}"
                print(
                    uvarid, sample_alias, sample_id, row_id, row_idx, row_label, col_id, col_idx, col_label,
                    *[getattr(row_info, col) for col in variant_entry_fields],
                    *[getattr(col_info, col) for col in variant_entry_fields if col not in common_fields],
                    sep='\t',
                    file=fout
                )

if __name__ == '__main__':
    main()