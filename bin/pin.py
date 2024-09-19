#!/usr/bin/env python3

from pathlib import Path
from typing import List, Set, Dict, Tuple
from tempfile import TemporaryDirectory
import subprocess
import argparse
import logging

## This is pin.py - pinning script for the DoBSeq pipeline.
# mads - 2024-07-03
#
# Requires:
# - bcftools
# - htslib

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class VcfMerger:
    def __init__(self, input_files: List[Path], output_prefix: str, output_folder: Path, tmpdir: TemporaryDirectory):
        self.input_files: List[Path] = input_files
        self.output_prefix: str = output_prefix
        self.output_folder: Path = output_folder
        self.tmpdir: TemporaryDirectory = tmpdir
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def drop_pl_tag(self) -> None:
        for input_file in self.input_files:
            cmd = ["bcftools", "annotate", "-Oz", "-x", "FORMAT/PL", "-o", str(self.tmpdir / (input_file.stem + "_annotated.vcf.gz")), str(input_file)]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                raise Exception("Error running bcftools annotate for individual.")

    def tabix(self, input_files: List[Path]) -> None:
        for input_file in input_files:
            cmd = ["tabix", str(input_file)]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                raise Exception(f"Error running bcftools tabix. Input file: {input_file}")

    def merge_files(self, force_samples: bool = False) -> None:
        merge_method = "AC:join,AF:join,AN:join,DP:join,FS:join,MLEAC:join,MLEAF:join,MQ:join,MQRankSum:join,QD:join,ReadPosRankSum:join,SOR:join"
        if force_samples:
            # Force merge of all samples AND avoid merge of multiallelic sites.
            cmd = ["bcftools", "merge", "--info-rules", merge_method, "-o", str(self.tmpdir / "tmp_merge.vcf"), "--merge", "none", "--force-samples"] + [str(input_file) for input_file in self.input_files]
        else:
            cmd = ["bcftools", "merge", "--info-rules", merge_method, "-o", str(self.output_folder / f"{self.output_prefix}.vcf.gz")] + [str(input_file) for input_file in self.input_files]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            raise Exception("Error running bcftools isec for DB.")

    def drop_genotypes_view(self) -> None:
        cmd = ["bcftools", "view", "--drop-genotypes", "-o", str(self.output_folder / f"{self.output_prefix}.vcf.gz"), str(self.tmpdir / "tmp_merge.vcf")]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            raise Exception("Error running bcftools view.")

    def merge_individuals(self) -> None:
        self.drop_pl_tag()
        self.input_files = [self.tmpdir / (input_file.stem + "_annotated.vcf.gz") for input_file in self.input_files]
        self.tabix(self.input_files)
        self.merge_files(force_samples=False)
        self.tabix([self.output_folder / f"{self.output_prefix}.vcf.gz"])
    
    def merge_matrix(self) -> None:
        self.merge_files(force_samples=True)
        self.drop_genotypes_view()

def filter_vcf(vcf_file: str, sample_variants: Set[str], output_file: Path) -> None:
    """
    Filters a VCF file based on allele specific variants.
    
    Parameters:
    vcf_file (str): Path to the input VCF file.
    sample_variants (set): A set of variants to filter.
    output_file (Path): Path to the output filtered VCF file.
    """
    variants = sample_variants.copy()
    with open(vcf_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
            else:
                chrom, pos, _, ref, alt = line.split('\t')[:5]
                
                if f"{chrom}:{pos}:{ref}:{alt}" in variants:
                    fout.write(line)
                    variants.remove(f"{chrom}:{pos}:{ref}:{alt}")
    if variants:
        raise ValueError("Not all sample variants were found in the VCF file.")


def get_variants_from_pool(variant_table: Path) -> Set[str]:
    """
    Takes a path to a variant table and returns a set of variants.
    Each variant is represented as a string formatted as "CHROM:POS:REF:ALT".
    
    Parameters:
    variant_table (Path): Path to a variant table.
    
    Returns:
    set: A set containing variants.
    """
    with open(variant_table, 'r') as fin:
        header = fin.readline().strip().split('\t')
        if header[0:5] != ['CHROM', 'POS', 'ID', 'REF', 'ALT']:
            raise ValueError("Variant table is in unexpected format.")
        
        # Construct a set of variants with the format "CHROM:POS:ALT"
        variants = set()
        for line in fin:
            var_info = line.strip().split('\t')
            CHROM,POS,_,REF,ALT = var_info[0:5]
            variants.add(f"{CHROM}:{POS}:{REF}:{ALT}")

    return variants

def decode(decodetable_path: Path, vartable_folder: Path, caller: str) -> Tuple[Dict[str, Tuple[int, int]], Dict[str, Tuple[str, str]], List[Path], List[Path]]:
    """
    Parses the decodetable file and returns a dictionary mapping sample IDs to their corresponding horizontal and vertical pool indices and pool ids.
    
    Parameters:
    decodetable_path (Path): Path to the decodetable file.
    
    Returns:
    sample_pool_map: A dictionary mapping sample IDs to their corresponding horizontal and vertical pool indices.
    sample_pool_ids: A dictionary mapping sample IDs to their corresponding horizontal and vertical pool ids.
    horizontal_vartables: A list of horizontal variant tables.
    vertical_vartables: A list of vertical variant tables.
    """

    sample_pool_map: Dict[str, Tuple[int, int]] = {}
    sample_pool_ids: Dict[str, Tuple[str, str]] = {}

    with open(decodetable_path, "r") as fin:
        horizontal_pools = set()
        vertical_pools = set()
        for line in fin:
            _, h, v = line.strip().split("\t")
            horizontal_pools.add(h)
            vertical_pools.add(v)
        horizontal_pools = sorted(list(horizontal_pools))
        vertical_pools = sorted(list(vertical_pools))

    with open(decodetable_path, "r") as fin:
        for line in fin:
            id, h, v = line.strip().split("\t")
            sample_pool_map[id] = (horizontal_pools.index(h), vertical_pools.index(v))
            sample_pool_ids[id] = (h,v)

    horizontal_vartables = [vartable_folder / Path(f"{pool}.{caller}.tsv") for pool in horizontal_pools]
    vertical_vartables = [vartable_folder / Path(f"{pool}.{caller}.tsv") for pool in vertical_pools]

    return sample_pool_map, sample_pool_ids, horizontal_vartables, vertical_vartables

def get_other_variants(pools, idx):
    return set().union(*(pools[i] for i in range(len(pools)) if i != idx))

def pin(vartable_folder: Path, decodetable_path: Path, caller: str) -> Dict[str, Dict[str, Set[str]]]:
    """
    Performs the actual pinpointing logic and returns a dictionary of variants per sample.
    
    Parameters:
    vartable_folder (Path): Path to the folder containing variant tables.
    decodetable_path (Path): Path to the decodetable file.
    caller (str): Caller name.
    
    Returns:
    A dictionary mapping sample IDs to their corresponding unique variants and all pinpointable variants.
    """
    sample_variants: Dict[str, Dict[str, Set[str]]] = {}

    sample_pool_map, _, horizontal_vartables, vertical_vartables = decode(decodetable_path, vartable_folder, caller)
    
    horizontal_pools = [set(get_variants_from_pool(vartable)) for vartable in horizontal_vartables]
    vertical_pools = [set(get_variants_from_pool(vartable)) for vartable in vertical_vartables]

    # Pinning logic:
    for sample, (h_idx, v_idx) in sample_pool_map.items():

        # Extract all variants from all other pools in each dimension
        other_h_variants = get_other_variants(horizontal_pools, h_idx)
        other_v_variants = get_other_variants(vertical_pools, v_idx)
        
        # Extract unique and pool-specific variants
        unique_h_variants = horizontal_pools[h_idx].difference(other_h_variants)
        unique_v_variants = vertical_pools[v_idx].difference(other_v_variants)
        unique_pins = unique_h_variants.intersection(unique_v_variants)

        # Extract all pinnable variants. Only required to be unique in one dimension
        unique_one_dimension_h = unique_h_variants.intersection(vertical_pools[v_idx])
        unique_one_dimension_v = unique_v_variants.intersection(horizontal_pools[h_idx])
        all_pins = unique_one_dimension_h.union(unique_one_dimension_v)

        sample_variants[sample] = {
            "unique_pins": unique_pins,
            "all_pins": all_pins
        }
    
    return sample_variants

def main(vartable_folder: Path, vcf_folder: Path, decodetable_path: Path, caller: str, results_folder: Path) -> None:

    logging.info("Hello - This is pin.py!")

    logging.info("Preparing input files")
    _, sample_pool_ids, _, _ = decode(decodetable_path, vartable_folder, caller)
    
    logging.info("Pinpointing variants to individuals")
    sample_variants = pin(vartable_folder, decodetable_path, caller)

    # Create a range of results files:
    # 1. Variants by chrom, pos, ref, alt in TSV format
    # 2. VCF files from each pool filtered by the pinnable sample variants
    # 3. Merged VCF files from each sample with combined information
    # 4. Merged VCF files from all samples with combined information without genotype/sample information
    # 5. Summary file with the number of unique and all pinnable variants per sample

    for v_type in ["unique_pins", "all_pins"]:
        logging.info(f"Processing {v_type.replace('_',' ')} by individual")
        output_folder = results_folder / v_type
        output_folder.mkdir(parents=True, exist_ok=True)
        
        for sample, pinnables in sample_variants.items():
            for pool in sample_pool_ids[sample]:
                filter_vcf(vcf_folder / f"{pool}.{caller}.vcf", pinnables[v_type], output_folder / f"{sample}_{pool}_{v_type}.vcf")

            merge_list = [output_folder / f"{sample}_{pool}_{v_type}.vcf" for pool in sample_pool_ids[sample]]
            with TemporaryDirectory() as tmpdir:
                process_vcf = VcfMerger(merge_list, f"{sample}_{v_type}", output_folder, Path(tmpdir))
                process_vcf.merge_individuals()
        
        logging.info(f"Merging {v_type.replace('_',' ')}")
        
        merge_list = [output_folder / f"{sample}_{v_type}.vcf.gz" for sample in sample_variants]
        with TemporaryDirectory() as tmpdir:
            process_vcf = VcfMerger(merge_list, f"{v_type}_merged", results_folder, Path(tmpdir))
            process_vcf.merge_matrix()

    logging.info("Summarizing results")
    with open(results_folder / "summary.tsv", 'w') as fout:
        for sample, pinnables in sample_variants.items():
            for v_type in ["unique_pins", "all_pins"]:
                print(sample, v_type, len(pinnables[v_type]), sep='\t', file=fout)
    
    logging.info("Creating lookup table")
    with open(results_folder / "lookup.tsv", 'w') as fout:
        for sample, pinnables in sample_variants.items():
            for pinnable in pinnables['all_pins']:
                shared = "shared_variant"
                if pinnable in pinnables['unique_pins']:
                    shared = "unique_variant"
                print(sample, pinnable, shared, sep='\t', file=fout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform pinning logic and variant filtering.')
    parser.add_argument('--vartable-folder', type=Path, default=Path("./"), help='Folder containing variant tables.')
    parser.add_argument('--vcf-folder', type=Path, default=Path("./"), help='Folder containing VCF files.')
    parser.add_argument('--decodetable', type=Path, default=Path("decodetable.tsv"), help='Path to the decodetable file.')
    parser.add_argument('--caller', type=str, default="GATK", help='Caller name.')
    parser.add_argument('--results-folder', type=Path, default=Path("results"), help='Folder to store the results.')
    
    args = parser.parse_args()
    
    main(args.vartable_folder, args.vcf_folder, args.decodetable, args.caller, args.results_folder)
