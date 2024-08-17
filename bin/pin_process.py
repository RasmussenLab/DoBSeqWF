#!/usr/bin/env python3

from pathlib import Path
from typing import List, Set, Dict, Tuple
from tempfile import TemporaryDirectory
from collections import namedtuple
import sys

## This is pin_process.py - script for processing pipeline output for easy downstream analysis.
# mads - 2024-08-17


Variant = namedtuple('Variant', ['varid','CHROM','POS','REF','ALT','is_snv','gene_target','coding_uncertain','gene_name','annotation','HGVSc','HGVSp','is_LOF','LOF_gene','CLNSIG','is_P','P_gene','is_LOFP','is_ACMG','ACMG_category','AD','DP_site','VAF'])

def get_variants(variant_table: Path, acmg: Tuple) -> Set[str]:
    """
    Takes a path to a variant table and returns a set of variants.
    Each variant is represented as a string formatted as "CHROM:POS:REF:ALT".
    
    Parameters:
    variant_table (Path): Path to a variant table.
    acmg (Tuple): A tuple containing ACMG gene information.
    
    Returns:
    set: A set containing variants.
    dict: A dictionary containing variant information.
    """
    with open(variant_table, 'r') as fin:
        header = fin.readline().strip().split('\t')
    
        if header[0:5] != ['CHROM', 'POS', 'ID', 'REF', 'ALT']:
            raise ValueError("Variant table is in unexpected format.")

        ANN_idx = header.index('ANN')
        GENE_idx = header.index('GENE')
        LOF_idx = header.index('LOF')
        CLNSIG_idx = header.index('CLNSIG')
        CLN_GENE_idx = header.index('GENEINFO')
        DP_site_idx = header.index('DP')
        for idx,element in enumerate(header):
            if '.AD' in element:
                AD_idx = idx
    
        # Construct a set of variants with the format "CHROM:POS:ALT"
        variants = set()
        variants_info = dict()
        for line in fin:
            var_info = line.strip().split('\t')
            
            # Basic info
            CHROM,POS,_,REF,ALT = var_info[0:5]
            varid = f"{CHROM}:{POS}:{REF}:{ALT}"
            is_snv = 1
            if len(ALT) > 1 or len(REF) > 1:
                is_snv = 0

            # Target region ID
            gene_target = var_info[GENE_idx]
            
            # Functional annotation:
            effects = var_info[ANN_idx].split(',')
            gene_name = ",".join([effect.strip().split('|')[3] for effect in effects])
            function = ",".join([effect.strip().split('|')[1] for effect in effects])
            HGVSc = ",".join([effect.strip().split('|')[9] for effect in effects])
            HGVSp = ",".join(['NA' if len(effect.strip().split('|')[10].strip()) == 0 else effect.strip().split('|')[10] for effect in effects])
            coding_uncertain = 0
            if len(effects) > 1:
                coding_uncertain = 1
            elif gene_name != gene_target:
                coding_uncertain = 1

            is_lof = 0
            lof_gene = 'NA'
            lof = var_info[LOF_idx]
            if lof != 'NA':
                is_lof = 1
                lof_gene = var_info[LOF_idx].split('|')[0].replace('(','')

            # ClinVar annotation:
            is_p = 0
            clnsig = var_info[CLNSIG_idx].strip()
            cln_gene = var_info[CLN_GENE_idx].split(':')[0].strip()
 
            if "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig:
                is_p = 1
            
            is_lofp = 0
            if is_p or is_lof:
                is_lofp = 1

            # ACMG state
            is_acmg = 0
            acmg_category = 'NA'
            for g in gene_name.split(','):
                if g in acmg:
                    is_acmg = 1
                    acmg_category = acmg[g][0]
            
            if lof_gene in acmg:
                is_acmg = 1
                acmg_category = acmg[lof_gene][0]

            if gene_target in acmg:
                is_acmg = 1
                acmg_category = acmg[gene_target][0]
            
            if cln_gene in acmg:
                is_acmg = 1
                acmg_category = acmg[cln_gene][0]
            
            DP_site = var_info[DP_site_idx].strip()
            AD = var_info[AD_idx]
            if AD.count(',') > 1:
                print(var_info)
                print('Multi-allelic sites are not split')
                sys.exit(1)

            AD = AD.split(',')[1].strip()
            assert CHROM.startswith('chr')
            assert POS.isdigit()
            assert REF.isalpha()
            assert ALT.isalpha()
            assert DP_site.isdigit()
            assert len(gene_target) > 0
            assert coding_uncertain in [0,1]
            assert len(gene_name) > 0
            assert len(function) > 0
            assert len(HGVSc) > 0
            assert len(HGVSp) > 0
            assert is_lof in [0,1]
            assert is_lofp in [0,1]
            assert len(lof_gene) > 0
            assert len(clnsig) > 0
            assert is_p in [0,1]
            assert len(cln_gene) > 0
            assert is_acmg in [0,1]
            assert len(acmg_category) > 0
            assert AD.isdigit()
            assert DP_site.isdigit()

            variants.add(varid)
            variants_info[varid] = Variant(
                varid,
                CHROM,
                POS,
                REF,
                ALT,
                is_snv,
                gene_target,
                coding_uncertain,
                gene_name,
                function,
                HGVSc,
                HGVSp,
                is_lof,
                lof_gene,
                clnsig,
                is_p,
                cln_gene,
                is_lofp,
                is_acmg,
                acmg_category,
                AD,
                DP_site,
                int(AD)/int(DP_site)
                )

    return variants, variants_info


def main():
    output = Path('output')
    vcftable = Path('data_links/vcftable.tsv')
    decode_table = Path('data_links/decodetable.tsv')
    pool_variant_tables = Path('data_links/variant_tables')
    pinpoint_variant_tables = Path('data_links/pinpoint_variant_tables')
    acmg_file = Path('data_links/acmg.tsv')
    acmg = {}
    samples = []
    decode = {}

    with open(acmg_file, 'r') as fin:
        fin.readline()
        for line in fin:
            try:
                info = line.split('\t')
                gene = info[0]
                disease = info[2]
                category = info[4]
            except IndexError:
                if len(line.strip()) == 0:
                    continue
                else:
                    print('error')
                    print(line)
                    sys.exit(1)
            acmg[gene] = (category,disease)

    with open(vcftable, 'r') as fin:
        for line in fin:
            sample_id = line.split()[0]
            samples.append(sample_id)
    
    pools = (set(),set())
    # row_pools = set()
    # column_pools = set()
    with open(decode_table, 'r') as fin:
        for line in fin:
            sample, pool1, pool2 = line.strip().split('\t')
            decode[sample] = (pool1,pool2)
            pools[0].add(pool1)
            pools[1].add(pool2)
    
    # Create combined pool variant table:
    with open(output / 'combined_pool_table.tsv', 'w') as fout:
        print('dim','pool','uvarid','varid','gene_target','CHROM','POS','REF','ALT','is_snv','is_acmg','is_lof','is_p','is_lofp','coding_uncertain','lof_gene','clnvar_gene','snpeff_gene','annotation','HGVSc','HGVSp','CLNSIG','acmg_category','AD','DP_site','VAF',file=fout,sep='\t')
        for n, dim in enumerate(['row','column']):
            for pool_id in pools[n]:
                filepath = pool_variant_tables / (f"{pool_id}.GATK.tsv")
                pool_variants, pool_info = get_variants(filepath,acmg)
                for variant in pool_variants:
                    print(
                        dim,
                        pool_id,
                        f"{pool_id}:{pool_info[variant].varid}",
                        pool_info[variant].varid,
                        pool_info[variant].gene_target,
                        pool_info[variant].CHROM,
                        pool_info[variant].POS,
                        pool_info[variant].REF,
                        pool_info[variant].ALT,
                        pool_info[variant].is_snv,
                        pool_info[variant].is_ACMG,
                        pool_info[variant].is_LOF,
                        pool_info[variant].is_P,
                        pool_info[variant].is_LOFP,
                        pool_info[variant].coding_uncertain,
                        pool_info[variant].LOF_gene,
                        pool_info[variant].P_gene,
                        pool_info[variant].gene_name,
                        pool_info[variant].annotation,
                        pool_info[variant].HGVSc,
                        pool_info[variant].HGVSp,
                        pool_info[variant].CLNSIG,
                        pool_info[variant].ACMG_category,
                        pool_info[variant].AD,
                        pool_info[variant].DP_site,
                        round(pool_info[variant].VAF, 6),
                        file=fout,
                        sep='\t')

    # Create combined pinpoint variant table:
    with open(output / 'combined_pinpoint_table.tsv', 'w') as fout:
        print('sample_id','gene_target','varid','CHROM','POS','REF','ALT','is_snv','is_acmg','is_lof','is_p','is_lofp','coding_uncertain','lof_gene','clnvar_gene','snpeff_gene','annotation','HGVSc','HGVSp','CLNSIG','acmg_category','column','AD_column','DP_site_column','VAF_column','row','AD_row','DP_site_row','VAF_row',file=fout,sep='\t')
        
        for sample_id in decode:
            col = decode[sample_id][0]
            filepath = pinpoint_variant_tables / (f"{sample_id}_{col}_all_pins.tsv")
            col_variants, col_info = get_variants(filepath,acmg)
            
            row = decode[sample_id][1]
            filepath = pinpoint_variant_tables / (f"{sample_id}_{row}_all_pins.tsv")
            row_variants, row_info = get_variants(filepath,acmg)
            
            assert col_variants == row_variants

            for variant in col_variants:
                assert col_info[variant].CHROM == row_info[variant].CHROM
                assert col_info[variant].POS == row_info[variant].POS
                assert col_info[variant].REF == row_info[variant].REF
                assert col_info[variant].ALT == row_info[variant].ALT
                assert col_info[variant].varid == row_info[variant].varid
                assert col_info[variant].gene_target == row_info[variant].gene_target
                assert col_info[variant].gene_name == row_info[variant].gene_name
                assert col_info[variant].coding_uncertain == row_info[variant].coding_uncertain
                assert col_info[variant].annotation == row_info[variant].annotation
                assert col_info[variant].HGVSc == row_info[variant].HGVSc
                assert col_info[variant].HGVSp == row_info[variant].HGVSp
                assert col_info[variant].is_LOF == row_info[variant].is_LOF
                assert col_info[variant].LOF_gene == row_info[variant].LOF_gene
                assert col_info[variant].CLNSIG == row_info[variant].CLNSIG
                assert col_info[variant].is_P == row_info[variant].is_P
                assert col_info[variant].P_gene == row_info[variant].P_gene
                assert col_info[variant].is_LOFP == row_info[variant].is_LOFP
                assert col_info[variant].is_ACMG == row_info[variant].is_ACMG
                assert col_info[variant].ACMG_category == row_info[variant].ACMG_category

                print(
                    sample_id,
                    col_info[variant].gene_target,
                    col_info[variant].varid,
                    col_info[variant].CHROM,
                    col_info[variant].POS,
                    col_info[variant].REF,
                    col_info[variant].ALT,
                    col_info[variant].is_snv,
                    col_info[variant].is_ACMG,
                    col_info[variant].is_LOF,
                    col_info[variant].is_P,
                    col_info[variant].is_LOFP,
                    col_info[variant].coding_uncertain,
                    col_info[variant].LOF_gene,
                    col_info[variant].P_gene,
                    col_info[variant].gene_name,
                    col_info[variant].annotation,
                    col_info[variant].HGVSc,
                    col_info[variant].HGVSp,
                    col_info[variant].CLNSIG,
                    col_info[variant].ACMG_category,
                    col,
                    col_info[variant].AD,
                    col_info[variant].DP_site,
                    round(col_info[variant].VAF, 6),
                    row,
                    row_info[variant].AD,
                    row_info[variant].DP_site,
                    round(row_info[variant].VAF, 6),
                    file=fout,
                    sep='\t'
                )


if __name__ == '__main__':
    main()