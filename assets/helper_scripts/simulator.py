import logging
import sys
import subprocess
import numpy as np
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from Bio.SeqIO.FastaIO import SimpleFastaParser

np.random.seed(12)

# 2024-02-15 mads
# Simulate all necessary input files for a DoBSeq matrix data analysis.
# Only single chromosome.

def simulate():
    """ Simulate single chromosome variation in a X * X dobseq matrix."""
    logger = logging.getLogger(__name__)
    
    id = "2x2"

    reference = project_dir / "assets/data/reference_genomes/small/small_reference.fna"
    output = project_dir / f"assets/data/simulated_data/sim_{id}"
    regions = [[1, 710], [720, 1300]]

    msize = 2       # Matrix dimension size (msize x msize)
    snv = 1         # Single nucleotide variation per individual
    cov = 35        # Coverage per allele
    err = 0         # Base error rate
    mut = 0         # Mutation fraction
    indel = 0       # Indel fraction
    ext = 0         # Indel extension

    logger.info('running simulation')

    with open(reference) as fin:
        header, ref_seq = next(SimpleFastaParser(fin))
        ref_seq_len = len(ref_seq)

    n_reads = int(ref_seq_len * cov / 150)

    # Create dobseq matrix:
    matrix = np.arange(0, msize**2).reshape(msize, msize)
    
    # Create mutated sequences:
    sequence_folder = output / "sequence_variations"
    sequence_folder.mkdir(parents=True, exist_ok=True)

    sample_folder = output / "pools"
    sample_folder.mkdir(parents=True, exist_ok=True)

    snvlist = output / "snvlist.tsv"
    with open(snvlist, 'w') as flist:
        for seed_n, i in enumerate(range(msize**2)):
            # Create sequence variations:
            with open(sequence_folder / f"sim_{i}.fna", "w") as fout:
                snvs = dict()
                for j in range(snv):
                    # Choose random region from regions:
                    region = regions[np.random.randint(0, len(regions))]
                    # Choose random index in region:
                    snv_idx = np.random.randint(region[0], region[1])
                    ref = ref_seq[snv_idx]
                    choices = ["A", "T", "C", "G"]
                    choices.remove(ref)
                    alt = np.random.choice(choices)
                    snvs[snv_idx] = (ref, alt)
                # Create mutated sequence:
                mut_seq = list(ref_seq)
                for idx, (ref, alt) in snvs.items():
                    mut_seq[idx] = alt
                    flist.write(f"sim_{i}\t{idx+1}\t{ref}\t{alt}\n")
                # Write mutated sequence to file:
                fout.write(f">{header}\n")
                fout.write("".join(mut_seq))
                fout.write("\n")
            
            # Create fastq files:
            # Random seed is changed to avoid duplicate reads.
            cmd = ["wgsim", "-S", str(seed_n+3), "-N", str(n_reads), "-1", "150", "-2", "150", "-e", str(err), "-r", str(mut), "-R", str(indel), "-X",
                   str(ext), sequence_folder / f"sim_{i}.fna", sample_folder / f"sim_{i}_1.fq", sample_folder / f"sim_{i}_2.fq"]
            return_value = subprocess.run(cmd, stdout=subprocess.DEVNULL)
            if return_value.returncode != 0:
                logging.error(f"Error running wgsim")
                sys.exit()
    
    pool_folder = output / "pools"
    pool_folder.mkdir(parents=True, exist_ok=True)

    # Combine into pools horizontally:
    for i in range(msize):
        samples = matrix[i, :]
        # Combine reads:
        for R in ["1", "2"]:
            files = [sample_folder / f"sim_{j}_{R}.fq" for j in samples]
            with open(pool_folder / f"B0_H{i}_{R}.fq", "w") as fout:
                for file in files:
                    with open(file) as fin:
                        for line in fin:
                            fout.write(line)
        
    # Combine pools vertically:
    for i in range(msize):
        samples = matrix[:, i]
        # Combine reads:
        for R in ["1", "2"]:
            files = [sample_folder / f"sim_{j}_{R}.fq" for j in samples]
            with open(pool_folder / f"B0_V{i}_{R}.fq", "w") as fout:
                for file in files:
                    with open(file) as fin:
                        for line in fin:
                            fout.write(line)

    # Gzip fastq files:
    
    cmd = ["gzip", "-f"] + list(pool_folder.glob("*.fq"))
    return_value = subprocess.run(cmd, stdout=subprocess.DEVNULL)
    if return_value.returncode != 0:
        logging.error(f"Error running gzip")
        sys.exit()

    # Create bed file with regions:
    bed_file = output / "target_calling.bed"
    with open(bed_file, "w") as fout:
        for region in regions:
            fout.write(f"{header}\t{region[0]}\t{region[1]}\n")
    
    # Create sampletable with pool ID, fastq1, fastq2:
    sampletable = output / "pooltable.tsv"
    with open(sampletable, "w") as fout:
        for i in ["H", "V"]:
            for j in range(msize):
                fout.write(f"B0_{i}{j}\t{pool_folder / f'B0_{i}{j}_1.fq.gz'}\t{pool_folder / f'B0_{i}{j}_2.fq.gz'}\n")
    
    # Crate decode table with sample ID, pool ID1, pool ID2:
    decodetable = output / "decodetable.tsv"
    with open(decodetable, "w") as fout:
        for i in range(msize):
            for j in range(msize):
                fout.write(f"sim_{i*msize+j}\tB0_H{i}\tB0_V{j}\n")

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    project_dir = Path(__file__).resolve().parents[2]
    load_dotenv(find_dotenv())

    simulate()
