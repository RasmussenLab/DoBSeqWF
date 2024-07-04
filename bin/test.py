#!/usr/bin/env python
import logging
import argparse
from pathlib import Path

# 2024-02-27 mads
# Simple comparison of DoBSeq test output.

def test():
    """ Comparison of DoBSeq test output """
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Compare test data')
    parser.add_argument('-t', type=Path, help='True SNVs')
    parser.add_argument('-v', type=Path, help='Pipeline variants')
    parser.add_argument('-m', type=str, default='pinpy', help='Pinpoint method')
    args = parser.parse_args()
    truth_set = Path(args.t)
    pipeline_set = Path(args.v)
    pinpoint_method = args.m

    logger.info('Comparing test data')
    sample_vars_truth = dict()
    sample_vars = dict()

    # Get truth_set data:
    with open(truth_set) as fin:
        for line in fin:
            sample_info = line.strip().split()
            if sample_info[0] not in sample_vars_truth:
                sample_vars_truth[sample_info[0]] = set()
            sample_vars_truth[sample_info[0]].add(":".join(sample_info[1:]))
    
    # Get pipeline data:
    if pinpoint_method == 'pilot':
        with open(pipeline_set) as fin:
            var_info = fin.readline().strip().split()
            var_id_idx = var_info.index('var.ID')
            fam_id_idx = var_info.index('FAMID.A')
            for line in fin:
                sample_info = line.strip().split()
                if sample_info[fam_id_idx] not in sample_vars:
                    sample_vars[sample_info[fam_id_idx]] = set()
                sample_vars[sample_info[fam_id_idx]].add(sample_info[var_id_idx])
    elif pinpoint_method == 'new':
        with open(pipeline_set) as fin:
            for line in fin:
                sample_info = line.strip().split()
                if sample_info[0] not in sample_vars:
                    sample_vars[sample_info[0]] = set()
                sample_vars[sample_info[0]].add(sample_info[1])
    
    # Compare data:
    test_pass = True
    for sample in sample_vars_truth:
        logging.info(f"Comparing {sample}")
        if sample not in sample_vars:
            logging.info(f"{sample} not in pipeline output")
            test_pass = False
        elif sample_vars_truth[sample] == sample_vars[sample]:
            logging.info(f"{sample} matches")
        else:
            logging.info(f"{sample} does not match")
            logging.info(f"Pipeline: {sample_vars[sample]}")
            logging.info(f"Truth: {sample_vars_truth[sample]}")
            logging.info(f"Difference: {sample_vars[sample] - sample_vars_truth[sample]}")
            test_pass = False
    
    if test_pass:
        logging.info('Test passed')
        open("test.passed", 'w').close()
    else:
        logging.info('Test failed')
        open("test.failed", 'w').close()


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    project_dir = Path(__file__).resolve().parents[1]

    test()
