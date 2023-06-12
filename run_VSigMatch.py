'''
This script identifies viral composition of a metagenomics sample given a viral sequence database
using an adapted kmer-signature-based and similairty approach.
'''

import argparse
from src.process_inputs_outputs import *
from src.viral_chracterization import *

parser = argparse.ArgumentParser(description='Identification of viral composition of a metagenomics sample given a viral sequence database using an adapted kmer-signature-based and similairty approach.')
parser.add_argument('--METAfastq', metavar='METAfastq', type=str, default="inputs/sampled25_SRR23022001.fastq", help='fastq file containing metagenomics sequences to chracterize')
parser.add_argument('--kmer', metavar='k', type=int, default=5, help='length of substring')
parser.add_argument('--topN', metavar='N', type=int, default=20, help='number of top virus names prioritized by the method')
parser.add_argument('-output_dir', metavar='output_dir', type=str, help='name of directory to save output files')
args = parser.parse_args()

def main(METAfastq, k, N, output_dir):
    # get the processed database
    taxo_numkmers_sigs = get_CRVDB_sigs()
    taxonomy_names = get_taxonames(taxo_numkmers_sigs)
    # process the input metagenomics sample
    meta_reads = readin_fastq(METAfastq)
    ## convert all input reads to numerical kmers
    keynum_list = reads_to_numkeys(meta_reads, k)
    ## compute abundance matrix
    numpy_abundance_matrix = compute_abundance_matrix(keynum_list, k, taxo_numkmers_sigs)
    # compute jaccard vector
    jaccard_vector = compute_jaccard_vector(keynum_list, k, taxo_numkmers_sigs)
    # prioritize viral taxonomy names by the weighted k-mer-abundance similarity score
    weighted_ave_res = characterize_sample(jaccard_vector,numpy_abundance_matrix, taxonomy_names)
    # output the results: (i) bar plot of the top N viral taxonomy
    make_barplot(result_data = weighted_ave_res, top_num = N, output_dir=output_dir)
    # output the results: (ii) supplementary metadata of the viruses
    output_annotation(result_data = weighted_ave_res, top_num = N, output_dir=output_dir)
    print("---------------- DONE -----------------")

main(args.METAfastq, args.kmer, args.topN, args.output_dir)
