'''
Functions to process inputs: 
(i) a fastq file containing short-read sequences in a metagenomics sample
(ii) a fasta file containing viral sequences obtained from Reference Viral DataBase (RVDB; lasted clustered DB 04/2023)

Author: Kewalin Samart
Created date: 06/08/2023
Last Modified: 06/11/2023
'''

import sys, os
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 

def readin_fastq(input_fastq):
    '''
    function to read in metagenomics sequences from an input fastq file
    @param input_fastq: path to an input fastq file (metagenomics sample)
    @returns meta_reads: a list of metagenomics reads 
    '''
    meta_reads = []
    nucleotides = ['A', 'T', 'C', 'G']

    with open(input_fastq, 'r') as fastq:
        lines = fastq.readlines()

    for line in lines:
        line = line.strip()
        if line.startswith(tuple(nucleotides)):
            meta_reads.append(line)
        else:
            next
    return meta_reads

def get_dbseqs(input_fasta):
    '''
    function to read in viral sequences from an input database fasta file
    @param input_fasta: path to an input database fasta file
    @returns viral_seqs: a dict with keys: tuple(genbank_id, description, taxonomy_name, category),
                         values: viral sequences
    '''
    viral_seqs = {}
    nucleotides = ['A', 'T', 'C', 'G', 'N']

    with open(input_fasta, 'r') as fasta:
        lines = fasta.readlines()

    whole_seq = "" 
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith(">"): 
            print(i+1)
            if i == 0:
                key = tuple(line.split(sep="|")[2:6]) # genbank_id, description, taxonomy_name, category
            else:
                if len(key) == 4:
                    viral_seqs[key] =  whole_seq.replace("N","")
                    key = tuple(line.split(sep="|")[2:6])
                    whole_seq = ""
        if line.startswith(tuple(nucleotides)):
            if whole_seq == "": 
                whole_seq = line
            else:
                whole_seq = whole_seq + line
    if len(key) == 4:
        viral_seqs[key] =  whole_seq.replace("N","") 

    return viral_seqs

def output_metaviral(viral_seqs, output_dir):
    '''
    function to output txt format of the viral_seqs derived from the database fasta file.
    @param viral_seqs: a dict with keys: tuple(genbank_id, description, taxonomy_name, category),
                       values: viral sequences
    @param output_dir: a path to a directory to store the output
    '''
    if isinstance(output_dir, str):
        output_file = open(output_dir + "DB_metadata.txt","w")
        output_file.write("genbank_id" + "\t" + "description" + 
                          "\t" + "taxonomy_name" + "\t" + "category" +
                            "\t" + "sequence" + "\n")
        meta_list = list(viral_seqs.keys())
        i = 0
        for meta in meta_list:
            if len(meta) == 4:
                output_file.write(str(meta[0]) + "\t" + str(meta[1]) + "\t" + 
                                str(meta[2]) + "\t" + str(meta[3]) + "\t" + 
                                str(viral_seqs[meta]) + "\n")   
            i = i + 1
    else:
        sys.exit("output_dir must be type of string")
    output_file.close()

def get_CRVDB_sigs():
    '''
    function to get the processed reference virus database (with viral signatures)
    '''
    dir = os.getcwd()
    taxo_numkmers_sigs = pd.read_csv("{}/database/example_processedDB_local.txt".format(dir), sep="\t")
    
    return taxo_numkmers_sigs

def make_barplot(result_data, output_dir, top_num = 20):
    '''
    function to output a bar plot of the sample's viral characterization
    @param result_data: a pandas dataframe with row names: viral taxonomy names, 
                        and the first column: score indicating kmer-signature-based similarity
    @param top_num: an integer indicating number of the top viral taxonomy names to plot
    '''
    sorted_res = result_data.sort_values(ascending=False)[1:top_num]
    sorted_res_df = sorted_res.to_frame()
    plt.figure(figsize=(10,6))
    ax = sns.barplot(x = sorted_res_df[0],y = sorted_res_df.index,
                    data = sorted_res_df, estimator = np.median, ci = 0, orient='h')
    plt.title('Viral characterization of the sample')
    ax.set(xlabel="score", ylabel="viral taxonomy")
    plt.savefig('{}/viral_char_barplot.png'.format(output_dir))
    plt.show()

def output_annotation(result_data, output_dir, top_num = 20):
    '''
    function to output the evolutionay information of top viral taxonomy prioritized by the method
    @param result_data: a pandas dataframe with row names: viral taxonomy names, 
                        and the first column: score indicating kmer-signature-based similarity
    @param top_num: an integer indicating number of the top viral taxonomy names to output the evolutionay information
    '''
    sorted_res = result_data.sort_values(ascending=False)[1:top_num]
    sorted_res_df = sorted_res.to_frame()
    top_taxonomy = list(sorted_res_df.index)
    dir = os.getcwd()
    metadf = pd.read_csv('{}/metadata/viralGenBank_meta.txt'.format(dir), sep='\t')
    meta_subset = metadf[metadf['taxonomy_name'].isin(top_taxonomy)]
    meta_subset =  meta_subset.drop('Unnamed: 0', axis=1)
    meta_subset.to_csv('{}/viral_annotation.txt'.format(output_dir), sep='\t')
