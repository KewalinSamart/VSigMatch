'''
Functions to compute kmer abundance matrix 

Author: Kewalin Samart
Created date: 06/10/2023
Last Modified: 06/11/2023
'''

import pandas as pd
import numpy as np

def get_kmers(input_seq, k):
    '''
    function to generate k-mers i.e., set of strings length of k
    @param input_seq: an input sequence from either the sample or database
    @param k: an integer indicating length of k-mer to create
    @returns kmer_list: a list of generated k-mers
    '''
    if k > len(input_seq):
        print("kmer must be shorter than the sequence length")
        
    num_kmers = len(input_seq) - k + 1
    # loop over the kmer start positions
    kmer_list = list()
    for i in range(num_kmers):
        # slice the string to get the kmer
        kmer = input_seq[i:i+k]
        kmer_list.append(kmer)

    return kmer_list

def get_numerical_rep(kmer):
    '''
    function to compute numerical representations of kmers
    @param kmer: string of nucleotides lenth of k 
    @returns key_val: kmer's numerical representation 
    '''
    kmer_char = list(reversed([*kmer]))
    kmer_len = len(kmer)
    key_val = 0
    char_val = 0
    if len(kmer) != 0:
        for i in reversed(range(kmer_len)):
            char = kmer_char[i]
            if char == "A":
                char_val = 0
            elif char == "C":
                char_val = 1
            elif char == "G":
                char_val = 2
            elif char == "T":
                char_val = 3
            key_val = key_val + char_val*(4**i) 
    else:
        key_val = None

    return key_val

def jaccard(A, B):
    '''
    function to compute jaccard similarity between two set of kmers
    @param A: a set of kmers
    @param B: a set of kmers
    @returns jaccard similarity of the two input sets
    '''
    intersection = len(list(set(A).intersection(B)))
    union = (len(A) + len(B)) - intersection

    return float(intersection) / union

def reads_to_numkeys(meta_reads, k):
    '''
    function to convert all input reads to numerical kmers
    @param mata_read: a list of reads in the given sample
    @param k: length of kmers 
    @returns keynum_list: a list of numerical kmers
    '''
    keynum_list = []
    for read in meta_reads:
        kmers = get_kmers(read, k) 
        keynums = [get_numerical_rep(kmer) for kmer in  kmers]
        keynum_list.extend(keynums)
    return keynum_list

def get_kmer_abundance(keynum_list):
    '''
    function to count unique kmers present in the sample.
    @param keynum_list: a list of numerical kmers
    @returns abundance_dict: a dictionary with keys: unique kmers and values: kmer counts
    '''
    freq_keynums = pd.Series(keynum_list).value_counts()
    abundance_dict = freq_keynums.to_dict()
    return abundance_dict

def get_taxonames(taxo_numkmers_sigs):
    '''
    function to get viral taxonomy names in the database
    @param taxo_numkmers_sigs: kmer signature table for all the viral taxonomy names in the database
    @return taxonomy_names: a list of viral taxonomy names
    '''
    taxonomy_names = list(taxo_numkmers_sigs['taxonomy_name'])
    return taxonomy_names

def compute_abundance_matrix(keynum_list, k, taxo_numkmers_sigs):
    '''
    function to compute abundance matrix
    @param keynum_list: a list of numerical kmers
    @param k: length of kmers
    @param taxo_numkmers_sigs: kmer signature table for all the viral taxonomy names in the database
    '''
    print("computing abundance matrix...")
    abundance_dict = get_kmer_abundance(keynum_list)
    unique_keynums = set(keynum_list)
    abundance_matrix = pd.DataFrame(list(unique_keynums), columns=["sample_kmer"])
    signatures = list(taxo_numkmers_sigs['kmers_sig_{}'.format(k)])
    taxonomy_names = get_taxonames(taxo_numkmers_sigs)

    for taxo_index in range(len(taxonomy_names)):
        abundance_sig = []
        signature = signatures[taxo_index]
        if isinstance(signature, str):
            signature = signature.split(",")
            num_signature = [eval(ele) for ele in signature]
            for keynum in unique_keynums:
                if keynum in num_signature:
                    count = abundance_dict[keynum]
                else:
                    count = 0
                abundance_sig.append(count)
        else:
            abundance_sig = [0] * len(unique_keynums)
        #abundance_matrix[taxonomy_names[taxo_index]] = abundance_sig
        abundance_matrix = pd.concat([abundance_matrix, pd.Series(abundance_sig)],axis = 1)
    abundance_matrix = abundance_matrix.drop('sample_kmer',axis=1)
    numpy_abundance_matrix = abundance_matrix.to_numpy()

    return numpy_abundance_matrix

def compute_jaccard_vector(keynum_list, k, taxo_numkmers_sigs):
    '''
    function to compare similarity between the matagenomics sample and all the viral signatures in the database
    @param keynum_list: a list of numerical kmers
    @param k: length of kmers
    @param taxo_numkmers_sigs: viral-taxonomy signature table (pandas dataframe)
    @returns jc_vector: jaccard vector (numpy array format)
    '''
    print("computing jaccard vector: similairty between the sample and viral sequences in database")
    jc_vector = []
    signatures = list(taxo_numkmers_sigs['kmers_sig_{}'.format(k)])
    for signature in signatures:
        if isinstance(signature, str):
            signature = signature.split(",")
            signature = [int(ele) for ele in signature]
            jc_score = jaccard(signature, keynum_list)
        else:
            jc_score = 0
        jc_vector.append(jc_score)
    return np.array(jc_vector)

def characterize_sample(jaccard_vector,numpy_abundance_matrix, taxonomy_names):
    '''
    function to calculate the weighted k-mer-abundance similarity scores for all the taxonomy
    @param jaccard_vector: (numpy array format)
    @param numpy_abundance_matrix: abundance matrix where rows: unique sample numerical kmers
    @param taxonomy_names: a list of viral taxonomy names
    '''
    print("computing the final result matrix and summarize the characterization result")
    combined_mat = jaccard_vector * numpy_abundance_matrix 
    normalizedData = (combined_mat-np.min(combined_mat))/(np.max(combined_mat)-np.min(combined_mat))
    combined_res = pd.DataFrame(normalizedData, columns = taxonomy_names)
    weighted_ave_res = combined_res.mean() 
    return  weighted_ave_res
