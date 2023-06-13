# VSigMatch: Virus Signature Matching Tool
## Objective
This program aims to characterize viral communities present in the sample using an adapted kmer-signature-based and similairty approach.

## Description
Given a metagenomics sample with short-read sequences in `.fastq` format, VSigMatch characterizes viral communities present in the sample using an adapted kmer-signature-based and similairty approach. First all the sequences in the sample are chopped up into substrings length of k (kmers). The set of sample kmers would then be compared against all the virus signatures in the processed database to get a vector of overall similarities. Based on kmer abundance in the sample and their presence in each virus signatures along with the similarity vector, VsigMatch score: weighted k-mer-abundance similarity is computed for every signature and used to prioritize the virus with high chances to be present in the sample.   

## Install dependencies
```
$ pip install -r requirements.txt
```
- [python ~= 3.11.3](https://www.python.org/downloads/release/python-3113/)
- [numpy~=1.24.3](https://numpy.org/install/)
- [pandas~=2.0.2](https://pandas.pydata.org/)
- [scipy~=1.10.1](https://scipy.org/install/)
- [matplotlib~=3.7.1](https://matplotlib.org/stable/users/installing/index.html)
- [seaborn~=0.12.2](https://seaborn.pydata.org/installing.html)

## Usage
```
Identification of viral composition of a metagenomics sample given a viral
sequence database using an adapted kmer-signature-based and similairty approach.

optional arguments:
  -h, --help            show this help message and exit
  --METAfastq METAfastq
                        fastq file containing metagenomics sequences to
                        chracterize
  --kmer k              length of substring
  --topN N              number of top virus names prioritized by the method
  -output_dir output_dir
                        name of directory to save output files
```

## Command line
To run this program in the command line interface type: 
```
$ python3 run_VSigMatch.py --METAfastq 'inputs/sampled25_SRR23022001.fastq' --kmer 5 --topN 20 -output_dir 'outputs'
```
or use below command to use all default parameters
```
$ python3 run_VSigMatch.py -output_dir 'outputs'
```
## Run unit tests
```
$ python3 tests/tests.py
```

## Input
- A [.fastq](https://www.zymoresearch.com/blogs/blog/fastq-file-format#:~:text=FASTQ%20format%20is%20a%20human,the%20sequencing%20platform%20flow%2Dcell.) file containing short-reads produced by metagenomics sequencing 

Example: `sampled25_SRR23022001.fastq`:
```
@SRR23022001.4.1 4 length=151
AAAGGTTTATACCTTCCCAGGTAACAAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGT
+SRR23022001.4.1 4 length=151
???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
@SRR23022001.10.1 10 length=151
AACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCT
+SRR23022001.10.1 10 length=151
???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
...
```
## Outputs
- **Output 1:** A Bar plot of the top 20 viral taxonomy names with highest VSigMatch scores

Example: `example_viral_char_barplot_k7.png`
![Viral charcterization of sampled25_SRR23022001.fastq](https://github.com/KewalinSamart/VSigMatch/blob/main/outputs/figures/example_viral_char_barplot_k7.png)

- **Output 2:** .txt file containing annotation of the top viruses: GenBank accession, Description, Taxonomy name

Example: `example_viral_annotation_k7.txt`
```
genbank_id	description	taxonomy_name
147	DQ136146.1	Tomato chlorosis virus isolate AT80/99 segment RNA2, complete sequence.	Tomato chlorosis virus
213	KX926428.1	Watermelon mosaic virus isolate WMV-Pg, complete genome.	Watermelon mosaic virus
301	KU754522.2	Lye Green virus isolate Dobs_PoolSeq4 VP1, VP2, VP3, putative glycoprotein, and putative polymerase genes, complete cds.	Lye Green virus
...
```

**Contact Info:** Kewalin Samart; kewalin.samart@cuanschutz.edu
