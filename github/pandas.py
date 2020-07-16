#!/usr/bin/env python
# coding: utf-8

# In[4]:


import subprocess, os, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', required=True, help='output prefix')
parser.add_argument('-c', '--code', default=11, type=int, help='genetic code to search with')
args= parser.parse_args()

def build_seq_dict(fasta_file):
    '''Parse a fasta file and store the sequences in a dict. Header is the key to the dict'''
    new_dict = {}
    with open(fasta_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            new_dict[id.split(' ')[0]] = seq
    return new_dict

nucl_seq_dict = build_seq_dict(args.input)
with open(args.output, 'w') as f:
    for key,val in nucl_seq_dict.items():
        f.write('{}:{}\n'.format(key,val))
f.close()