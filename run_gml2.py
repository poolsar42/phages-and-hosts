#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess, os, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


# In[3]:


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', required=True, help='output prefix')
parser.add_argument('--gms2_path', default='/home/emil/Desktop/Sci/gms/gms2.pl', help='location of prodigal executable (default=PATH)')
#parser.add_argument('-c', '--code', default=11, type=int, help='genetic code to search with')

args= parser.parse_args()


# In[4]:


def build_seq_dict(fasta_file):
    '''Parse a fasta file and store the sequences in a dict. Header is the key to the dict'''
    new_dict = {}
    with open(fasta_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            new_dict[id.split(' ')[0]] = seq
    return new_dict


# In[7]:


def run_prodigal_w_specific_code(fasta_test, output_prefix):
    '''Run prodigal on the selected fasta with the selected genetic code, parse the results and store them in a table'''
    #Run prodigal, throw output into /dev/null because we don't need it and it clogs the termina
    with open(os.devnull) as out_handle:
        subprocess.call([args.gms2_path, '--seq', fasta_test, '--genome-type', 'auto'
                            '--output', output_prefix + '_' + '.faa'], stdout=out_handle, stderr=out_handle)

    nucl_seq_dict = build_seq_dict(fasta_test)

    #initialize table
    out_handle = open(output_prefix + ' results', 'w')
    out_handle.write('\t'.join(['contig_id', 'prot_id', 'start', 'end', 'length(bp)', 'strand', 'start_codon', 'end_codon', 'partial',
                                'rbs_motif', 'rbs_spacer', 'gc_content', 'aa_seq']) + '\n')

    #parse header information and store info in table
   # with open('.'.join(fasta_test.split('.')[-1]) + '_' + str(code) + '.faa') as in_handle:
    #    for id, seq in SimpleFastaParser(in_handle):
    #        splitId = id.split('#')
    #        prot_id = splitId[0].strip(' ')
    #        seq_id = '_'.join(prot_id.split('_')[:-1])
    #        start = int(splitId[1].strip(' ')) - 1
    #        end = int(splitId[2].strip(' '))
    #        length = end - start
    #        strand = int(splitId[3].strip(' '))
    #        add_info = splitId[-1].split(';')
    #        partial = True if '1' in add_info[1] else False
    #        start_codon = 'None'
    #        end_codon = 'None'
    #        if strand == 1:
    #            start_codon =  nucl_seq_dict[seq_id][start:start + 3]
    #            end_codon = nucl_seq_dict[seq_id][end -3:end]
    #        else:
    #            start_codon = Seq(nucl_seq_dict[seq_id][end-3 :end], generic_dna).reverse_complement()
    #            end_codon = Seq(nucl_seq_dict[seq_id][start :start + 3], generic_dna).reverse_complement()

    #        out_handle.write('\t'.join(list(map(str, [seq_id, prot_id, start+1, end, length, strand, start_codon, end_codon, partial,
     #                                                 add_info[3].split('=')[-1], add_info[4].split('=')[-1], add_info[5].split('=')[-1],
     #                                                 seq]))) + '\n')
   # out_handle.close()

