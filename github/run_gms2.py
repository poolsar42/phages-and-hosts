import subprocess, os, argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', required=True, help='output prefix')
parser.add_argument('--gms2_path', default='gms2.pl', help='location of prodigal executable (default=PATH)')

args= parser.parse_args()



def run_prodigal_w_specific_code(fasta_test, output_prefix):
    with open(os.devnull) as out_handle:
        system_requiry = ['perl', args.gms2_path, '--seq', fasta_test + '.fasta', '--genome-type', 'auto', '--output', output_prefix + '.OUT']
        subprocess.call(system_requiry,stdout=out_handle, stderr=out_handle)

run_prodigal_w_specific_code(args.input, args.output)


