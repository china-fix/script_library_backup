'''this program is use to rename the contigs in fasta format, for the error meet during annotaion with prokka
[15:02:25] Contig ID must <= 37 chars long: NODE_9_length_100086_cov_12.2232_ID_17
[15:02:25] Please rename your contigs OR try '--centre X --compliant' to generate clean contig
script just remove all the old >names and add new name with NODE_#'''

import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord


def parse_args():
    parser=argparse.ArgumentParser(description="welcome to use Xiao_fei_Robot")
    parser.add_argument('--original_fasta', required=True, type=str, metavar='FILENAME', help='the original fasta file for renaming')
    parser.add_argument('--output', default="xiao_renamed_contigs", type=str, metavar='FILENAME', help='output file name')
    return parser.parse_args()

def main():
    args = parse_args()
    seq_records = SeqIO.parse(args.original_fasta, "fasta")
    n = 1
    new_seq_records =[]
    for seq_record in seq_records:
        seq_record.id = args.original_fasta.split('.',1)[0] + "_" + str(n)
        seq_record.name = '' #str(n)
        seq_record.description = ''#str(n)
        new_seq_records.append(seq_record)
        n = n + 1
    SeqIO.write(new_seq_records, args.output, "fasta")

if __name__ == '__main__':
    sys.exit(main())



