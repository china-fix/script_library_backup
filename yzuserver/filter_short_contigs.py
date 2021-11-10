'''this program is help to filter and remove the short contigs'''

import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord


def parse_args():
    parser=argparse.ArgumentParser(description="welcome to use Xiao_fei_Robot")
    parser.add_argument('--original_fasta', required=True, type=str, metavar='INT', help='the original fasta file for filtering')
    parser.add_argument('--length', default=200, type=int, metavar='FILENAME', help='contigs whose length short than this num will be drop, default is 200')
    parser.add_argument('--output', default="xiao_filtered_contigs", type=str, metavar='FILENAME', help='output file name')
    return parser.parse_args()

def main():
    args = parse_args()
    seq_records = SeqIO.parse(args.original_fasta, "fasta")
    new_seq_records =[]
    drop_seq_records =[]
    for seq_record in seq_records:
        if len(seq_record.seq) > args.length:
            new_seq_records.append(seq_record)
        else:
            drop_seq_records.append(seq_record)
    SeqIO.write(new_seq_records, args.output, "fasta")
    SeqIO.write(drop_seq_records, args.output + ".drop", "fasta")

if __name__ == '__main__':
    sys.exit(main())



