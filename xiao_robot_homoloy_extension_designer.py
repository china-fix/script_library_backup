"""this is a robot used to design homology extensions (H1 and H2) primer part. can be use in lamda red mutation system.
written by Xiao Fei china-fixing@hotmail.com"""

import sys
import argparse
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--reference', required=True, type=str, metavar='FILENAME', help="they fasta reference (1 contig) sequence you want to design with")
    parser.add_argument('--del_start', required=True, type=int, metavar='delete_mutation_start', help="the start position (include the position i) you want to delete in the reference")
    parser.add_argument('--del_end', required=True, type=int, metavar='delete_mutation_end', help="the end position (include the position itself) you want to delete in the reference")
    parser.add_argument('--homo_length', default=40, type=int, metavar='homology_extension_length', help="the homology extension length you want, default is 40")
    #parser.add_argument('--output', default='homo_sequences', type=str, metavar='FILENAME', help="Output filename, default name is homo_sequences")
    return parser.parse_args()

def main():
    args = parse_args()
    seq_record = SeqIO.read(args.reference, "fasta")
    seq_upstream = seq_record [args.del_start -1 - args.homo_length : args.del_start - 1]
    seq_upstream.id = seq_upstream.id+"_up_"+str(args.del_start)
    seq_downstream = seq_record[args.del_end : args.del_end + args.homo_length]
    seq_downstream.id = seq_downstream.id+"_down_"+str(args.del_start)
    seq_downstream_rc = seq_downstream
    seq_downstream_rc.id = seq_downstream_rc.id + "_rc"
    seq_downstream_rc.seq = seq_downstream_rc.seq.reverse_complement()
    print(seq_upstream)
    print(seq_downstream_rc)
    print("----------------------------------------------------------------------")
    print("primer_F_homo_extension(5'->3'): %s" % str(seq_upstream.seq))
    print("primer_R_homo_extension(5'->3'): %s" % str(seq_downstream_rc.seq))
    print("you can just copy or you can check the information in the output file, enjoy!")
# about the output, if needed you can add by yourself!


if __name__ == '__main__':
    sys.exit(main())