"doing BLAST online with the query fasta file of query sequences and save the blast output in xlm"

from Bio import SeqIO
from Bio.Blast import NCBIWWW
import sys
import argparse
from Bio import SeqRecord
from Bio import Seq
from Bio.Blast import NCBIXML

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS', required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blast")
    return parser.parse_args()

def main():
    args=parse_args()
    for seq_record in SeqIO.parse(args.CDS,"fasta"):
        result_handle = NCBIWWW.qblast("blastn","nt",seq_record.format("fasta"))

        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            # Do something with blast_record
            try:
                title = blast_record.alignments[0].title
                print(title + "---FIX---" + seq_record.description)
                #print(seq_record.description)
            except IndexError:
                print("no match"+"---FIX---"+ seq_record.description)
                #print(seq_record.description)

            #for alignment in blast_record.alignments:
                #title=alignment.title
                #print(title)
    



       # with open("my_blast.xml","a") as out_handle:
           # out_handle.write(result_handle.read())
 
if __name__ == '__main__':
    sys.exit(main())
