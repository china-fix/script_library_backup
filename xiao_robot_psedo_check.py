import sys
import argparse
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

"""This robot is used to check the CDS fasta file psedo or not
Xiao 13-10-2019"""

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--list', required=True, type=str, metavar='FILENAME', help="gene name list that your want to check")
    parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="sample names of group 1, each line one name")
    parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="sample names of group 2, each line one name")
    parser.add_argument('--output', required=True, type=str, metavar='FILENAME', help="Output filename")
    return parser.parse_args()

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (bacteria table)."""
    return SeqRecord (seq = nuc_record.seq.translate(cds=True, table= "Bacterial"), \
                     id = nuc_record.id, \
                     description = "translation of CDS, using bacteria table")

def main():
    args=parse_args()
    with open(args.list) as key_list:
        proteins=[]
        for key in key_list.read().splitlines():
            for nuc_record in SeqIO.parse("temp_" + key, "fasta"):
                try:
                    protein = make_protein_record(nuc_record)
                    proteins.append(protein)
                except:
                    pass
                    #print(nuc_record.id)
            SeqIO.write(proteins, key + ".p", "fasta")
            proteins =[]   
            try:
                proteins_dict = SeqIO.index(key + ".p", "fasta")
            except ValueError:
                #pass
                proteins_dict = {record.id:record for record in SeqIO.parse(key + ".p",'fasta')}
                #d_seq  ={'%s%s%i'%(record.id.strip(),separator,i):record for i,record in enumerate(SeqIO.parse(key + ".p",'fasta'))} 
                #SeqIO.write([SeqRecord(record.seq ,id=record.id,name=record.name,description=record.description) for key_1,record in d_seq.items()],'%s_clean'%(key + ".p"),'fasta')
                #proteins_dict = SeqIO.index(key + ".p", "fasta")
            #print(list(proteins_dict.keys()))
            with open(args.GROUP_1) as group_1:
                group_1_positive_sum = 0
                group_1_sum = 0
                for name_1 in group_1.read().splitlines():
                    group_1_positive_sum += 1 
                    if name_1 +'---FIX---' + key in list(proteins_dict.keys()):
                        group_1_sum += 1
                    #else:
                        #print(name_1 +'---FIX---' + key)
            with open(args.GROUP_2) as group_2:
                group_2_positive_sum = 0
                group_2_sum = 0
                for name_2 in group_2.read().splitlines():
                    group_2_positive_sum += 1 
                    if name_2 +'---FIX---' + key in proteins_dict.keys():
                        group_2_sum += 1
                    #else:
                     #   print(name_2 +'---FIX---' + key)
            print(key + "---FIX---" + str(group_1_sum) + "---FIX---" + str(group_1_positive_sum) + "---FIX---" + str(group_2_sum) + "---FIX---" + str(group_2_positive_sum), file=open(args.output, "a"))

  
if __name__ == '__main__':
    sys.exit(main())
