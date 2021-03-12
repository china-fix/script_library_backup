'''
This robot can help to locally blastx the CDS.fasta to the reference epitops and filter the orginal CDS.fasta to matched_CDS.fasta and un_matched_CDS.fasta  
'''

import subprocess
import sys
import argparse
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import copy


def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS', required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blastx")
    parser.add_argument('--REF', required=True, type=str, metavar='FILENAME', help="the complete reference epitopes fasta format you want to blast to")
    parser.add_argument('--CUTOFF', default=0.9, type=float, metavar='DEFAULT 0.9', help="the lowest similarity value which classify as matched")
    parser.add_argument('--CUTOFF_1', default= 2, type=float, metavar='DEFAULT 2', help="the max value of the unmatched amino acids")
    parser.add_argument('--RULE_1', action='store_const', const=True, metavar='use the CUTOFF_1 rule', help="optional, use the unmatched amino acids number as cutoff")
    parser.add_argument('--OUT', default="xiao_robot_match_classify_CDS", type=str, metavar='outputname', help="Output name")
    parser.add_argument('--ANN', action='store_const', const=True, metavar='use the annotation funtion (pluse with --ANN_FILE together)', help="optional, use the annotation function to output a table with description of epitopes")
    parser.add_argument('--ANN_FILE', default="epitope_table_Gammaproteobacteria.csv", type=str, metavar='epitope annotation file name', help="epitope annotation file name")
    return parser.parse_args()

#1. blastx and get the xml
#2. parse the xml and get the matched and un_matched list
#3. extract and output the matched.fasta and un_matched.fasta

def doing_blast(query, reference):
    subprocess.run(["blastx", "-query", query, "-subject", reference, "-outfmt", "5", "-out", "temp_blast.xml"], check=True)


def filter_matching(cutoff):
    matched_list=[]
    temp_blast_handle = open("temp_blast.xml")
    blast_records = NCBIXML.parse(temp_blast_handle)
    for blast_record in blast_records:
        # Do something with blast_record
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            if  hsp.identities / alignment.length >= cutoff : # and hsp.expect <= 1e-10 : #and hsp.query_start == 1 and hsp.query_end == blast_record.query_letters:
                matched_result=blast_record.query
                matched_list.append(matched_result)

                print(alignment.title.split()[0]+'\t'+blast_record.query
                        , file=open("temp_cross_result.csv_1", "a"))
            else:
                pass
    return matched_list

def filter_matching_1(cutoff):
    matched_list=[]
    temp_blast_handle = open("temp_blast.xml")
    blast_records = NCBIXML.parse(temp_blast_handle)
    for blast_record in blast_records:
        # Do something with blast_record
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            if  alignment.length - hsp.identities <= cutoff : # and hsp.expect <= 1e-10 : #and hsp.query_start == 1 and hsp.query_end == blast_record.query_letters:
                matched_result=blast_record.query
                matched_list.append(matched_result)

                print(alignment.title.split()[0]+'\t'+blast_record.query
                        , file=open("temp_cross_result.csv_1", "a"))
            else:
                pass
    return matched_list


def main():
    args = parse_args()
    doing_blast(args.CDS, args.REF)
    print("##########################################")
    print("blasx and get the xml passed")

    if args.RULE_1:
        matched_list = filter_matching_1(args.CUTOFF_1)
        print("using RULE_1")
    else:
        matched_list = filter_matching(args.CUTOFF)
        print("using default RULE")
    matched_list = list(set(matched_list))
    print("##########################################")
    print("parse the xml and get the matched and un_matched list passed")

    
    seq_records=SeqIO.parse(args.CDS,"fasta")
    matched_seq_records=[]
    un_matched_seq_records=[]
    un_matched_list = []
    main_list=[]
    for seq_record in seq_records:
        main_list.append(seq_record.description)
        for matched_name in matched_list:
            if matched_name == seq_record.description:
                new_seq_record = copy.deepcopy(seq_record)
                matched_seq_records.append(new_seq_record)
    
    for main_name in main_list:
        if main_name not in matched_list:
            un_matched_list.append(main_name)     
    #print(un_matched_list)
    #print(seq_records)
    seq_records=SeqIO.parse(args.CDS,"fasta")
    for seq_record in seq_records:
        for un_matched_name in un_matched_list:
            if  un_matched_name == seq_record.description:
                new_seq_record = copy.deepcopy(seq_record)
                un_matched_seq_records.append(new_seq_record)

    
    SeqIO.write(matched_seq_records,"matched_" + args.OUT, "fasta")
    SeqIO.write(un_matched_seq_records,"un_matched_" + args.OUT, "fasta")
    print("##########################################")
    print("extract and output the matched and un_matched fasta file passed")

    print("hi xiao, you filtering work succeed! cheers!")
    subprocess.run(["mv", "temp_cross_result.csv_1", "temp_cross_result.csv"], check=True)

    if args.ANN:
        import pandas as pd
        EPITOPE_ANN_DF=pd.read_csv(args.ANN_FILE, skiprows=1 )
        EPITOPE_LINK=pd.read_csv('temp_cross_result.csv', sep='\t', names=['Epitope_ID','LOCUS_TAG'])
        EPITOPE_LINK_W = EPITOPE_LINK.merge(EPITOPE_ANN_DF, how='left', left_on='Epitope_ID', right_on='Epitope ID')
        EPITOPE_LINK_W.to_csv('temp_cross_result_ann.csv')
        print("##########################################")
        print("annotation done, output in temp_cross_result_ann.csv")
    else:
        print("##########################################")
        print("annotation not selected.")

    

if __name__ == '__main__':
    sys.exit(main())
