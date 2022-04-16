'''
This robot can help to locally blastn the query.fasta to the reference genomes and extract the matched subject seq to matched seq.fa
20220415---update remove the temp_blast.xml file, handle it in memory
-----------uodate tblastn as optional function for protein 
'''
import glob
import os
import subprocess
import sys
import argparse
import io
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq




def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--Q_SEQ', required=True, type=str, metavar='FILENAME', help="the query fasta filename you want to blast")
    parser.add_argument('--R_FOLDER', required=True, type=str, metavar='FILENAME', help="the complete reference genomes folder you want to blast to")
    parser.add_argument('--CUTOFF', default=0.9, type=float, metavar='DEFAULT 0.9', help="the lowest similarity value which classify as matched")
    parser.add_argument('--BLAST', default="blastn", type=str, metavar='command name', help="can be blastn(default) or tblastn")
    parser.add_argument('--OUT', default="DEFAULT_NAME", type=str, metavar='FILENAME', help="Output fasta name, DEFAULT_NAME as default")
    return parser.parse_args()

#1. blast and get the xml
#2. parse the xml and get the matched and un_matched list
#3. check if matched or not
#4. if matched, extract and output matched.fasta



def doing_blast(query, reference, tabname, blastcmd):
    blast_run_xml = subprocess.run([blastcmd, "-query", query, "-subject", reference, "-outfmt", "5"], check=True, capture_output=True)
    xml = blast_run_xml.stdout
    xml_obj = io.BytesIO(xml)
    subprocess.run([blastcmd, "-query", query, "-subject", reference, "-outfmt", "7", "-out", tabname + ".tab"], check=True)
    os.system("mkdir --parents ./temp_tab; mv *.tab ./temp_tab/.")
    return xml_obj
    



def filter_matching(cutoff, xml_obj):
    seqs=[]
    blast_records = NCBIXML.parse(xml_obj)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            if  hsp.identities / blast_record.query_letters >= cutoff and hsp.expect <= 1e-10 :   
                seq = SeqRecord(Seq(hsp.sbjct), id=alignment.hit_def, description='XIAO_ROBOT')
                seqs.append(seq)
            else:
                pass
    return seqs


def extract_matched_seq(query, reference,tabname, blastcmd, cutoff):
    print(' hey xiao i am doing blast for '+ tabname)
    xml_obj = doing_blast(query,reference,tabname, blastcmd)
    print(" hey xiao blast finished")
    seqs = filter_matching(cutoff, xml_obj)
    print(' hey xiao seq extracted for ' + tabname)
    return seqs


def main():
    args = parse_args()
    
    references=glob.glob('./'+args.R_FOLDER+'/*')
    seqs_sum=[]
    for reference in references:
        tabname=reference.split("/")[2]
        seqs=extract_matched_seq(args.Q_SEQ, reference, tabname, args.BLAST, args.CUTOFF)
        seqs_sum.extend(seqs)
    SeqIO.write(seqs_sum, args.OUT, "fasta")
    print('ohhh! Xiao job done!')

    # os.system("mkdir --parents ./temp_tab; mv *.tab ./temp_tab/.")    

        


    


if __name__ == '__main__':
    sys.exit(main())
