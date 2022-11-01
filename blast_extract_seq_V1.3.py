'''
This robot can help to locally blastn the query.fasta to the reference genomes and extract the matched subject seq to matched seq.fa
20220415---update remove the temp_blast.xml file, handle it in memory
-----------update tblastn as optional function for protein
20220508---update to output a .XIAO file record the summary of match or no-match things
-----------update to output DEFAULT_NAME.file_name which is the same as DEFAULT_NAME but the naming system is based on the file name of reference genome
20221027---update to add group annotation function to replace the Fixed "XIAO_ROBOT" annotation
20221031---upgrade V1.3 to add function to caculate the extracted seq types distribution
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
import pandas as pd



def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--Q_SEQ', required=True, type=str, metavar='FILENAME', help="the query fasta filename you want to blast")
    parser.add_argument('--R_FOLDER', required=True, type=str, metavar='FILENAME', help="the complete reference genomes folder you want to blast to")
    parser.add_argument('--CUTOFF', default=0.9, type=float, metavar='DEFAULT 0.9', help="the lowest similarity value which classify as matched")
    parser.add_argument('--BLAST', default="blastn", type=str, metavar='command name', help="can be blastn(default) or tblastn")
    parser.add_argument('--OUT', default="DEFAULT_NAME", type=str, metavar='FILENAME', help="Output fasta name, DEFAULT_NAME as default")
    parser.add_argument('--GROUP', default="XIAO_ROBOT", type=str, metavar='NAME', help="seqs header infomation in the output fasta file, XIAO_ROBOT as default")
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
    



def filter_matching(cutoff, xml_obj, tabname, seq_description):
    seqs=[]
    seqs_fname=[]
    blast_records = NCBIXML.parse(xml_obj)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            if  hsp.identities / blast_record.query_letters >= cutoff and hsp.expect <= 1e-10 :   
                seq = SeqRecord(Seq(hsp.sbjct), id=alignment.hit_def, description=seq_description)
                seqs.append(seq)
                seq_fname = SeqRecord(Seq(hsp.sbjct), id=tabname, description=seq_description)
                seqs_fname.append(seq_fname)
            else:
                pass
    return seqs, seqs_fname


def extract_matched_seq(query, reference,tabname, blastcmd, cutoff, seq_description):
    print(' hey xiao i am doing blast for '+ tabname)
    xml_obj = doing_blast(query,reference,tabname, blastcmd)
    print(" hey xiao blast finished")
    seqs, seqs_fname = filter_matching(cutoff, xml_obj, tabname, seq_description)
    print(' hey xiao seq extracted for ' + tabname)
    return seqs, seqs_fname

def analysis_seq_types(file_name,group_name):
    cmd = 'seqkit fx2tab '+file_name+'.f_name'+' > '+file_name+'.f_name.tab'
    os.system(cmd)

    TAB_DF=pd.read_csv(file_name+'.f_name.tab', sep='\t', names=['assembly_name','seq'], index_col=False)
    df=TAB_DF.groupby('seq').count().sort_values(by='assembly_name', ascending=0)
    # df['seq-type'] = [i+1 for i in range(len(df))]
    df.reset_index(inplace=True)
    df.rename(columns={'assembly_name':'count' }, inplace=True)
    df['frequency']=df['count']/df['count'].sum()
    df['group']=group_name
    return df
    

def main():
    args = parse_args()
    
    references=glob.glob(args.R_FOLDER+'*')
    seqs_sum=[]
    seqs_sum_fname=[]
    for reference in references:
        tabname=reference.split("/")[-1]
        seqs, seqs_fname =extract_matched_seq(args.Q_SEQ, reference, tabname, args.BLAST, args.CUTOFF, args.GROUP)
        if not seqs:
            print(tabname + "\tunmatch", file=open(args.OUT+".XIAO", "a") )
        else:
            print(tabname + "\tmatch", file=open(args.OUT+".XIAO", "a") )
        seqs_sum.extend(seqs)
        seqs_sum_fname.extend(seqs_fname)
    SeqIO.write(seqs_sum, args.OUT, "fasta")
    SeqIO.write(seqs_sum_fname, args.OUT+".f_name", "fasta")
    df=analysis_seq_types(args.OUT, args.GROUP)
    df.to_csv(args.OUT+"_"+args.GROUP+".csv")

    print('ohhh! Xiao job done!')

    # os.system("mkdir --parents ./temp_tab; mv *.tab ./temp_tab/.")    

        


    


if __name__ == '__main__':
    sys.exit(main())
