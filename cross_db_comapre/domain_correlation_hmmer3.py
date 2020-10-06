"""
Hi there,
This is a hmmer3 dependent python script to correlate the proteins (references proteins/query proteins) according the pfam database (or other HMM_profile)
1. hmmscan the proteins to get the hmmer3-domtab
2. parse the hmmer3-domatab
3. correlated the two hmmer3-tabs according the pfam ACC key

Enjoy!
Xiao Fei (china-fixing@hotmail.com)
PS. if you are debuging the program or want record everying in the screen, use command: python3 xxx.py xxxxxxx 2>&1 | tee screen.log
"""

import sys
import subprocess
import argparse
from Bio import SearchIO
import pandas as pd
#from shutil import move
#import vcf
#import os

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao's robot")
    parser.add_argument('--Q_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as query proteins")
    parser.add_argument('--R_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as reference proteins")
    parser.add_argument('--HMMDB_PATH', required=True, type=str, metavar='PATH', help="the path of hmm profile database. eg. /home/fix/downloads/pfam.hmmlib")
    parser.add_argument('--OUTPUT', default="domain_cor_hmmer3", type=str, metavar='directory', help="Output directory name")
    #parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    parser.add_argument('--PFAM', action='store_const', const=True, metavar='add --cut_tc option in hmmscan specific for pfam_db', help="optional, specific designed for pfam_db")
    parser.add_argument('--CPU', default=1, type=int, metavar='the cpu num you want to use', help="the cup number you want to use during hmmscan process, default is 1")
    parser.add_argument('--HIT_EVALUE', default=0.001, type=int, metavar='evalue thread', help="the evalue thread use during extracting hmm profile in domain-table, default is 0.001")
    parser.add_argument('--BIT_SCORE', default=10, type=int, metavar='bitscore thread', help="the bitscore thread use during extracting hmm profile in domain-table, default is 10")
    return parser.parse_args()

### parse the hmmer3-domtab to panda object ### 
def domtab_to_pandatab(input, HIT_EVALUE, BIT_SCORE):
    hmmscan_list = []
    domtab = SearchIO.parse(input, 'hmmscan3-domtab')
    for queryresult in domtab:
        for hit in queryresult:
            if hit.evalue <= HIT_EVALUE and hit.bitscore >= BIT_SCORE:
                hmmscan_list.append(((hit.accession, hit.description, hit.id, hit.seq_len),(queryresult.id, queryresult.description, queryresult.seq_len, hit.bitscore, hit.evalue, hit.bias,len(hit))))
            else:
                pass
    #from list to dictornary
    d = dict()
    for x,y in hmmscan_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y]           
    #write the dictornary to tab csv temp file
    print("ACC\tDESCRIPTION\tHIT_ID\tHIT_SEQ_LEN\tHIT_PROTEINS(name;description;seq_len;bitscore;evalue;bias;domain_num)\tHIT_PROTEINS\tHIT_PROTEINS_NUM", file=open("temp.csv", "a"))
    for key in d.keys():
            proteins_elements_str = ';'.join(str(elem) for elem in d[key])
            proteins_element_str = ';'.join(elem[0] for elem in d[key])
            print(key[0]+'\t'+key[1]+'\t'+key[2]+'\t'+str(key[3])+'\t'+proteins_elements_str +'\t' +proteins_element_str +'\t' + str(len(d[key])), file=open("temp.csv", "a"))
    # temp.csv:  ACC   DESCRIPTION   HIT_ID   HIT_SEQ_LEN   HIT_PROTEINS   HIT_PROTEINS_SIMPLE   HIT_PROTEINS_NUM
    # read the temp.csv and out put the panda object
    pandatab = pd.read_csv("temp.csv", sep='\t')
    subprocess.run(["rm", "temp.csv"], check=True)
    return pandatab

### parse the hmmer3-domtab to dictornary (designed for the Q_R_cross_analysis) ###
def domtab_to_dictornary(input, HIT_EVALUE, BIT_SCORE):
    hmmscan_list = []
    domtab = SearchIO.parse(input, 'hmmscan3-domtab')
    for queryresult in domtab:
        for hit in queryresult:
            if hit.evalue <= HIT_EVALUE and hit.bitscore >= BIT_SCORE:
                hmmscan_list.append(((hit.accession, hit.description, hit.id, hit.seq_len),(queryresult.id, queryresult.description, queryresult.seq_len, hit.bitscore, hit.evalue, hit.bias,len(hit))))
            else:
                pass
    #from list to dictornary
    d = dict()
    for x,y in hmmscan_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y]     
    return d

### analysis the query and reference dictornary to crosse analysis the functional related proteins and caculate the delta bitscore ###
def Q_R_cross_analysis(dic_Q, dic_R, out_diretory):
    print('AAC\tDESCRIPTION\tDOMAIN_ID\tDOMAIN_LEN\tQ_ID\tQ_BITSCORE\tQ_DOMAIN_NUM\tR_ID\tR_BITSCORE\tR_DOMAIN_NUM\tDELTA_SCORE\tDETAIL_Q(name;description;seq_len;bitscore;evalue;bias;domain_num)\tDETAIL_R(name;description;seq_len;bitscore;evalue;bias;domain_num)'
                    , file=open("./"+out_diretory+"/cross_result.csv", "a"))
    for key_Q in dic_Q.keys():
        if key_Q in dic_R.keys():
            for elem_Q in dic_Q[key_Q]:
                for elem_R in dic_R[key_Q]:
                    # ACC   DESCRIPTION   HIT_ID   HIT_SEQ_LEN   Q_ID   Q_BITSCORE   Q_DOMAIN_NUM   R_ID   R_BITSCORE   R_DOMAIN_NUM   DELTA_SCORE   DETAIL_Q()   DETAIL_R()
                    ## DETAIL (name;description;seq_len;bitscore;evalue;bias;domain_num)
                    print(key_Q[0]+'\t'+key_Q[1]+'\t'+key_Q[2]+'\t'+str(key_Q[3])+'\t'+elem_Q[0]+'\t'+str(elem_Q[3])+'\t'+str(elem_Q[6])+'\t'+elem_R[0]+'\t'+str(elem_R[3])+'\t'+str(elem_R[6])+'\t'+str(elem_Q[3]-elem_R[3])+'\t'+str(elem_Q)+'\t'+str(elem_R)
                    , file=open("./"+out_diretory+"/cross_result.csv", "a"))
        else:
            pass


'''### subset pd into subset pd according pd['DB']
def subset (pdobject):
    DBs=['CDD','COILS','GENE3D','HAMAP','MOBIDB_LITE','PFAM','PIRSF','PRINTS','PROSITE_PATTERNS','PROSITE_PROFILES','SFLD','SMART','SUPERFAMILY','TIGRFAM']
    pds=[]
    for DB in DBs:
        subset_pd = pdobject[pdobject['DB'] == DB]
        subset_pd_record = [DB,subset_pd]
        pds.append(subset_pd_record)
    return pds'''


### extract pd['HIT_PROTEINS'] column to list according proteins no duplicate
def hit_proteins_sum(pdobject, column_name, proteins_txt):
    protein_list=[]
    column_list=pdobject[column_name].tolist()
    for proteins in column_list:
        for protein in proteins.split(';'):
            if protein in protein_list:
                pass
            else:
                protein_list.append(protein)
    with open(proteins_txt, "w") as f:
        for s in protein_list:
            f.write(str(s) +"\n")




def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUTPUT], check=True)

    ### hmmscan -o ./human_readable_hmmscan.txt --noali --cpu 2 --domtblout test.tbl --cut_tc /home/fix/Downloads/deltaBS.hmmlib protein.fasta ###
    if args.PFAM:
        subprocess.run(["hmmscan", "-o", "./"+args.OUTPUT+"/human_readable_hmmscan.R", "--noali", "--cpu", str(args.CPU), "--domtblout", "./"+args.OUTPUT+"/temp-R.hmmscan3-domtab", "--cut_tc", args.HMMDB_PATH, args.R_PROTEINS], check=True)
        print("----\n" "hi xiao, hmmscan R_PROTEINS command passed\n" "----")

        subprocess.run(["hmmscan", "-o", "./"+args.OUTPUT+"/human_readable_hmmscan.Q", "--noali", "--cpu", str(args.CPU), "--domtblout", "./"+args.OUTPUT+"/temp-Q.hmmscan3-domtab", "--cut_tc", args.HMMDB_PATH, args.Q_PROTEINS], check=True)
        print("----\n" "hi xiao, hmmscan Q_PROTEINS command passed\n" "----")
    else:
        subprocess.run(["hmmscan", "-o", "./"+args.OUTPUT+"/human_readable_hmmscan.R", "--noali", "--cpu", str(args.CPU), "--domtblout", "./"+args.OUTPUT+"/temp-R.hmmscan3-domtab",  args.HMMDB_PATH, args.R_PROTEINS], check=True)
        print("----\n" "hi xiao, hmmscan R_PROTEINS command passed\n" "----")

        subprocess.run(["hmmscan", "-o", "./"+args.OUTPUT+"/human_readable_hmmscan.Q", "--noali", "--cpu", str(args.CPU), "--domtblout", "./"+args.OUTPUT+"/temp-Q.hmmscan3-domtab",  args.HMMDB_PATH, args.Q_PROTEINS], check=True)
        print("----\n" "hi xiao, hmmscan Q_PROTEINS command passed\n" "----")

    ### get the pandatabs ###
    pandatab_R = domtab_to_pandatab(args.OUTPUT+"/temp-R.hmmscan3-domtab", args.HIT_EVALUE, args.BIT_SCORE)
    pandatab_Q = domtab_to_pandatab(args.OUTPUT+"/temp-Q.hmmscan3-domtab", args.HIT_EVALUE, args.BIT_SCORE)

    ### merge pandatabs and output ###
    pandatab_RQ= pd.merge( pandatab_R, pandatab_Q, how= 'outer', on=['ACC', 'DESCRIPTION', 'HIT_ID', 'HIT_SEQ_LEN'], suffixes=['_R','_Q'])
    pandatab_RQ_inner= pd.merge( pandatab_R, pandatab_Q, how= 'inner', on=['ACC', 'DESCRIPTION', 'HIT_ID', 'HIT_SEQ_LEN'], suffixes=['_R','_Q'])
    pandatab_RQ.to_csv("./"+args.OUTPUT+"/domain_correlation.csv", sep='\t') 
    pandatab_RQ_inner.to_csv("./"+args.OUTPUT+"/domain_correlation_inner.csv", sep='\t') 
    

    ### get the inner sum protein list and output
    subprocess.run(["mkdir", "./"+args.OUTPUT+"/HITPROTEINS"], check=True)
    '''pandatab_RQ_inner_subsets=subset(pandatab_RQ_inner)'''
    '''for pandatab_RQ_inner_subset in pandatab_RQ_inner_subsets:'''
    hit_proteins_sum(pandatab_RQ_inner,'HIT_PROTEINS_R', "./"+args.OUTPUT+"/HITPROTEINS"+"/"+"_HITPROTEINS_R")
    hit_proteins_sum(pandatab_RQ_inner,'HIT_PROTEINS_Q', "./"+args.OUTPUT+"/HITPROTEINS"+"/"+"_HITPROTEINS_Q")



    ### get the dictornary ###
    dic_R = domtab_to_dictornary(args.OUTPUT+"/temp-R.hmmscan3-domtab", args.HIT_EVALUE, args.BIT_SCORE)
    dic_Q = domtab_to_dictornary(args.OUTPUT+"/temp-Q.hmmscan3-domtab", args.HIT_EVALUE, args.BIT_SCORE)

    ### analysis the dictornary ###
    Q_R_cross_analysis(dic_Q, dic_R, out_diretory=args.OUTPUT)
    #subprocess.run(["mv", "temp_cross_result.csv","./"+args.OUTPUT+"/cross_result.csv"], check=True)
    print("****\nhi xiao, the correlation analysis was done! Enjoy!\n****")

    ### backup the input query_proteins and reference_proteins
    subprocess.run(["mkdir", "./"+args.OUTPUT+"/INPUT_BACKUP"], check=True)
    subprocess.run(["cp", args.R_PROTEINS, args.Q_PROTEINS, "./"+args.OUTPUT+"/INPUT_BACKUP/"], check=True)
    print("****\nhi xiao, I aslo backup the input files in INPUT_BACUP folder, enjoy!\n****")

if __name__ == '__main__':
    sys.exit(main())