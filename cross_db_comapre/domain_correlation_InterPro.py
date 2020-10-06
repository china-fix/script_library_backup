"""
Hi there,
This is a InterProScan dependent python script to correlate the proteins (references proteins/query proteins) according the different databases.
1. InterProScan the proteins to get the hmmer3-domtab
2. parse the XML file
3. correlated the two xml according the ACC&DB key

Enjoy!
Xiao Fei (china-fixing@hotmail.com)
PS. if you are debuging the program or want record everying in the screen, use command: python3 xxx.py xxxxxxx 2>&1 | tee screen.log
"""

import sys
import subprocess
import argparse
from Bio.SearchIO.InterproscanIO import interproscan_xml
import pandas as pd
#import os

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao's robot")
    parser.add_argument('--Q_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as query proteins")
    parser.add_argument('--R_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as reference proteins")
    parser.add_argument('--INTERPROSCAN_BIN', required=True, type=str, metavar='PATH', help="the path of hmm profile InterProScan. eg. /home/fix/export_bin/InterProScan/interproscan-5.46-81.0/interproscan.sh")
    parser.add_argument('--OUTPUT', default="domain_cor_InterPro", type=str, metavar='directory', help="Output directory name")
    #parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    #parser.add_argument('--DETAIL_DIC', action='store_const', const=True, metavar='get detailed hmms dictornary', help="optional, get more detailed hmms dictorynary")
    #parser.add_argument('--CPU', default=1, type=int, metavar='the cpu num you want to use', help="the cup number you want to use during hmmscan process, default is 1")
    #parser.add_argument('--HIT_EVALUE', default=0.001, type=int, metavar='evalue thread', help="the evalue thread use during extracting hmm profile in domain-table, default is 0.001")
    #parser.add_argument('--BIT_SCORE', default=10, type=int, metavar='bitscore thread', help="the bitscore thread use during extracting hmm profile in domain-table, default is 10")
    return parser.parse_args()

### parse the scan.xml to panda object ### 
def scan_to_pandatab(input):
    scan_list = []
    scan_results = interproscan_xml.InterproscanXmlParser(input)
    for queryresult in scan_results:
        for hit in queryresult:
            scan_list.append(((hit.id, hit.attributes['Target'],hit.description),(hit.query_id,len(hit.hsps))))
            
    #from list to dictornary
    d = dict()
    for x,y in scan_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y]           
    #write the dictornary to tab csv temp file
    print("ACC\tDB\tDESCRIPTION\tHIT_PROTEINS(name,domain_num)\tHIT_PROTEINS\tHIT_PROTEINS_NUM", file=open("temp.csv", "a"))
    for key in d.keys():
            proteins_elements_str = ';'.join(str(elem) for elem in d[key])
            proteins_element_str = ';'.join(elem[0] for elem in d[key])
            print(key[0]+'\t'+key[1]+'\t'+key[2]+'\t'+proteins_elements_str +'\t' +proteins_element_str +'\t' + str(len(d[key])), file=open("temp.csv", "a"))
    # temp.csv:  ACC   DB   DESCRIPTION   HIT_PROTEINS   HIT_PROTEINS_SIMPLE   HIT_PROTEINS_NUM
    # read the temp.csv and out put the panda object
    pandatab = pd.read_csv("temp.csv", sep='\t')
    #hitproteins_list= pandatab['HIT_PROTEINS']
    subprocess.run(["rm", "temp.csv"], check=True)
    return pandatab

### parse the scan.xml to dictornary (designed for the Q_R_cross_analysis) ###
def scan_to_dictornary(input):
    scan_list = []
    scan_results = interproscan_xml.InterproscanXmlParser(input)
    for queryresult in scan_results:
        for hit in queryresult:
            a=[]
            for hsp in hit.hsps:
                try:
                    a.append(hsp.bitscore)
                except AttributeError:
                    a.append(None)
            try:
                sum_bitscore=sum(a)
            except TypeError:
                sum_bitscore=None


            scan_list.append(((hit.id, hit.attributes['Target'],hit.description),(hit.query_id,sum_bitscore,len(hit.hsps))))
            
    #from list to dictornary
    d = dict()
    for x,y in scan_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y]   
    return d
    

### analysis the query and reference dictornary to crosse analysis the functional related proteins and caculate the delta bitscore ###
def Q_R_cross_analysis(dic_Q, dic_R, out_diretory):
    print('AAC\tDB\tDESCRIPTION\tQ_ID\tQ_BITSCORE\tQ_DOMAIN_NUM\tR_ID\tR_BITSCORE\tR_DOMAIN_NUM\tDELTA_SCORE'
                    , file=open("./"+out_diretory+"/cross_result.csv", "a"))
    for key_Q in dic_Q.keys():
        if key_Q in dic_R.keys():
            for elem_Q in dic_Q[key_Q]:
                for elem_R in dic_R[key_Q]:
                    
                    ###### ACC   DB   DESCRIPTION   Q_ID   Q_BITSCORE   Q_DOMAIN_NUM   R_ID   R_BITSCORE   R_DOMAIN_NUM   DELTA_SCORE   
                    try:
                        print(key_Q[0]+'\t'+key_Q[1]+'\t'+key_Q[2]+'\t'+elem_Q[0]+'\t'+str(elem_Q[1])+'\t'+str(elem_Q[2])+'\t'+elem_R[0]+'\t'+str(elem_R[1])+'\t'+str(elem_R[2])+'\t'+str(elem_Q[1]-elem_R[1])
                        , file=open("./"+out_diretory+"/cross_result.csv", "a"))
                    
                    except TypeError:
                        print(key_Q[0]+'\t'+key_Q[1]+'\t'+key_Q[2]+'\t'+elem_Q[0]+'\t'+str(elem_Q[1])+'\t'+str(elem_Q[2])+'\t'+elem_R[0]+'\t'+str(elem_R[1])+'\t'+str(elem_R[2])+'\t'+str(None)
                        , file=open("./"+out_diretory+"/cross_result.csv", "a"))


        else:
            pass

### subset pd into subset pd according pd['DB']
def subset (pdobject):
    DBs=['CDD','COILS','GENE3D','HAMAP','MOBIDB_LITE','PFAM','PIRSF','PRINTS','PROSITE_PATTERNS','PROSITE_PROFILES','SFLD','SMART','SUPERFAMILY','TIGRFAM']
    pds=[]
    for DB in DBs:
        subset_pd = pdobject[pdobject['DB'] == DB]
        subset_pd_record = [DB,subset_pd]
        pds.append(subset_pd_record)
    return pds



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


    ### /home/fix/export_bin/InterProScan/interproscan-5.46-81.0/interproscan.sh -b ./interpro/test -i protein.fasta -iprlookup -pa ###
    subprocess.run([args.INTERPROSCAN_BIN, "-b", "./"+args.OUTPUT+"/R_scan", "-i", args.R_PROTEINS, "-iprlookup", "-pa", "-goterms"], check=True)
    print("----\n" "hi xiao, InterProScan R_PROTEINS command passed\n" "----")

    subprocess.run([args.INTERPROSCAN_BIN, "-b", "./"+args.OUTPUT+"/Q_scan", "-i", args.Q_PROTEINS, "-iprlookup", "-pa", "-goterms"], check=True)
    print("----\n" "hi xiao, InterProScan Q_PROTEINS command passed\n" "----")

    ### get the pandatabs ###
    pandatab_R = scan_to_pandatab(args.OUTPUT+"/R_scan.xml")
    pandatab_Q = scan_to_pandatab(args.OUTPUT+"/Q_scan.xml")

    ### merge pandatabs and output ###
    pandatab_RQ= pd.merge( pandatab_R, pandatab_Q, how= 'outer', on=['ACC', 'DB', 'DESCRIPTION'], suffixes=['_R','_Q'])
    pandatab_RQ_inner= pd.merge( pandatab_R, pandatab_Q, how= 'inner', on=['ACC', 'DB', 'DESCRIPTION'], suffixes=['_R','_Q'])
    pandatab_RQ.to_csv("./"+args.OUTPUT+"/domain_correlation.csv", sep='\t') 
    pandatab_RQ_inner.to_csv("./"+args.OUTPUT+"/domain_correlation_inner.csv", sep='\t') 
    
    ### get the inner sum protein list and output
    subprocess.run(["mkdir", "./"+args.OUTPUT+"/HITPROTEINS"], check=True)
    pandatab_RQ_inner_subsets=subset(pandatab_RQ_inner)
    for pandatab_RQ_inner_subset in pandatab_RQ_inner_subsets:
        hit_proteins_sum(pandatab_RQ_inner_subset[1],'HIT_PROTEINS_R', "./"+args.OUTPUT+"/HITPROTEINS"+"/"+pandatab_RQ_inner_subset[0]+"_HITPROTEINS_R")
        hit_proteins_sum(pandatab_RQ_inner_subset[1],'HIT_PROTEINS_Q', "./"+args.OUTPUT+"/HITPROTEINS"+"/"+pandatab_RQ_inner_subset[0]+"_HITPROTEINS_Q")
    
    
    
    ### get the dictornary ###
    dic_R = scan_to_dictornary(args.OUTPUT+"/R_scan.xml")
    dic_Q = scan_to_dictornary(args.OUTPUT+"/Q_scan.xml")

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