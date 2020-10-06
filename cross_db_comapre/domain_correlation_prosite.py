"""
Hi there,
This is a ps_scan dependent python script to correlate the proteins (references proteins/query proteins) according the scanProsite database.
1. ps_scan the proteins to get the pff_table (tabular format listing bounding positions on the sequence and the profile)
2. parse the pff_table
3. correlated the two pff_tables according the prosite database ACC key

Enjoy!
Xiao Fei (china-fixing@hotmail.com)
PS. if you are debuging the program or want record everying in the screen, use command: python3 xxx.py xxxxxxx 2>&1 | tee screen.log
"""

import sys
import subprocess
import argparse
from Bio.ExPASy import Prosite
#from Bio import SearchIO
import pandas as pd
#from shutil import move
#import vcf
import os

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao's robot")
    parser.add_argument('--Q_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as query proteins")
    parser.add_argument('--R_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file as reference proteins")
    parser.add_argument('--PROSITE_DB_PATH', required=True, type=str, metavar='PATH', help="the path of prosite.dat profile. eg. /home/fix/downloads/prosite.dat")
    parser.add_argument('--EVALUATOR_PATH', required=True, type=str, metavar='PATH', help="the path of evaluator.dat profile (for evaluate the pattern). eg. /home/fix/downloads/evaluator.dat")
    parser.add_argument('--ps_scan_PATH', required=True, type=str, metavar='PATH', help="the folder path of ps_scan.pl program. eg. /home/fix/ps_scan/")
    parser.add_argument('--OUTPUT', default="domain_cor_prosite", type=str, metavar='directory', help="Output directory name")
    #parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    parser.add_argument('--ONLY_DOM', action='store_const', const=True, metavar='ONLY_DOMMAIN_DB', help="optional,only ps_scan the protein domain profiles")
    parser.add_argument('--ONLY_PAT', action='store_const', const=True, metavar='ONLY_PATTERN_DB', help="optional,only ps_scan the protein pattern profiles")
    parser.add_argument('--SIMPLE_ANN', action='store_const', const=True, metavar='ONLY SHOW THE MATCH ACC', help="optional, in the ouput only show the match profile ID, without detail description")
    #parser.add_argument('--CPU', default=1, type=int, metavar='the cpu num you want to use', help="the cup number you want to use during hmmscan process, default is 1")
    #parser.add_argument('--HIT_EVALUE', default=0.001, type=int, metavar='evalue thread', help="the evalue thread use during extracting hmm profile in domain-table, default is 0.001")
    #parser.add_argument('--BIT_SCORE', default=10, type=int, metavar='bitscore thread', help="the bitscore thread use during extracting hmm profile in domain-table, default is 10")
    return parser.parse_args()



### parse the pff.csv to dictornary object ###
def pff_to_dic(input,simple):
    pff_pd=pd.read_csv(input, sep='\t', names=range(9))   ##read csv in to dataframe
    pff_array=pff_pd.to_numpy()   ##dataframe to array
    ## array to list
    pff_list=[]
    ### q_name[0]  q_from[1]  q_to[2]    ACC[3] r_from[4]  r_to[5]    score[6]   normarized_score[7]    level_tag[8]
    ### cysB	1	59	PS50931	1	-1	80077	25.505	0
    #print(simple)
    if simple is not None:
        for hit in pff_array:
            #print(hit[6])
            if str(hit[6]) == 'nan':
                ACC=hit[3]+'###pattern'
            else:
                ACC=hit[3]+'###profile'
            pff_list.append((ACC,(hit[0],(hit[7],hit[8]))))
    else:
        for hit in pff_array:
            ACC=hit[3]
            pff_list.append((ACC,(hit[0],(hit[7],hit[8]))))
    ### from list to dic
    ### d ={ACC:[(q_name,(normarized_score, level_tag)),()...]...}
    d = dict()
    for x,y in pff_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y] 
    return d

### rebuild the dictornary to make new_dic###
def dic_to_dic2(d):
    d_new =dict()
    for key in d.keys():
        hit_proteins=d[key]
        d2 =dict()
        for x,y in hit_proteins:
            if x in d2:
                d2[x].append(y)
            else:
                d2[x] = [y]
        scores=[]
        for d2_key in d2.keys():
            for protein_scores in d2[d2_key]:
                if protein_scores[0] == 'nan':
                    scores.append(0)
                else:
                    scores.append(protein_scores[0])
            if key in d_new:
                d_new[key].append([d2_key, sum(scores), len(scores)])
            else:
                d_new[key] =[[d2_key, sum(scores), len(scores)]]
    return d_new



### parse dictornary to pandas object ###   
def dic_to_pandas(d,temp_out):   
    ### make a dic inside d[key]
    print('ACC\tHIT_PROTEINS\tHIT_PROTEINS_NUM',file=open(temp_out,"a"))
    for key in d.keys():
        hit_proteins=d[key]
        d2 =dict()
        for x,y in hit_proteins:
            if x in d2:
                d2[x].append(y)
            else:
                d2[x] = [y]
            proteins_elements_str = ';'.join(str(elem) for elem in d2.keys())
            proteins_num =len(d2.keys())
        print(key+'\t'+proteins_elements_str+'\t'+ str(proteins_num),file=open(temp_out, "a"))
    pandatab = pd.read_csv(temp_out, sep='\t')
    #subprocess.run(["rm", temp_out], check=True)
    return pandatab



### analysis the query and reference new_dictornary to crosse analysis the functional related proteins and caculate the delta bitscore ###
def Q_R_cross_analysis(dic2_Q, dic2_R, out_diretory):
    print('AAC\tQ_ID\tQ_SCORE\tQ_DOMAIN_NUM\tR_ID\tR_SCORE\tR_DOMAIN_NUM\tDELTA_SCORE'
                    , file=open("./"+out_diretory+"/cross_result.csv", "a"))
    for key_Q in dic2_Q.keys():
        if key_Q in dic2_R.keys():
            for elem_Q in dic2_Q[key_Q]:
                for elem_R in dic2_R[key_Q]:




                    # ACC   Q_ID   Q_SCORE   Q_DOMAIN_NUM   R_ID   R_SCORE   R_DOMAIN_NUM   DELTA_SCORE
                    print(key_Q+'\t'+elem_Q[0]+'\t'+str(elem_Q[1])+'\t'+str(elem_Q[2])+'\t'+elem_R[0]+'\t'+str(elem_R[1])+'\t'+str(elem_R[2])+'\t'+str(elem_Q[1]-elem_R[1])
                    , file=open("./"+out_diretory+"/cross_result.csv", "a"))
        else:
            pass

### parse the prosite.dat and output pandas object ###
def prosite_to_pandas(input,temp_out):
    # ACC[accession]   TYPE[type]   NAME[name]   DESCRIPTION[description]
    print('ACC\tTYPE\tNAME\tDESCRIPTION',file=open(temp_out,"a"))
    with open(input) as handle:
        prosite_db=Prosite.parse(handle)
        for prosite_record in prosite_db:
            ACC = prosite_record.accession
            TYPE = prosite_record.type
            NAME = prosite_record.name
            DESCRIPTION = prosite_record.description
            print(ACC+'\t'+TYPE+'\t'+NAME+'\t'+DESCRIPTION,file=open(temp_out,"a"))
    pandatab = pd.read_csv(temp_out, sep='\t')
    return pandatab







def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUTPUT], check=True)
    cmd= "export PATH="+args.ps_scan_PATH+":$PATH"
    #print(cmd)
    os.system(cmd)
    #export PATH=~/opt/bin:$PATH
    
    ###  get the pff tables
    ###  /home/fix/Downloads/ps_scan/ps_scan.pl -d /home/fix/Downloads/ps_scan/prosite.dat -b /home/fix/Downloads/ps_scan/evaluator.dat -o pff -s -r?   /home/fix/Desktop/custom_pfamdb/CDS.fasta.p ###
    if args.ONLY_PAT:  #patten profile scan '-r'
        sub1=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s','-r', args.Q_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/Q_pff.csv', 'w') as file:
            file.write(sub1.stdout)
        sub2=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s','-r', args.R_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/R_pff.csv', 'w') as file:
            file.write(sub2.stdout)
    elif args.ONLY_DOM:   #domain profile scan  '-m'
        sub1=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s','-m', args.Q_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/Q_pff.csv', 'w') as file:
            file.write(sub1.stdout)
        sub2=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s','-m', args.R_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/R_pff.csv', 'w') as file:
            file.write(sub2.stdout)
    else:   #else, scan both pattern and domain
        print("----\nhi xiao, domain or pattern profile was not defined, so i will scan both\n----")
        sub1=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s', args.Q_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/Q_pff.csv', 'w') as file:
            file.write(sub1.stdout)
        sub2=subprocess.run(['ps_scan.pl','-d', args.PROSITE_DB_PATH, '-b', args.EVALUATOR_PATH, '-o','pff', '-s', args.R_PROTEINS], capture_output=True, check=True, text=True)
        #print(sub1.stdout)
        with open("./"+args.OUTPUT+'/R_pff.csv', 'w') as file:
            file.write(sub2.stdout)
    print("----\n" "hi xiao, ps_scan command passed and get the pff tables\n" "----")

    ### get the pandatabs ###
    pff_dic_Q = pff_to_dic("./"+args.OUTPUT+'/Q_pff.csv', args.SIMPLE_ANN)
    pandatab_Q = dic_to_pandas(pff_dic_Q,"./"+args.OUTPUT+'/Q_temp_pandatab.csv')

    pff_dic_R = pff_to_dic("./"+args.OUTPUT+'/R_pff.csv', args.SIMPLE_ANN)
    pandatab_R = dic_to_pandas(pff_dic_R,"./"+args.OUTPUT+'/R_temp_pandatab.csv')
    

    ### merge pandatabs and output ###
    pandatab_RQ= pd.merge( pandatab_R, pandatab_Q, how= 'outer', on='ACC', suffixes=['_R','_Q'])
    pandatab_RQ_inner= pd.merge( pandatab_R, pandatab_Q, how= 'inner', on='ACC', suffixes=['_R','_Q'])
    pandatab_RQ.to_csv("./"+args.OUTPUT+"/domain_correlation.csv", sep='\t') 
    pandatab_RQ_inner.to_csv("./"+args.OUTPUT+"/domain_correlation_inner.csv", sep='\t') 
    

    
    ### get the new_dictornary ###
    pff_dic_Q_new = dic_to_dic2(pff_dic_Q)
    pff_dic_R_new = dic_to_dic2(pff_dic_R)

    ### analysis the dictornary ###
    Q_R_cross_analysis(pff_dic_Q_new, pff_dic_R_new, out_diretory=args.OUTPUT)
    #subprocess.run(["mv", "temp_cross_result.csv","./"+args.OUTPUT+"/cross_result.csv"], check=True)



    ### add the annotation information to the results ###
    if args.SIMPLE_ANN is None:    
        # get the prosite pandatab
        pandatab_prosite=prosite_to_pandas(args.PROSITE_DB_PATH,"./"+args.OUTPUT+'/temp_prosite_db.csv')
        # merge
        pandatab_RQ_ann= pd.merge(pandatab_RQ, pandatab_prosite, how='left', on='ACC')
        pandatab_RQ_inner_ann= pd.merge(pandatab_RQ_inner, pandatab_prosite, how='left', on='ACC')
        pandatab_RQ_ann.to_csv("./"+args.OUTPUT+"/domain_correlation_ann.csv", sep='\t') 
        pandatab_RQ_inner_ann.to_csv("./"+args.OUTPUT+"/domain_correlation_inner_ann.csv", sep='\t') 
    else:
        print("----\nhi xiao, remeber you have chosen the SIMPLE_ANN mode...so, simply annotated with ACC number,you can use this number search online by yourself!\n----")


    print("****\nhi xiao, the correlation analysis was done! Enjoy!\n****")


    ### backup the input query_proteins and reference_proteins
    subprocess.run(["mkdir", "./"+args.OUTPUT+"/INPUT_BACKUP"], check=True)
    subprocess.run(["cp", args.R_PROTEINS, args.Q_PROTEINS, "./"+args.OUTPUT+"/INPUT_BACKUP/"], check=True)
    print("****\nhi xiao, I aslo backup the input files in INPUT_BACUP folder, enjoy!\n****")

if __name__ == '__main__':
    sys.exit(main())