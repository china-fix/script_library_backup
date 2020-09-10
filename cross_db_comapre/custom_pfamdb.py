"""
Hi there,
This is a hmmer3 dependent python script to constract your own customed pfam profiles database according the proteins you feed to the total pfam databse
1. 
2. 
3. 
4. 

Enjoy!
Xiao Fei (china-fixing@hotmail.com)
PS. if you are debuging the program or want record everying in the screen, use command: python3 xxx.py xxxxxxx 2>&1 | tee screen.log
"""

import sys
import subprocess
import argparse
from Bio import SearchIO
from shutil import move
#import vcf
import os


def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao's robot")
    parser.add_argument('--CUSTOM_PROTEINS', required=True, type=str, metavar='FILENAME', help="the user selected proteins file to help construct the custom_pfam databse")
    parser.add_argument('--HMMDB_PATH', required=True, type=str, metavar='PATH', help="the path of hmm profile database. eg. /home/fix/downloads/pfam.hmmlib")
    parser.add_argument('--OUTPUT', default="custom_pfamdb_out", type=str, metavar='directory', help="Output directory name")
    #parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    parser.add_argument('--DETAIL_DIC', action='store_const', const=True, metavar='get detailed hmms dictornary', help="optional, get more detailed hmms dictorynary")
    parser.add_argument('--CPU', default=1, type=int, metavar='the cpu num you want to use', help="the cup number you want to use during hmmscan process, default is 1")
    parser.add_argument('--HIT_EVALUE', default=0.001, type=int, metavar='evalue thread', help="the evalue thread use during extracting hmm profile in domain-table, default is 0.001")
    parser.add_argument('--BIT_SCORE', default=10, type=int, metavar='bitscore thread', help="the bitscore thread use during extracting hmm profile in domain-table, default is 10")
    return parser.parse_args()

def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUTPUT], check=True)

    ### hmmscan -o ./human_readable_hmmscan.txt --noali --cpu 2 --domtblout test.tbl --cut_tc /home/fix/Downloads/deltaBS.hmmlib protein.fasta ###
    subprocess.run(["hmmscan", "-o", "./"+args.OUTPUT+"/human_readable_hmmscan.txt", "--noali", "--cpu", str(args.CPU), "--domtblout", "./"+args.OUTPUT+"/temp.hmmscan3-domtab", "--cut_tc", args.HMMDB_PATH, args.CUSTOM_PROTEINS], check=True)
    print("----\n" "hi xiao, hmmscan command passed\n" "----")

    ### parse the hmmer3-domtab ###
    hmmscan_list =[]
    dom_table = SearchIO.parse('./'+args.OUTPUT+'/temp.hmmscan3-domtab', 'hmmscan3-domtab')
    for queryresult in dom_table:
        for hit in queryresult:
            if hit.evalue <= args.HIT_EVALUE and hit.bitscore >= args.BIT_SCORE:
                hmmscan_list.append((hit.accession, queryresult.id ))
            else:
                pass
    print("----\n" "hi xiao, i have read and selected the information you wanted form temp.hmmscan3-domtab\n" "----")
    #print(hmmscan_list)


    ### from list to dictornary ###
    d = dict()
    for x,y in hmmscan_list:
        if x in d:
            d[x].append(y)
        else:
            d[x] = [y]
    print("----\n" "hi xiao, dictornary done, key is the name of each protein you inputed\n" "----")    
    #print(d)


    ### use hmmfetch and dictornay-key to extract the pfam profile one by one and add custom protein information to each###
    for key in d.keys():
        subprocess.run(["hmmfetch", "-o", "./"+args.OUTPUT+"/"+key+".hmm", args.HMMDB_PATH, key], check=True)
        proteins_elements_str = ','.join(str(elem) for elem in d[key])
        temp = open('temp', 'w')
        with open("./"+args.OUTPUT+"/"+key+".hmm", 'r') as f:
            for line in f:
                if line.startswith('DESC'):
                    line = line.strip() + '---'+proteins_elements_str+'\n'
                temp.write(line)
        temp.close()
        move('temp', "./"+args.OUTPUT+"/"+key+".hmm")
    print("----\n""hi xiao, the customed hmm files was done\n""----")

    ### cat the hmm file together, make the final hmm library 'OUPUT.hmmlib' ###
    subprocess.run(["mkdir", "./"+args.OUTPUT+"/hmms"], check=True)
    for key in d.keys():
        cmd = "cat "+"./"+args.OUTPUT+"/"+key+".hmm "+">> "+"./"+args.OUTPUT+"/"+args.OUTPUT+".hmmlib"
        print(cmd)
        os.system(cmd)
        subprocess.run(["mv", "./"+args.OUTPUT+"/"+key+".hmm", "./"+args.OUTPUT+"/hmms/"], check=True)
    print("====\n""Dear Xiao, the customed hmmlib was done! Enjoy!\n""====")
        
    ### save the dictronary ###
    print("ACC\tHIT_PROTEINS", file=open("./"+args.OUTPUT+"/hmms_dictorynary", "a"))
    for key in d.keys():
        proteins_elements_str_1 = ','.join(str(elem) for elem in d[key])
        print(key+'\t'+proteins_elements_str_1, file=open("./"+args.OUTPUT+"/hmms_dictorynary", "a"))


    ### optional, read temp.hmmscan3-domtab and output more information as designed###
    if args.DETAIL_DIC:
        hmmscan_full =[]
        dom_table = SearchIO.parse('./'+args.OUTPUT+'/temp.hmmscan3-domtab', 'hmmscan3-domtab')
        for queryresult in dom_table:
            for hit in queryresult:
                if hit.evalue <= args.HIT_EVALUE and hit.bitscore >= args.BIT_SCORE:
                    hmmscan_full.append(((hit.accession, hit.description), (queryresult.id, hit.bitscore,hit.evalue, hit.bias)))
                else:
                    pass
        #from list to dictornary
        d2 = dict()
        for x,y in hmmscan_full:
            if x in d2:
                d2[x].append(y)
            else:
                d2[x] = [y]   
        #print(d2)
        #save the detailed dictornary
        print("ACC\tDESCRIPTION\tHIT_PROTEINS_NUM\tHIT_PROTEINS(name,bitscore,evalue,bias)", file=open("./"+args.OUTPUT+"/hmms_detailed_dictorynary", "a"))
        for key in d2.keys():
            proteins_elements_str_2 = ';'.join(str(elem) for elem in d2[key])
            print(key[0]+'\t'+key[1]+'\t'+str(len(d2[key]))+'\t'+proteins_elements_str_2, file=open("./"+args.OUTPUT+"/hmms_detailed_dictorynary", "a"))
        print("****\nHi xiao, as your suggestion, the detailed dictornary was added!\n****")
    else:
        pass



        
if __name__ == '__main__':
    sys.exit(main())