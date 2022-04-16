####
#This script is designed to get the metadata of Biosample ID,written by fix,any question,please sent email to fix(china-fixing@hotmmail.com)
##The metadata contain 1.isolation source, 2.geographic location, 3.collection data, 4.collected by, 5.collected by, 6.specific host 7.host 8.platform 
####

import os,sys
from Bio import Entrez
import argparse
import subprocess
import time
Entrez.email = "China-fixing@hotmail.com"

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--BIOSAMPLE_LIST', required=True, type=str, metavar='FILENAME', help="the Biosample id list you want to search, eg. accid\tbiosampleid")
    #parser.add_argument('--STRAIN_LIST', required=True, type=str, metavar='FILENAME', help="the strain list you want to extract from the fasta files (each file comtain all the CDS of a strain)")
    #parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    #parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    #parser.add_argument('--FAST', action='store_const', const=True, metavar='FAST ALIGNMENT WITH MAFFT', help="this command help to use MAFFT to do a fast alignment")
    #parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    parser.add_argument('--OUT', default="get_Biosample_meta", type=str, metavar='directory', help="Output file/directory prefix name")
    return parser.parse_args()







def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUT], check=True)
    print('biosample_acc'+"\t"+'biosample_id'+"\t"+'collection_date'+"\t"+'geo_loc_name'+"\t"+'host'+"\t"+'isolation_source'+"\t"+'strain'+"\t"+'organismname', file=open("./"+args.OUT+"/"+args.OUT+".tab", "w") ) 
#input Biosample IDs
    with open(args.BIOSAMPLE_LIST) as biosample_list:
            for biosample in biosample_list.read().splitlines():
                biosample_acc = biosample.split("\t")[0]
                biosample_id = biosample.split("\t")[1]
                
                tries = 10
                for i in range(tries):
                    while True:
                        try:
                            search_handle = Entrez.esearch(db="Biosample", term=biosample_id+"[Accession]") #esearch to get id for efetch
                        except:
                            if i < tries -1:
                                print("hey xiao,maybe esearch network problem, i will retry!")
                                time.sleep(60)
                                continue
                            else:
                                raise
                        break
                
                search_record = Entrez.read(search_handle)
                search_handle.close()
                try:
                    efetch_id = search_record['IdList'][0]
                except IndexError:
                    efetch_id = biosample_id
                
                tries = 10
                for i in range(tries):
                    while True:
                        try:
                            efetch_handle=Entrez.efetch(db="Biosample", id=efetch_id, rettype="xml", retmode="text") #efetch the meta data
                        except:
                            if i < tries -1:
                                print("hey xiao,maybe efetch network problem, i will retry!")
                                time.sleep(60)
                                continue
                            else:
                                raise
                        break
                
                efetch_record = efetch_handle.read()
                try:
                    efetch_record = efetch_record.decode("utf-8")
                except AttributeError:
                    #print("no decode problem")
                    pass
                efetch_handle.close()
                
                collection_date = "ND"
                geo_loc_name = "ND"
                host ="ND"
                isolation_source = "ND"
                strain = "ND"
                organismname = "ND"
                
                if 'display_name="collection date">' in efetch_record:
                    collection_date = efetch_record.split('display_name="collection date">')[1].split("</Attribute>")[0]
                if 'display_name="geographic location">' in efetch_record:
                    geo_loc_name = efetch_record.split('display_name="geographic location">')[1].split("</Attribute>")[0]
                if 'display_name="host">' in efetch_record:
                    host = efetch_record.split('display_name="host">')[1].split("</Attribute>")[0]
                if 'display_name="isolation source">' in efetch_record:
                    isolation_source = efetch_record.split('display_name="isolation source">')[1].split("</Attribute>")[0]
                if 'display_name="strain">' in efetch_record:
                    strain = efetch_record.split('display_name="strain">')[1].split("</Attribute>")[0]
                if '<OrganismName>' in efetch_record:
                    organismname = efetch_record.split('<OrganismName>')[1].split("</OrganismName>")[0]

                print(biosample_acc+"\t"+biosample_id+"\t"+collection_date+"\t"+geo_loc_name+"\t"+host+"\t"+isolation_source+"\t"+strain+"\t"+organismname, file=open("./"+args.OUT+"/"+args.OUT+".tab", "a") ) 




if __name__ == '__main__':
    sys.exit(main())