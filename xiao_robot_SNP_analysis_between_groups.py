"""this is a robot used to analysis the SNP between different strain groups. Tying to find the conserved SNP cases between the groups. 
1. extract specific CDS from strains and combine
2. Alignment with PRANK
3. SNP-sites
4. get score"""

import sys
import subprocess
import argparse
from Bio import SeqIO
import copy
import vcf

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS_LIST', required=True, type=str, metavar='FILENAME', help="the CDS list you want to extract")
    parser.add_argument('--STRAIN_LIST', required=True, type=str, metavar='FILENAME', help="the strain list you want to extract from the fasta files (each file comtain all the CDS of a strain)")
    parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one group")
    parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one group")
    #parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    parser.add_argument('--OUT', default="xiao_robot_SNP_analysis_between_groups", type=str, metavar='directory', help="Output directory name")
    return parser.parse_args()

def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUT], check=True)

    #_1. extract specific CDS from strains and combine according different CDS
    with open(args.CDS_LIST) as cds_list:
        for cds in cds_list.read().splitlines():
            new_seq_records=[]
            with open(args.STRAIN_LIST) as strain_list:
                for strain in strain_list.read().splitlines():
                    seq_records=SeqIO.parse(strain,"fasta")
                    for seq_record in seq_records:
                        if cds == seq_record.description.split('---FIX---')[1]:
                            new_seq_record = copy.deepcopy(seq_record)
                            new_seq_record.id = seq_record.description.split(" ")[1]  # modifiy the fasta title for the downstream use
                            new_seq_record.description = ""                          #
                            new_seq_records.append(new_seq_record)
                        else:
                            pass
            SeqIO.write(new_seq_records,"temp_" + cds.split(' ')[0], "fasta")
            #subprocess.run(["mv", "temp_" + cds.split(' ')[0], "./"+args.OUT], check=True)
    print("step 1. extract specific CDS from strains and combine according different CDS passed")

    #_2.Alignment with PRANK
    with open(args.CDS_LIST) as cds_list:
        for cds in cds_list.read().splitlines():
            subprocess.run(["prank", "-d=" + "temp_" + cds.split(' ')[0], "-o=" + "temp_" + cds.split(' ')[0]], check=True) # 3.7 python can change to capture_output=True
            subprocess.run(["rm", "temp_" + cds.split(' ')[0]], check=True)
            #subprocess.run(["mv", "temp_" + cds.split(' ')[0] + ".best.fas", "./"+args.OUT], check=True)
    print("#############################################################")
    print("step 2.Alignment with PRANK is passed")

    #_3.SNP calling with snp-sites, output the vcf files
    with open(args.CDS_LIST) as cds_list:
        vcf_filter_list=[]
        vcf_unfilter_list=[]
        for cds in cds_list.read().splitlines():
            return_code = subprocess.run(["snp-sites", "-v" , "-otemp_" + cds.split(' ')[0] + ".vcf", "temp_" + cds.split(' ')[0] + ".best.fas"]).returncode
            if float(return_code) == 0:
                subprocess.run(["rm", "temp_" + cds.split(' ')[0] + ".best.fas"], check=True)
                subprocess.run(["mv", "temp_" + cds.split(' ')[0] + ".vcf", "./"+args.OUT], check=True)
                vcf_filter_list.append("temp_" + cds.split(' ')[0] + ".vcf")
            elif float(return_code) == 1:
                print("in "+ cds.split(' ')[0])
                subprocess.run(["rm", "temp_" + cds.split(' ')[0] + ".best.fas"], check=True)
                vcf_unfilter_list.append("temp_" + cds.split(' ')[0] + ".vcf")
            else:
                print(return_code)
                raise Exception("xiao robot meet some erros at the snp-sites command in step 3")
    '''for vcf_filter in vcf_filter_list:
        print(vcf_filter[5:-4], file=open("temp_vcf_filter_list", "a") )  
    for vcf_unfilter in vcf_unfilter_list:
        print(vcf_unfilter[5:-4], file=open("temp_vcf_unfilter_list", "a") ) '''
    print("#############################################################")
    print("step 3.SNP calling with snp-sites, output the vcf files passed")


    #_4.read and analysis the vcf files and get the scores
    # the function of this part is to get compare each SNP between each group and if the average difference is more than 0.75, add score as 1,then sum the scores in each CDS
    for vcf_filter in vcf_filter_list:
        with open("./"+args.OUT + "/"+ vcf_filter) as VCF_FILTER:
            vcf_reader = vcf.Reader(VCF_FILTER)
            analysis_score = 0
            for vcf_record in vcf_reader: 
                with open(args.GROUP_1) as group_1:
                    group_score_1 = 0
                    group_1_num = 0
                    for name_1 in group_1.read().splitlines():
                        name_1 = name_1.split("ashed_")[1] 
                        name_1_follow = vcf_filter[5:-4]
                        score_1 = vcf_record.genotype(name_1 + "---FIX---" + name_1_follow)["GT"]
                        group_score_1 += float(score_1)
                        group_1_num += 1
                with open(args.GROUP_2) as group_2:
                    group_score_2 = 0
                    group_2_num = 0
                    for name_2 in group_2.read().splitlines():
                        name_2 = name_2.split("ashed_")[1] 
                        name_2_follow = vcf_filter[5:-4]
                        score_2 = vcf_record.genotype(name_2 + "---FIX---" + name_2_follow)["GT"]
                        group_score_2 +=float(score_2)
                        group_2_num +=1
                final_score = abs(group_score_1/group_1_num - group_score_2/group_2_num)
                if final_score > 0.75:
                    analysis_score += 1
            #print(vcf_filter[5:-4] + "---FIX---" + str(analysis_score))
            print(vcf_filter[5:-4] + "---FIX---" + str(analysis_score), file=open(args.OUT+".XIAO", "a"))
    for vcf_unfilter in vcf_unfilter_list:
        print(vcf_unfilter[5:-4] + "---FIX---NO_SNP", file=open(args.OUT+".XIAO", "a") )  
    subprocess.run(["rm", "-r", "./"+args.OUT], check=True)
    print("step 4.read and analysis the vcf files and get the scores passed")
    print("please check the file called " + args.OUT + ".XIAO")


if __name__ == '__main__':
    sys.exit(main())