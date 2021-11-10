#!/bin/bash
#SBATCH --job-name=epi-report               #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=nullarbor.err               #错误输出
#SBATCH --output=nullarbor.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数

#this is a script handling the nullarbor runing  with clean reads. china-fixing@hotmail.com by Xiao Fei

start_dir=$(pwd)
mkdir nullarbor_working

dir=$(ls ./input)
for folder in $dir
do
        echo -e ""$folder"\t"$start_dir"/input/"$folder"/"$folder"_1.clean.fq.gz\t"$start_dir"/input/"$folder"/"$folder"_2.clean.fq.gz" >> $start_dir/nullarbor_working/working.tab
done
nullarbor.pl --name xiao-epi-report --ref $start_dir/REF.gbff --input $start_dir/nullarbor_working/working.tab --outdir $start_dir/nullarbor_working/OUTPUT --cpus 23 --force --run --assembler shovill --taxoner kraken2
