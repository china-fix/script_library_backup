#!/bin/bash
#SBATCH --job-name=serovar-report               #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=seqsero.err               #错误输出
#SBATCH --output=seqsero.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数

#this is a script handling the seqsero2 runing  with clean reads. china-fixing@hotmail.com by Xiao Fei

start_dir=$(pwd)
mkdir seqsero_working

dir=$(ls ./input)
cd seqsero_working/
for folder in $dir
do
SeqSero2_package.py -p 24 -t 2 -i "$start_dir"/input/"$folder"/"$folder"_1.clean.fq.gz "$start_dir"/input/"$folder"/"$folder"_2.clean.fq.gz -d "$folder"_serovar
done
