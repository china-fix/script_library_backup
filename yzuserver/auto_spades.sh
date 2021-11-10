#!/bin/bash
#SBATCH --job-name=spades                #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=spades.err               #错误输出
#SBATCH --output=spades.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数

#this is a script handling the spades denovo assemblying with clean reads. china-fixing@hotmail.com by Xiao Fei

start_dir=$(pwd)
mkdir assemblys

dir=$(ls ./input)
for folder in $dir
do
	cd $start_dir/input/$folder
	/data2/qcLi/xiaofei/export_bin/SPAdes-3.15.2-Linux/bin/spades.py  -1 "$folder"_1.clean.fq.gz  -2 "$folder"_2.clean.fq.gz  -o $folder.spades  -t 24  --cov-cutoff auto --careful
	cp ./$folder.spades/scaffolds.fasta $start_dir/assemblys/$folder.fasta
done