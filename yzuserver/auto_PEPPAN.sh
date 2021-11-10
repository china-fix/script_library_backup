#!/bin/bash
#SBATCH --job-name=PEPPAN               #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=sbatch.err               #错误输出
#SBATCH --output=sbatch.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数


#annotation then running the PEPPAN pipeline 20211110
#input directory put the contigs of samples
#env nullarbor_only or BiO

start_dir=$(pwd)
mkdir prokka_ann
mkdir PEPPAN_working
mkdir $start_dir/PEPPAN_working/gffs

dir=$(ls ./input)
for fasta in $dir
do
	prokka --outdir $start_dir/prokka_ann/$fasta -cpus 23 --prefix $fasta ./input/$fasta
done	

for fasta in $dir
do
	cp $start_dir/prokka_ann/$fasta/$fasta.gff  $start_dir/PEPPAN_working/gffs/.
done

cd $start_dir/PEPPAN_working
PEPPAN -p OUT  $start_dir/PEPPAN_working/gffs/*.gff -t 23 
PEPPAN_parser -g $start_dir/PEPPAN_working/OUT.PEPPAN.gff -p PAR_OUT -t -a 98 -c 
