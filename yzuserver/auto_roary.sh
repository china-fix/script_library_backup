#!/bin/bash
#SBATCH --job-name=roary                #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=roary.err               #错误输出
#SBATCH --output=roary.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数

roary -f ./roary_output -e -s -n -cd 95 -p 24 -g 100000 ./input/*.gff
