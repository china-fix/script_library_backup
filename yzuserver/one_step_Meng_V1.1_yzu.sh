#!/bin/bash
#SBATCH --job-name=one_step_Meng                #指定脚本名
#SBATCH --partition=hqueue            #指定队列名
#SBATCH --error=one_step_Meng.err               #错误输出
#SBATCH --output=one_step_Meng.out          #结果输出
#SBATCH -N 1                                   #指定节点数
#SBATCH -n 24                                 #指定核数


#this is a main one step script combined with different steps. china-fixing@hotmail.com by Xiao Fei
#20210801 modified to add traits preparing for pyseer by Xiao

start_dir=$(pwd)
cp -r ./input ./input_assemblys
mkdir temp_gffs

###using the filter_short_contigs.py
dir=$(ls ./input)
for folder in $dir
do
	cd $start_dir/input_assemblys/$folder
	mkdir temp_drop
	for j in *
	do
		python /data2/qcLi/xiaofei/export_bin/fix_kits/filter_short_contigs.py --original_fasta $j --output $j.clean --length 500  ##here we set to 500bp to drop
		mv $j.clean.drop ./temp_drop/$j.drop
	done

### rename contigs to avoiding prokka annotation errors, rename_contigs.py
	mkdir temp_rename_save
	for i in *.clean
	do
		python3 /data2/qcLi/xiaofei/export_bin/fix_kits/rename_contigs.py --original_fasta $i --output $i.rename
		cp $i.rename ./temp_rename_save/.
	done

### annotation with prokka
	mkdir $folder.gffs
	for i in *.rename
	do
		prokka  --outdir temp_fix --prefix mygenome  --cpus 24 --force $i
		cp temp_fix/mygenome.gff $folder.gffs/$i.gff
		cp temp_fix/mygenome.gff $start_dir/temp_gffs/$i.gff
		rm -r temp_fix
	done
	cp -r ./$folder.gffs $start_dir/temp_gffs/.
done



###runing roary
#conda activate Bio
cd $start_dir
roary -f ./roary_output -e -s -n -cd 95 -p 24 -g 100000 ./temp_gffs/*.gff

###prepare the trait table
mkdir temp_traits
dir=$(ls ./input)
for folder in $dir
do
        ls ./input_assemblys/$folder/temp_rename_save > ./temp_traits/$folder.list
	ls ./input_assemblys/$folder/temp_rename_save > ./temp_traits/$folder.csv
done

echo > ./temp_traits/traits.header
echo -e "samples\t"> ./temp_traits/traits.header_pyseer
dir_1=$(ls ./input)
for folder_1 in $dir_1
do
	sed -i "s/$/,$folder_1/" ./temp_traits/traits.header
	sed -i "s/$/\t$folder_1/" ./temp_traits/traits.header_pyseer
	dir_2=$(ls ./input)
	for folder_2 in $dir_2
	do
		if [ "$folder_2" == "$folder_1" ]; then
			sed -i "s/$/,1/" ./temp_traits/$folder_2.csv
		else
			sed -i "s/$/,0/" ./temp_traits/$folder_2.csv
		fi
	done
done
cat ./temp_traits/*.csv > ./temp_traits/temp_all_csv
sed 's/,/\t/g' temp_all_csv > temp_all_tsv
cat ./temp_traits/traits.header ./temp_traits/temp_all_csv > ./temp_traits/final_traits
cat ./temp_traits/traits.header_pyseer ./temp_traits/temp_all_tsv > ./temp_traits/final_traits_pyseer



###runing scoary
#conda activate Bio-scoary
scoary -o ./scoary_output -u --threads 24 -t ./temp_traits/final_traits  -g ./roary_output/gene_presence_absence.csv

                                                 
