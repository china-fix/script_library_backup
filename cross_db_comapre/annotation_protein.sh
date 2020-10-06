list=$(cat list.txt)
for locus in $list
do
	grep $locus Salmonella_LT2_locus_annotation.csv >> out_annotation.csv 
done
