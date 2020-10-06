i#!/bin/bash
#remeber use bash only


echo Hello, Xiao! Please input the reference_proteins name! Tak!
read refname

q_files=$(ls *_Q)
for q_file in $q_files
do
	q_file_name=$(echo $q_file | sed -r 's/_HITPROTEINS_Q//')
        list=$(cat $q_file)
        $null > $q_file.Rtab
	for locus in $list
        do
               # grep "STM210"  -xq CDD_HITPROTEINS_Q
		if grep $locus -xq Sal-LT2-VF.txt
		then
			(echo -e "$locus\t$q_file_name\tYes\t$refname")  >>  "$q_file.Rtab"
		else
			(echo -e "$locus\t$q_file_name\tNo\t$refname")  >> "$q_file.Rtab"
		fi
        done
done
rm Interpro_all.Rtab.$refname
cat *.Rtab > Interpro_all.Rtab.$refname
echo done! Please check the file called Interpro_all.Rtab
