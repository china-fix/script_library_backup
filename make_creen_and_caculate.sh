list=$(cat species.list)
for locus in $list
do
	hmmdb=$(grep $locus species-seed_database.list |cut -d , -f2)
	python /home/fix/Desktop/Duo/domain_correlation_hmmer3.py --Q_PROTEINS Query_proteins/"$locus.fa.p" --R_PROTEINS Reference_proteins/"VFDB_setA_NO_$locus.fas" --HMMDB_PATH /home/fix/Downloads/HMMS_DB/"$hmmdb"/"$hmmdb.hmms" --OUTPUT "$locus"_output --CPU 6
	total_num_before=$(grep -c ">" Query_proteins/"$locus.fa.p")
	verified_num_before=$(grep --regexp="$" --count Verified_lists/"$locus.vflist.before")
	total_num_after=$(grep --regexp="$" --count "$locus"_output/HITPROTEINS/_HITPROTEINS_Q)
	list_1=$(cat "$locus"_output/HITPROTEINS/_HITPROTEINS_Q)
	for locus_1 in $list_1
	do
		grep $locus_1 Verified_lists/"$locus.vflist.before" >> "$locus"_output/HITPROTEINS/experiment_confirmed_after
	done
	verified_num_after=$(grep --regexp="$" --count "$locus"_output/HITPROTEINS/experiment_confirmed_after)
	echo "total_number_before: $total_num_before" >> "report.$locus.txt"
	echo "verified_number_before: $verified_num_before" >> "report.$locus.txt"
	echo "total_number_after: $total_num_after" >> "report.$locus.txt"
	echo "verified_number_after: $verified_num_after" >> "report.$locus.txt"
done
