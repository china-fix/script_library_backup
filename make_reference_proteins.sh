list=$(cat species.list)
for locus in $list
do
	seqkit grep -n -v -r -p "\[$locus " VFDB_setA_pro.fas > ./Reference_proteins/"VFDB_setA_NO_$locus.fas"
done
