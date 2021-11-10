start_dir=$(pwd)
mkdir assembly_collect_dir
mkdir assembly_collect_file
tax_ids=$(ls ./output)
for tax_id in $tax_ids
do
	echo $tax_id
	mkdir ./assembly_collect_dir/$tax_id
	assembly_folders=$(ls $start_dir/output/$tax_id/genbank/bacteria/)
	for assembly_folder in $assembly_folders
	do
		cp $start_dir/output/$tax_id/genbank/bacteria/$assembly_folder/*_genomic.fna.gz ./assembly_collect_dir/$tax_id/.
		rm ./assembly_collect_dir/$tax_id/*_from_genomic.fna.gz
	done
done

for tax_id in $tax_ids
do
	assembly_files=$(ls ./assembly_collect_dir/$tax_id)
	for assembly_file in $assembly_files
	do
		cp ./assembly_collect_dir/$tax_id/$assembly_file ./assembly_collect_file/"$tax_id"_"$assembly_file"
		gunzip ./assembly_collect_file/"$tax_id"_"$assembly_file"
	done
done
