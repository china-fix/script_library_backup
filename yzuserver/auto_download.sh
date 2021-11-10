dir=$(cat list.txt)
for folder in $dir
do
	echo $folder
	ncbi-genome-download -s genbank -F gff,all -t $folder -p 10 -o ./output/$folder bacteria
done	
