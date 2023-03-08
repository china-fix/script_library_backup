

rm -r prokka_ann
rm -r PIRATE_working
rm -r wash_input

start_dir=$(pwd)
mkdir prokka_ann
mkdir PIRATE_working
mkdir $start_dir/PIRATE_working/gffs
mkdir wash_input


ls input/ | parallel --verbose -j 40 "awk '/^>/{{print \">contig_\" ++i; next}}{{print}}' input/{} > wash_input/{}"


ls wash_input/ | parallel --verbose -j 40 "prokka --prefix {} --cpus 6  --outdir  $start_dir/prokka_ann/{} wash_input/{}"


ls wash_input/ | parallel --verbose -j 40 "cp $start_dir/prokka_ann/{}/{}.gff  $start_dir/PIRATE_working/gffs/."


cd $start_dir/PIRATE_working
# PEPPAN -p OUT  $start_dir/PIRATE_working/gffs/*.gff -t 50 
PIRATE -i $start_dir/PIRATE_working/gffs -o $start_dir/PIRATE_working -r -t 40 
# PEPPAN_parser -g $start_dir/PIRATE_working/OUT.PEPPAN.gff -p PAR_OUT -t -a 98 -c 

echo "hi xiao job done!"