
# This script performs the following tasks (xiao.fei@sund.ku.dk):

#     Removes the directories "prokka_ann" and "PEPPAN_working" (if they exist)
#     Creates the directories "prokka_ann", "PEPPAN_working", and "PEPPAN_working/gffs"
#     Renames the contig names in the fasta files in the "input" directory using the "awk" command and saves the output in the "wash_input" directory
#     Runs Prokka on the fasta files in the "wash_input" directory using the "parallel" command and saves the output in the "prokka_ann" directory
#     Copies the GFF files from the "prokka_ann" directory to the "PEPPAN_working/gffs" directory
#     Changes the current directory to "PEPPAN_working"
#     Runs the PEPPAN program on the GFF files in the "PEPPAN_working/gffs" directory using the "-p OUT" option and saves the output in the "OUT.PEPPAN.gff" file
#     Runs the PEPPAN_parser program on the "OUT.PEPPAN.gff" file and saves the output in the "PAR_OUT" directory with a 98% identity threshold and circular genome detection enabled.

rm -r prokka_ann
rm -r PEPPAN_working
rm -r wash_input

start_dir=$(pwd)
mkdir prokka_ann
mkdir PEPPAN_working
mkdir $start_dir/PEPPAN_working/gffs
mkdir wash_input


ls input/ | parallel --verbose -j 40 "awk '/^>/{{print \">contig_\" ++i; next}}{{print}}' input/{} > wash_input/{}"


ls wash_input/ | parallel --verbose -j 40 "prokka --prefix {} -cpus 6 --outdir $start_dir/prokka_ann/{} wash_input/{}"


ls wash_input/ | parallel --verbose -j 40 "cp $start_dir/prokka_ann/{}/{}.gff  $start_dir/PEPPAN_working/gffs/."


cd $start_dir/PEPPAN_working
PEPPAN -p OUT  $start_dir/PEPPAN_working/gffs/*.gff -t 50 
PEPPAN_parser -g $start_dir/PEPPAN_working/OUT.PEPPAN.gff -p PAR_OUT -t -a 98 -c 
