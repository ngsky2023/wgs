#awk '{print $2"\t"$3"\t"$5"\t"$6"\t"$1}' 01.mutation.txt >input.txt
#ln -s /share/work1/wangrr/DB/.bgdata/ /home/hanwj4457/
export PATH="/share/work1/wangrr/local/miniconda3/bin$PATH"
perl 01.get_mutation_for_OncodriveFML.pl
/share/work1/wangrr/local/miniconda3/bin/oncodrivefml --input input.txt --elements gc19_pc.cds.bed --sequencing wgs --configuration oncodrivefml_v2.conf --output gc19_pc_cds --cores 8
