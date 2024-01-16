tumor=$1
normal=$2
meerkat="/share//work3/wangrr/local/Meerkat/scripts"
mkdir -p tumor normal
cd normal
ln -s $normal normal.bam
ln -s $normal.bai normal.bam.bai
perl $meerkat/pre_process.pl -b normal.bam -s 20 -k 1500 -q 15 -l 0  &
cd -
cd tumor
ln -s $tumor tumor.bam
ln -s $tumor.bai tumor.bam.bai
perl $meerkat/pre_process.pl -b tumor.bam -s 20 -k 1500 -q 15 -l 0 
wait

mv tumor.blacklist.gz tumor.blacklist.gz.bak
ln -s ../normal/normal.blacklist.gz tumor.blacklist.gz
perl $meerkat/meerkat.pl -b tumor.bam -s 20 -d 5 -p 3 -o 1 -m 0 -l 0 -W /share/public/software/bwa-0.7.15/ -F /share/work3/wangrr/DB/hg19/directory/
perl $meerkat/mechanism.pl -b tumor.bam -R /share/work3/wangrr/DB/hg19/rmsk.txt 
cd -

perl $meerkat/somatic_sv.pl -i tumor/tumor.variants -o somatica.variants -F ./normal/ -l 1000 -R /share/work3/wangrr/DB/hg19/rmsk.txt
perl $meerkat/somatic_sv.pl -i somatica.variants -o somaticb.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -n 1 -D 5 -Q 10 -B normal/normal.bam -I normal/normal.isinfo -K normal/normal.blacklist.gz
perl $meerkat/somatic_sv.pl -i somaticb.variants -o somaticc.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -u 1 -Q 10 -B normal/normal.bam -K normal/normal.blacklist.gz
perl $meerkat/somatic_sv.pl -i somaticc.variants -o somaticd.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -f 1 -Q 10 -B normal/normal.bam -K normal/normal.blacklist.gz
perl $meerkat/somatic_sv.pl -i somaticd.variants -o somatice.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -e 1 -D 5 -Q 10 -B tumor/tumor.bam -I tumor/tumor.isinfo -K tumor/tumor.blacklist.gz
perl $meerkat/somatic_sv.pl -i somatice.variants -o somaticf.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -z 1
perl $meerkat/somatic_sv.pl -i somaticf.variants -o somaticg.variants -R /share/work3/wangrr/DB/hg19/rmsk.txt -d 40 -t 20
perl $meerkat/meerkat2vcf.pl -i somaticg.variants -o somaticg.vcf -H head -F /share/work3/wangrr/DB/hg19/directory/
perl $meerkat/fusions.pl -i somaticg.variants -G /share/work3/wangrr/DB/hg19/refGene.txt
perl $meerkat/filter_fusions.pl -i somaticg.fusions -o somatich_fusion -m 4 -s 20000
grep -v -w no_impact somatich_fusion | grep gene-gene | grep -v del | grep -v tandem > fusion.txt
perl $meerkat/discon.pl -d 5 -Q 10 -i somaticg.variants -o somaticg.rp.variants -B tumor/tumor.bam -C tumor/tumor.clusters -I tumor/tumor.isinfo -K tumor/tumor.blacklist.gz -S samtools 

rm -r */*.bam */*.sai */*.sam */*/
