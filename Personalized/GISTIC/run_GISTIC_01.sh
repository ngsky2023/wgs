sed 1d 01.cn.txt |sort -k 1,1 -k 2,2 -k 3n,3 |perl remove_overlap_segments.pl - segmentationfile.txt
zcat segmentationfile.txt.gz | cut -f 2-4  | grep -w -v Chr | perl -lane 'print "p$F[0]_$F[1]\t$F[0]\t$F[1]";print "p$F[0]_$F[2]\t$F[0]\t$F[2]"'  |sort -k 2,2 -k 3n,3 | uniq >markersfile1

zcat GermlineHetPon.hg19.vcf.gz|cut -f1,2 |grep -v '#'|sed s/chr//g|awk '{print "p"$1"_"$2"\t"$1"\t"$2}'|sort -k 2,2 -k 3n,3 | uniq >markersfile2
cat markersfile1 markersfile2|sort -k 2,2 -k 3n,3 | uniq >markersfile.txt
mkdir `pwd`/result
#-armpeel 1 染色体臂分离

/share/work1/hanwj4457/softwore/GISTIC2/gistic2 -b `pwd`/result -seg  segmentationfile.txt -refgene /share/work1/hanwj4457/softwore/GISTIC2/refgenefiles/hg19.mat -conf 0.95 -mk markersfile.txt -smallmem 1 -broad 1 -armpeel 1 -savegene 1 -gcm extreme
