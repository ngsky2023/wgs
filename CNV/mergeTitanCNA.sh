dir=$1
sample=$2
cat $dir/chr*/$sample.chr*.seg | head -1 > $sample.seg
cat $dir/chr*/$sample.chr*.seg | grep -v sample >> $sample.seg
cat $dir/chr*/$sample.chr*.segs.txt | head -1 > $sample.segs.txt
cat $dir/chr*/$sample.chr*.segs.txt | grep -v Sample >> $sample.segs.txt
