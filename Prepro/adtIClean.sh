for ss in 002/*.adapter.fastq.gz
do
	fq=$(basename $ss | sed 's/adapter/clean/')
	echo "/share/public/software/pigz/2.3.4/pigz -c -d -p 4 $ss | /share/public/software/pigz/2.3.4/pigz -p 4 >> $fq"
done
