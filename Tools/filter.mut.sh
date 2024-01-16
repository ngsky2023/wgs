#!/bin/bash -
#===============================================================================
#
#          FILE: filter.sh
#
#         USAGE: ./filter.sh
#
#   DESCRIPTION: 
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Wang Ruiru (wangrr), wangruiru@berrygenomics.com
#  ORGANIZATION: Berry Genomics
#       CREATED: 12/12/2019 05:15:26 PM
#      REVISION:  ---
#===============================================================================

set -o nounset                                  # Treat unset variables as an error

xls=$1
bam=$2

sed '1d' $xls | cut -f 2-6 | grep ^chr | sort -u > pos

echo -n '' > cmd.sh
cat $bam | while read bb
do
	sample=$(basename $bb | sed 's/\..*$//')
	echo "/share/work1/wangrr/local/simple/bin/mutationInfo -F 1024 -pos pos -file $bb > $sample.info.txt" >> cmd.sh
done
n=$(wc -l cmd.sh | cut -d ' ' -f 1)
qsub -cwd -l vf=1g -q oncrd.q -sync y -t 1-$n -tc $n -N POS /share/work1/wangrr/local/simple/bin/runtask.sh cmd.sh
perl /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/newadd/count/m.pl *.info.txt > m.txt

perl /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/newadd/overlap/filter.pl -i $xls > filter.xls
