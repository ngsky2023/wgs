source /path/wgs/bashrc
python3 ./wgs/WGS.py -i fastq.txt  -p pair.txt -o `pwd`/run -c config.txt
nohup make -j 50 -f `pwd`/run/.run.job -k -s 2>nohup.out 
