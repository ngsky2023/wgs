#!/bin/bash

task=$1
line=1

if [ $2 ]
then 
	line=$2
else
	line=1
fi

echo '[Start]' `date +"%Y-%m-%d %H:%M:%S"`
cmd=$(head -n $[SGE_TASK_ID * line] $task | tail -n $line)
echo '[CMD start]'
echo "$cmd"
echo '[CMD end]'
echo "$cmd" | sh
excode=$?
if [ $excode == 0 ]
then
	echo '[Success]' `date +"%Y-%m-%d %H:%M:%S"`
else
	echo '[Failed]' `date +"%Y-%m-%d %H:%M:%S"`
fi
exit $excode
