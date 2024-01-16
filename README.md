WGS Pipeline
===

## Summary
  The WGS Pipeline is used to analyze human whome genome data for somatic SNV/INDEL, somatic CNV, somatic SV, Germline mutation, as well as other downstream analysis such as cancer driver gene prediction and mutational signature.     

## Software Requirements

python>=3.5
perl=5.24.1
java=1.8.0_162
R=3.4.0
fastp=0.13.1
split=8.28
pigz=2.3.4
bwa=0.7.12
samtools=1.4
picard=2.12.1
sambamba=0.6.8
gatk=3.8
gatk=4.0.4
strelka=2.8.4
bcftools=1.4
bedtools=2.19.1
OptiType=1.3.1
msisensor=v0.2
manta=1.2.2
sequenza-utils=2.1.9999b0
annovar=2017Jul16
transvar=2.3.4.20161215
vep=ensembl-tools-release-85
pvacseq=4.0.10

## Download and Install

git clone -b master git@github.com:ngsky2023/wgs.git 

## Usage

### 1.Help

    python3 wgs/WGS.py --help

```
 -h, --help            show this help message and exit
  -i LIST, --list LIST  fastq file list
  -b BAMLIST, --bamList BAMLIST
                        bam file list
  -p PAIR, --pair PAIR  pair file list
  -n NUMBERJOB, --numberJob NUMBERJOB
                        the max synchronic job number
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -a, --noanalysis      no analysis
  -c CONFIG, --config CONFIG
                        config file
  -L BED, --bed BED     One genomic intervals over which to operate
  -s STEP, --step STEP  run step snv,cnv,sv,ana
```

### 2.Example
    
    2.1 python3 wgs/WGS.py -i fastq.txt  -p pair.txt -o  run -c config.txt

    2.2 nohup make -j 50 -f run/.run.job -k -s 2>nohup.out
