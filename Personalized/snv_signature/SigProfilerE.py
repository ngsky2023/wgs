#!/share/work1/wangrr/local/miniconda3/bin/python

from SigProfilerExtractor import sigpro as sig
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='')
parser.add_argument('-o', '--outdir', default='./SigProfilerExtractor', help='')
parser.add_argument('-f', '--format', default='matrix', help='[matrix vcf]')

args = parser.parse_args()


if __name__ == '__main__':
    sig.sigProfilerExtractor(args.format, args.outdir, args.input, refgen="GRCh37", genome_build = 'GRCh37', startProcess=1, endProcess=15, totalIterations=1000, init="alexandrov-lab-custom", cpu=20, mtype = "default",exome = False, penalty=0.05, resample = True, wall= False, gpu=False)
