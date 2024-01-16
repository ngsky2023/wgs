#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class Preprocess(makeJob):
    def __init__(self, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.r1 = {}
        self.dataDir = {}
    
    def cutData(self, txt):
        dat = []
        mark = {}
        out = {}
        if isinstance(txt, list):
            dat = txt
        else:
            dat = self.readfile(txt)

        fq = {}
        for dd in dat:
            sample, r1, r2, ku, lane = dd
            if sample not in self.r1:
                self.r1[sample] = []
            self.r1[sample].append(r1)

            if sample not in mark:
                mark[sample] = {}
            if sample not in out:
                out[sample] = {}
            if ku not in mark[sample]:
                mark[sample][ku] = {}
            if ku not in out[sample]:
                out[sample][ku] = {}
            if sample not in fq:
                fq[sample] = {}
            if ku not in fq[sample]:
                fq[sample][ku] = {}
            if sample not in self.dataDir:
                self.dataDir[sample] = "%s%s/" % (self.outdir, sample)

            sp = "%s_%s_%s" % (sample, ku, lane)
            dr = "%s%s/%s/%s/" % (self.outdir, sample, ku, lane)
            prx = dr + sp

            mark[sample][ku][lane] = self.jcmd('fastp', {'r1':r1, 'r2':r2, 'prx':prx}, sp, order=self.order, wd=dr)
            r1 = dr + '*.' + sp + '.R1.fq.gz'
            r2 = dr + '*.' + sp + '.R2.fq.gz'
            out[sample][ku][lane] = [dr, r1, r2, sample, ku, lane]
            with open(dr+'rawData.txt', 'w') as of:
                print("\t".join((sample, r1, r2, ku, lane)), file=of)
            fq[sample][ku][lane] = dr+'rawData.txt'

        return mark, fq

    def readfile(self, txt):
        dat = []
        with open(txt, 'r') as fh:
            for line in fh:
                sample, r1, r2, ku, lane = line.strip().split("\t")
                dat.append([sample, r1, r2, ku, lane])
        return dat

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-l', '--list', help='fastq file list', required=True)
    parser.add_argument('-o', '--outdir', help='output dir', default='./')
    parser.add_argument('-c', '--config', help='config file', default='')
    args = parser.parse_args()
    makefile = open('data.job', 'w')
    sys.stdout = makefile
    Preprocess(args.list, args.outdir)
    makefile.close 



