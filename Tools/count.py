#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class cnv(makeJob):
    def __init__(self, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.mark = {}
    
    def run(self, bam, markd):
        dat = {}
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)

        mark = {}
        mantaDir = {}
        for sample in dat.keys():
            mark[sample] = {}
            bam = dat[sample]
            sp = "%s" % (sample)
            dr = "%s%s/" % (self.outdir, sample)
            prx = dr + sp
            if sample in markd:
                order = markd[sample]

            mark[sample] = self.jcmd('count', {'bam':bam, 'sample':sample, 'dir':dr}, sp, order=order, wd=dr)
        self.mark = mark

    def readbam(self, txt):
        bam = {}
        with open(txt, 'r') as fh:
            for line in fh:
                sample, ch, f = line.strip().split("\t")
                if sample not in bam:
                    bam[sample] = {}
                bam[sample][ch] = f
        return bam



