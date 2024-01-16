#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class detect(makeJob):
    def __init__(self, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.rmark = {}
        self.cnmark = {}
        self.cnvmark = {}
        self.xls = {}
        self.segs = {}
    
    def manta(self, txt, bam, markd):
        mapping, dat = [], {}
        if isinstance(txt, list):
            mapping = txt
        else:
            mapping = self.readfile(txt)
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)

        mark = {}
        mantaDir = {}
        for pp in mapping:
            sample = pp[0]
            if sample == '' or sample == '.':
                continue

            normal = ''
            if len(pp) == 2:
                normal = pp[1]
            else:
                continue
            mark[sample] = {}
            Tbam = dat[sample]
            Nbam = dat[normal]
            sp = "%s" % (sample)
            dr = "%ssv/manta/%s/" % (self.outdir, sample)
            prx = dr + sp
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            mark[sample] = self.jcmd('manta', {'Tbam':Tbam, 'Nbam':Nbam, 'dir':dr, 'num':10}, sp, order=order, wd=dr)
            mark[sample] = self.jcmd('SVanno', {'in':dr+'somaticSV.vcf', 'prx':prx}, sp, order=mark[sample], wd=dr)
            mantaDir[sample] = dr
        return mark, mantaDir

    def readfile(self, txt):
        mapping = []
        with open(txt, 'r') as fh:
            for line in fh:
                block = line.strip().split("\t")
                mapping.append(block)
        return mapping

    def readbam(self, txt):
        bam = {}
        with open(txt, 'r') as fh:
            for line in fh:
                sample, ch, f = line.strip().split("\t")
                if sample not in bam:
                    bam[sample] = {}
                bam[sample][ch] = f
        return bam



