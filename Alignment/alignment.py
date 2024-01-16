#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class aln(makeJob):
    def __init__(self, txt, markd, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.mark = {}
        self.alnDir = {}
        self.run(txt, markd)
    
    def run(self, dat, markd):
        self.bam = self.outdir+'rawBam.list'
        with open(self.bam, 'w') as bfh:
            for sample in dat.keys():
                self.mark[sample] = {}
                self.alnDir[sample] = "%s%s/" % (self.outdir, sample)
                for ku in dat[sample].keys():
                    self.mark[sample][ku] = ''
                    for lane in dat[sample][ku].keys():
                        order = self.order
                        if sample in markd and ku in markd[sample] and lane in markd[sample][ku]:
                            order = markd[sample][ku][lane]
                        txt = dat[sample][ku][lane]
                        sp = sample+'_'+ku+'_'+lane
                        dr = "%s%s/%s/%s/" % (self.outdir, sample, ku, lane)
                        prx = dr + sp
                        mark = self.jcmd('bwaFor', {'txt':txt, 'prx':prx}, sp, order=order, wd=dr)
                        self.mark[sample][ku] += mark + ' '
                        print("%s\t%s.*.sort.bam\t%s\t%s"%(sample, prx, ku, lane), file=bfh)


