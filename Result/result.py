#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class over(makeJob):
    def __init__(self, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)

    def run(self, txt, bed, markd):
        mapping = []
        if isinstance(txt, list):
            mapping = txt
        else:
            mapping = self.readfile(txt)
        if bed == '':
            bed = self.soft['bed']

        mark = ''
        sp = 'one'
        dr = "%sresult/" % (self.outdir)
        prx = dr + sp
        order = markd

        mark = self.jcmd('result', {'dir':self.outdir, 'bed':bed}, sp, order=order, wd=dr)
        return mark
       
    def readfile(self, txt):
        mapping = []
        with open(txt, 'r') as fh:
            for line in fh:
                block = line.strip().split("\t")
                mapping.append(block)
        return mapping


