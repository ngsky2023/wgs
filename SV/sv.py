#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class somatic(makeJob):
    def __init__(self, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.rmark = {}
        self.cnmark = {}
        self.cnvmark = {}
        self.xls = {}
        self.segs = {}
    
    def run(self, txt, bam, markd):
        mapping, dat = [], {}
        if isinstance(txt, list):
            mapping = txt
        else:
            mapping = self.readfile(txt)
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)

        for pp in mapping:
            sample = pp[0]
            if sample == '' or sample == '.':
                continue
            normal = ''
            order = ''
            if sample in markd:
                order = markd[sample]
            if len(pp) == 2:
                normal = pp[1]
            else:
                continue

            self.rmark[sample] = {}
            self.xls[sample] = {}
            
            for ch in self.soft['humanChr'].split(' '):
                if sample in markd:
                    if ch in markd[sample]:
                        order = markd[sample][ch]
                if normal in markd:
                    if ch in markd[normal]:
                        order += ' '+markd[normal][ch]

                sp = "%s.%s" % (sample, ch)
                dr = "%ssv/%s/%s/" % (self.outdir, sample, ch)
                prx = dr + sp
                self.rmark[sample][ch] = self.jcmd('Meerkat', {'normal':dat[normal][ch], 'tumor':dat[sample][ch]}, sp, order=order, wd=dr)
                self.xls[sample][ch] = "%ssomaticg.variants" % (dr)

    def mantaGermline(self, bam, markd):
        dat = {}
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)

        mark = {}
        mantadir = {}
        for sample in dat.keys():
            mark[sample] = {}
            bam = dat[sample]
            sp = "%s" % (sample)
            dr = "%ssv/manta/%s/" % (self.outdir, sample)
            prx = dr + sp
            if sample in markd:
                order = markd[sample]

            mark[sample] = self.jcmd('mantaGermline', {'bam':bam, 'dir':dr, 'num':8, 'prx':prx, 'sample':sample}, sp, order=order, wd=dr)
#            mark[sample] = self.jcmd('svanno', {'in':dr+'somaticsv.vcf', 'prx':prx}, sp, order=mark[sample], wd=dr)
            mantadir[sample] = dr
        return mark, mantadir


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
            order = ''
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            mark[sample] = self.jcmd('manta', {'sample':sample, 'Tbam':Tbam, 'Nbam':Nbam, 'dir':dr, 'num':10}, sp, order=order, wd=dr)
            mark[sample] = self.jcmd('SVanno', {'in':dr+'somaticSV.vcf', 'prx':prx}, sp, order=mark[sample], wd=dr)
            mantaDir[sample] = dr
        return mark, mantaDir
    def geneCNV(self, bed, txt, bam, markd):
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
        for pp in mapping:
            sample = pp[0]
            if sample == '' or sample == '.':
                continue

            normal = ''
            if len(pp) == 2:
                normal = pp[1]
            else:
                continue
            self.cnmark[sample] = {}

            Tbam = dat[sample]
            Nbam = dat[normal]
            sp = sample
            dr = "%scnv/%s/" % (self.outdir, sample)
            prx = dr + sp
            order = ''
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            self.cnvmark[sample] = self.jcmd('geneCNV', {'Ts':sample, 'Ns':normal, 'Tbam':Tbam, 'Nbam':Nbam, 'bed':bed}, sample, order=order, wd=dr)
            self.segs[sample] = "%s/Somatic_CNV_Summary.txt" % (dr)


    def cnv(self, bed, txt, bam, markd):
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
        CP = {}
        for pp in mapping:
            sample = pp[0]
            if sample == '' or sample == '.':
                continue

            normal = ''
            if len(pp) == 2:
                normal = pp[1]
            else:
                continue
            self.cnmark[sample] = {}

            chrmark = ''
            for ch in bed.keys():
                if ch == 'chrY' or ch == 'chrM':
                    continue
                cf = bed[ch]
                if cf == ch:
                    cf = ''
                Tbam = dat[sample]
                Nbam = dat[normal]
                sp = "%s.%s" % (sample, ch)
                dr = "%scnv/%s/%s/" % (self.outdir, sample, ch)
                prx = dr + sp
                order = ''
                if sample in markd:
                    order = markd[sample]
                if normal in markd:
                    order += ' '+markd[normal]

                self.cnmark[sample][ch] = self.jcmd('TitanCNA', {'bed':cf, 'Tbam':Tbam, 'Nbam':Nbam, 'chr':ch, 'sample':sample}, sp, order=order, wd=dr)
                chrmark += ' '+self.cnmark[sample][ch]
            #self.cnvmark[sample] = self.jcmd('CNVresult', {'sample':sample, 'dir':"%scnv/%s/" % (self.outdir, sample)}, sample, order=chrmark, wd="%scnv/%s/" % (self.outdir, sample))
            self.cnvmark[sample] = self.jcmd('CNVresult', {'T':dat[sample], 'N':dat[normal], 'ts':sample, 'ns':normal}, sample, order=chrmark, wd="%scnv/%s/" % (self.outdir, sample))
            self.segs[sample] = "%scnv/%s/purple_%s/%s.purple.cnv.somatic.tsv" % (self.outdir, sample, sample, sample)
            CP[sample] = "%scnv/%s/purple_%s/%s" % (self.outdir, sample, sample, sample)
        return CP

                
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



