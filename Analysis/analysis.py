#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class analy(makeJob):
    def __init__(self, outdir='./', beforeJob='', config='', bed=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.bed = bed

    def msi(self, txt, bam, markd):
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
        bed = self.bed
        if bed != '':
            bed = '-e ' + bed
        for pp in mapping:
            sample = pp[0]
            if sample == '' or sample == '.':
                continue
            normal = ''
            if len(pp) == 2:
                normal = pp[1]
            else:
                continue
            Tbam = dat[sample]
            Nbam = dat[normal]
            sp = "%s" % (sample)
            dr = "%smsi/%s/" % (self.outdir, sample)
            prx = dr + sp
            order = ''
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            mark[sample] = self.jcmd('msi', {'tumor':Tbam, 'normal':Nbam, 'prx':prx, 'bed':bed}, sp, order=order, wd=dr)
        return mark
       
    def purity(self, txt, bam, markd):
        mapping, dat = [], {}
        if isinstance(txt, list):
            mapping = txt
        else:
            mapping = self.readfile(txt)
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)
        bed = self.soft['bed']
        if self.bed != '':
            bed = self.bed

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
            Tbam = dat[sample]
            Nbam = dat[normal]
            sp = "%s" % (sample)
            dr = "%spurity/%s/" % (self.outdir, sample)
            prx = dr + sp
            order = ''
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            mark[sample] = self.jcmd('purity', {'tumor':Tbam, 'normal':Nbam, 'sample':sample, 'bed':bed}, sp, order=order, wd=dr)
            CP[sample] = prx
        return mark, CP

    def ccf(self, snv, cnv, purity, markd):
        mark = {}
        for sample in snv.keys():
            sp = "%s" % (sample)
            dr = "%sCCF/%s/" % (self.outdir, sample)
            prx = dr + sp
            order = ''
            if sample in markd:
                order = markd[sample]
            purityf = purity[sample] + '.purple.purity.tsv'
            segf = purity[sample] + '_segments.txt'
            if self.bed == '':
                segf = cnv[sample]
            mark[sample] = self.jcmd('CCF', {'sample':sample, 'snv':snv[sample], 'cnv':segf, 'purity':purityf, 'prx':prx}, sp, order=order, wd=dr)
        return mark

    def antigen(self, snv, vcf, hlaDir, markd):
        mark = {}
        for sample in snv.keys():
            sp = "%s" % (sample)
            dr = "%santigen/%s/" % (self.outdir, sample)
            prx = dr + sp
            order = ''
            if sample in markd:
                order = markd[sample]
            mark[sample] = self.jcmd('pvacseq', {'ann':snv[sample], 'in':vcf[sample], 'sample':sample, 'hladir':hlaDir[sample]}, sp, order=order, wd=dr)
        return mark


    def hla(self, bam, markd):
        dat = {}
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)

        mark = {}
        drr = {}
        for sample in dat.keys():
            order = ''
            if sample in markd:
                order = markd[sample]
            sp = "%s" % (sample)
            dr = "%sHLA/%s/" % (self.outdir, sample)
            mark[sample] = self.jcmd('hlaTyping', {'bam':dat[sample], 'dir':dr}, sp, order=order, wd=dr)
            drr[sample] = dr
        return mark, drr

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



