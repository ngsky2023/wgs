#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class bam(makeJob):
    def __init__(self, bed='', outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.mark = {}
        self.rmdupBam = {}
        self.realignBam = {}
        self.bam = {}
        self.bed = {}
        if bed != '':
            self.cutBed(bed)
        else:
            for ch in self.soft['humanChr'].split(' '):
                self.bed[ch] = ch

    def cutBed(self, bed):
        bedf = {}
        with open(bed, 'r') as fh:
            for line in fh:
                zone = line.strip().split("\t")
                if zone[0] not in bedf:
                    bedf[zone[0]] = []
                bedf[zone[0]].append(zone)
        bedDir = self.mkdir(self.outdir+'bed')
        for chr in bedf:
            self.bed[chr] = bedDir + chr + '.bed'
            with open(self.bed[chr], 'w') as fh:
                for zz in bedf[chr]:
                    print("\t".join(map(str, zz)), file=fh)
    
    def merge(self, dat, markd={}, order=''):
        if order == '':
            order = self.order
        mark = {}
        for sample in dat.keys():
            dr = "%smerge/%s/" % (self.outdir, sample)
            prx = dr + sample
            if sample in markd:
                order = ' '.join(markd[sample].values())
            bbam = dat[sample].values()
            bamf = ' '.join(bbam)
            if len(bbam) > 1:
                mark[sample] = self.jcmd('mergeBam', {'in':bamf, 'prx':prx}, sample, order=order, wd=dr)
                self.bam[sample] = prx + '.bam'
            else:
                mark[sample] = order
                self.bam[sample] = bamf

        return mark

    def rmdup(self, markd, dat='', order=''):
        if order == '':
            order = self.order
        mark = {}
        for sample in dat.keys():
            mark[sample] = {}
            self.rmdupBam[sample] = {}
            for ku in dat[sample].keys():
                mark[sample][ku] = {}
                dr = "%srmdup/%s/%s/" % (self.outdir, sample, ku)
                sp = "%s.%s" % (sample, ku)
                prx = "%s%s" % (dr, sp)
                bam = ' '.join(dat[sample][ku])
                if sample in markd and ku in markd[sample]:
                    order = markd[sample][ku]
                mark[sample][ku] = self.jcmd('rmdup', {'in':bam, 'sample':sp}, sp, order=order, wd=dr)
                self.rmdupBam[sample][ku] = prx + '.rmdup.bam'
        return mark
   
    def QC(self, bam, fq, dataDir, markd={}, order='', bed=''):
        mark = {}
        for sample in bam.keys():
            dr = "%sQC/%s/" % (self.outdir, sample)
            prx = dr + sample
            if sample in markd:
                order = markd[sample]
            if bed == '':
                mark[sample] = self.jcmd('bamQC', {'dataDir': dataDir[sample], 'bam':bam[sample], 'prx':prx}, sample, order=order, wd=dr)
            else:
                mark[sample] = self.jcmd('bedQC', {'dataDir': dataDir[sample], 'bam':bam[sample], 'dir':dr, 'bed':bed}, sample, order=order, wd=dr)
        return mark

    def cut(self, markd, dat='', mapping='', order=''):
        if order == '':
            order = self.order
        if dat == '':
            dat = self.bam
        self.chrBam = {}
        mark = {}
        kk = self.ctype(mapping)
        for sample in dat.keys():
            mark[sample] = {}
            self.chrBam[sample] = {}
            bam = dat[sample]
            if sample in markd:
                order = markd[sample]
            dr = "%schr/%s/" % (self.outdir, sample)
            typec = ''
            if sample in kk:
                typec = kk[sample]

            for ch in self.bed.keys():
                sp = "%s.%s" % (sample, ch)
                prx = dr + sp
                mark[sample][ch] = self.jcmd('cutBam', {'in':bam, 'chr':ch, 'prx':prx, 'type':typec}, sp, order=order, wd=dr)
                self.chrBam[sample][ch] = prx + '.bam'
        return mark

    def ctype(self, mapping):
        k = {}

        if mapping != '':
            for pp in self.readpair(mapping):
                k[pp[0]] = 'T'
                k[pp[1]] = 'N'
        return k

    def readpair(self, txt):
        mapping = []
        with open(txt, 'r') as fh:
            for line in fh:
                block = line.strip().split("\t")
                mapping.append(block)
        return mapping

    def realign(self, markd, dat='', order=''):
        if dat == '':
            dat = self.chrBam
        if order == '':
            order = self.order
        mark = {}
        for sample in dat.keys():
            mark[sample] = {}
            dr = "%srealign/%s/" % (self.outdir, sample)
            self.realignBam[sample] = {}
            for ch in dat[sample].keys():
                cf = self.bed[ch]
                sp = "%s.%s" % (sample, ch)
                prx = dr + sp
                bam = dat[sample][ch]
                if sample in markd:
                    if ch in markd[sample]:
                        order = markd[sample][ch]
                mark[sample][ch] = self.jcmd('realign', {'in':bam, 'prx':prx, 'chr':cf}, sp, order=order, wd=dr)
                self.realignBam[sample][ch] = prx + '.over.bam'
        return mark

    def readfile(self, txt):
        dat = {}
        with open(txt, 'r') as fh:
            for line in fh:
                block = line.strip().split("\t")
                if len(block) == 4:
                        block.append('aa')
                sample, bam, ku, lane, sb = block
                if sample not in dat:
                    dat[sample] = {}
                    dat[sample][ku] = []
                elif ku not in dat[sample]:
                    dat[sample][ku] = []
                dat[sample][ku].append(bam)
        return dat

    def dataTo(self, idat):
        dat = {}
        for sample in idat.keys():
            for ch in idat[sample].keys():
                dat[sample].append(idat[sample][ch])
        return dat





