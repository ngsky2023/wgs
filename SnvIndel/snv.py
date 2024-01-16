#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
from pipeline import makeJob

class somatic(makeJob):
    def __init__(self, txt, bam, bed, markd={}, outdir='./', beforeJob='', config=''):
        makeJob.__init__(self, beforeJob, outdir, config)
        self.mark = {}
        self.snvMark = {}
        self.annMark = {}
        self.annot = {}
        self.annotDir = {}
        self.ipt = {}
        self.vcf = {}
        self.mvcf = {}
        self.vcfDir = {}
        mapping, dat = [], {}
        if isinstance(bam, dict):
            dat = bam
        else:
            dat = self.readbam(bam)
        if txt:
            if isinstance(txt, list):
                mapping = txt
            else:
                mapping = self.readfile(txt)
            self.run(mapping, dat, bed, markd)
        else:
            self.germline(dat, markd)

        self.annotation(txt)
        self.result = {}
        self.merge(txt)

    def germline(self, dat, markd={}):
        order = self.order
        for sample in dat.keys():
            if sample not in self.mark:
                self.mark[sample] = {}
            if sample not in self.vcf:
                self.vcf[sample] = {}
            if sample not in self.snvMark:
                self.snvMark[sample] = ''
            self.vcfDir[sample] = dr = "%scall/%s/" % (self.outdir, sample)

            for ch in dat[sample].keys():
                if sample in markd:
                    if ch in markd[sample]:
                        order = markd[sample][ch]
                bam = dat[sample][ch]
                sp = "%s.%s" % (sample, ch)
                dr = "%scall/%s/%s/" % (self.outdir, sample, ch)
                prx = dr + sp
                self.mark[sample][ch] = self.jcmd('gatk4HaplotypeCaller', {'bam':bam, 'prx':prx, 'chr':ch}, sp, order=order, wd=dr)
                #self.mark[sample][ch] = self.jcmd('HaplotypeCaller', {'bam':bam, 'prx':prx, 'chr':ch}, sp, order=order, wd=dr)
                self.snvMark[sample] += " " + self.mark[sample][ch]
                self.vcf[sample][ch] = prx+".raw.vcf"

    def acmg(self, mapping, datdir, markd={}):
        mapping = self.readfile(mapping)

        bed = self.soft['gbed'];
        rt = self.mkdir("%sgerm/" % (self.outdir))
        cc = {}
        with open(bed, 'r') as bh:
            for ll in bh:
                llb = ll.strip().split("\t")
                cc[llb[0]] = 1

        bh = open(rt + 'bam.txt', 'w')
        sh = open(rt + 'sample', 'w')
        mark = ''
        for pp in mapping:
            sample = pp[1]
            if sample == '' or sample == '.':
                continue
            for ch in cc.keys():
                print("%s\t%s/%s/%s.%s.over.bam" % (sample, datdir, sample, sample, ch), file=bh)
                mark += ' ' + markd[sample][ch]
            print(".\t%s" % (sample), file=sh)
        bh.close()
        sh.close()
        mark = self.jcmd('acmg', {'txt':rt + 'bam.txt', 'bed':bed}, 'germ', order=mark, wd=rt)
        return mark
    
    def run(self, mapping, dat, bed, markd={}):
        order = self.order
        for pp in mapping:
            #-I {tumor} -tumor {ts} -I {normal} -normal {ns}
            sample = pp[0]
            if sample == '' or sample == '.':
                continue
            normal = ''
            if sample not in self.mark:
                self.mark[sample] = {}
            if sample not in self.vcf:
                self.vcf[sample] = {}
            if sample not in self.snvMark:
                self.snvMark[sample] = ''

            self.vcfDir[sample] = dr = "%scall/%s/" % (self.outdir, sample)
            for ch in dat[sample].keys():
                cf = bed[ch]
                inapt = ''
                if sample in markd:
                    if ch in markd[sample]:
                        order = markd[sample][ch]
                if len(pp) == 2:
                    normal = pp[1]
                    inapt = "-I %s -tumor %s -I %s -normal %s" % (dat[sample][ch], sample, dat[normal][ch], normal)
                    if normal in markd:
                        if ch in markd[normal]:
                            order += ' ' + markd[normal][ch]
                else:
                    inapt = "-I %s -tumor %s" % (dat[sample][ch], sample)
                sp = "%s.%s" % (sample, ch)
                dr = "%scall/%s/%s/" % (self.outdir, sample, ch)
                prx = dr + sp
                self.mark[sample][ch] = self.jcmd('Mutect2', {'in':inapt, 'prx':prx, 'chr':cf, 'bam':"%s %s"%(dat[sample][ch], dat[normal][ch])}, sp, order=order, wd=dr)
                #self.mark[sample][ch] = self.jcmd('MuTect2', {'tumor':dat[sample][ch], 'normal':dat[normal][ch], 'prx':prx, 'chr':cf}, sp, order=order, wd=dr)

                self.snvMark[sample] += " " + self.mark[sample][ch]
                self.vcf[sample][ch] = prx+".filter.vcf"

    def strelka(self, txt, bam, markd, mantaDir):
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
            mark[sample] = ''
            Tbam = dat[sample]
            Nbam = dat[normal]
            sp = "%s" % (sample)
            dr = "%sstrelka/%s/" % (self.outdir, sample)
            prx = dr + sp
            if sample in markd:
                order = markd[sample]
            if normal in markd:
                order += ' '+markd[normal]

            markk = self.jcmd('strelka', {'sample': sample, 'Tbam':Tbam, 'Nbam':Nbam, 'dir':dr, 'num':20, 'mantaDir':mantaDir[sample]}, sp, order=order, wd=dr)
            mark[sample] = self.jcmd('annovar', {'in':dr+'somatic.snvs.vcf', 'prx':prx+'.snps', 'sample':sample, 'sapl':'TUMOR'}, sp, order=markk, wd=dr)
            mark[sample] += ' ' + self.jcmd('annovar', {'in':dr+'somatic.indels.vcf', 'prx':prx+'.indels', 'sample':sample, 'sapl':'TUMOR'}, sp+'_indel', order=markk, wd=dr)
        return mark

    def annotation(self, txt):
        for sample in self.vcf.keys():
            if sample not in self.mark:
                self.mark[sample] = {}
            if sample not in self.annot:
                self.annot[sample] = {}
            if sample not in self.annMark:
                self.annMark[sample] = ''

            self.annotDir[sample] = "%sannotation/%s/" % (self.outdir, sample)
            for ch in self.vcf[sample].keys():
                sp = "%s.%s" % (sample, ch)
                dr = "%sannotation/%s/%s/" % (self.outdir, sample, ch)
                prx = dr + sp
                order = self.mark[sample][ch]
                sapl = sample
                if txt:
                    sapl = 'TUMOR'
                self.mark[sample][ch] = self.jcmd('annovar', {'in':self.vcf[sample][ch], 'prx':prx, 'sample':sample, 'sapl':sapl}, sp, order=order, wd=dr)
                self.annMark[sample] += ' ' + self.mark[sample][ch]
                self.mark[sample][ch] += ' ' + self.jcmd('wgsa', {'in':self.vcf[sample][ch], 'prx':prx}, sp, order=order, wd=dr)
                self.annot[sample][ch] = prx + '.snp.gz'

    def interpretation(self):
        for sample in self.vcf.keys():
            vcf = ' '.join(self.vcf[sample].values())
            order = ' '.join(self.mark[sample].values())
            annot = ' '.join(self.annot[sample].values())
            dr = "%sinterpretation/%s/" % (self.outdir, sample)
            prx = dr + sample
            self.mark[sample] = self.jcmd('pcgr', {'in':vcf, 'prx':prx, 'annot':annot, 'sample':sample, 'outdir':dr}, sample, order=order, wd=dr)
            self.ipt[sample] = prx + '.pcgr.snvs_indels.tiers.tsv'

    def merge(self, txt):
        for sample in self.vcf.keys():
            dr = "%sresult/%s/" % (self.outdir, sample)
            annDir = "%sannotation/%s/" % (self.outdir, sample)
            prx = dr + sample
            if txt:
                order = self.annMark[sample]
                #self.mark[sample] = self.jcmd('mergeSNP', {'vdir':self.vcfDir[sample], 'adir':self.annotDir[sample], 'sample':sample}, sample, order=order, wd=dr) + ' ' + ' '.join(self.mark[sample].values())
                self.mark[sample] = self.jcmd('mergeSNP', {"sample": sample, 'vdir':self.vcfDir[sample], 'adir':self.annotDir[sample], 'sample':sample}, sample, order=order, wd=dr)

            else:
                order = self.snvMark[sample]
                self.mark[sample] = self.jcmd('merge4Germline', {'vdir':self.vcfDir[sample], 'adir':self.annotDir[sample], 'sample':sample, 'prx':prx, 'annDir':annDir}, sample, order=order, wd=dr) + ' ' + ' '.join(self.mark[sample].values())
            self.result[sample] = dr + sample + '.final.xls'
            self.mvcf[sample] = dr + sample + '.vcf'

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



if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-b', '--bam', help='bam file list', required=True)
    parser.add_argument('-p', '--pair', help='pair file list', required=True)
    parser.add_argument('-o', '--outdir', help='output dir', default='./')
    parser.add_argument('-c', '--config', help='config file', default='')
    args = parser.parse_args()
    makefile = open('snv.job', 'w')
    sys.stdout = makefile
    somatic(args.pair, args.bam, outdir=args.outdir)
    makefile.close 



