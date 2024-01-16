#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys, argparse, subprocess, shutil
#import data, alignment, bamDeal, snv, sv, analysis, result
from SnvIndel import snv
from SV import sv
from Analysis import analysis
from Prepro import data
from Alignment import alignment,bamDeal
from Result import result

parser = argparse.ArgumentParser(description='genome reseq analysis pipeline')

ingroup = parser.add_mutually_exclusive_group(required=True)
ingroup.add_argument('-i', '--list', help='fastq file list')
ingroup.add_argument('-b', '--bamList', help='bam file list')

parser.add_argument('-p', '--pair', help='pair file list')
parser.add_argument('-n', '--numberJob', help='the max synchronic job number', default='50')
parser.add_argument('-o', '--outdir', help='output directory', default='./')
parser.add_argument('-a', '--noanalysis', help='no analysis', action="store_true")
#parser.add_argument('-r', '--norun', help='no run the job, only generate shell script and directory', action="store_true")
parser.add_argument('-c', '--config', help='config file', default='')
parser.add_argument('-L', '--bed', help='One genomic intervals over which to operate', default='')
parser.add_argument('-s', '--step', help='run step snv,cnv,sv,ana', default='snv,cnv,sv,ana')
args = parser.parse_args()

if args.bed != '':
    if not os.path.exists(args.bed):
        print("%s is not exists." % (args.bed), file=sys.stderr)
        sys.exit(1)
    else:
        args.bed = os.path.realpath(args.bed)

def mergeMark(mm):
    mo = {}
    for mi in mm:
        for key in mi.keys():
            if key in mo:
                mo[key] += ' ' + mi[key]
            else:
                mo[key] = mi[key]
    return mo
def mergeMarkL(omark):
    end = ''
    for emk in omark:
        if isinstance(emk, str):
            end += ' ' + emk
        else:
             for k in emk.keys():
                if isinstance(emk[k], dict):
                    for ch in emk[k].keys():
                        end += ' ' + emk[k][ch]
                else:
                    end += ' ' + emk[k]
           
    return end



args.outdir = os.path.abspath(args.outdir)
obase = os.path.basename(args.outdir)
if args.outdir[-1] != '/':
    args.outdir += '/'

mft = '.'+obase+'.job'
makefile = open(mft, 'w')
stdout = sys.stdout
sys.stdout = makefile
print("all : end\n")


if args.list:
    dat = data.Preprocess(outdir=args.outdir+'data', config=args.config)
    mark, fq = dat.cutData(args.list)
    aln = alignment.aln(fq, mark, outdir=args.outdir+'alignment', config=args.config)

rmdupBam = {}
realignBam = {}
if args.bamList:
    with open(args.bamList, 'r') as bamfh:
        for ll in bamfh:
            bb = ll.strip().split("\t")
            if len(bb) == 2:
                rmdupBam[bb[0]] = bb[1]
            elif len(bb) == 3:
                realignBam[bb[0]][bb[1]] = bb[2]

deal = bamDeal.bam(bed = args.bed, outdir=args.outdir+'bam', config=args.config)

qcMark = ''
bamMark = {}
if len(realignBam) > 0:
    deal.realignBam = realignBam
    deal.bam = rmdupBam
elif len(rmdupBam) > 0:
    deal.bam = rmdupBam
    chrMark = deal.cut('')
    mark = deal.realign(chrMark)
elif args.list:
    mark = deal.rmdup(aln.mark, deal.readfile(aln.bam))
    bamMark = deal.merge(deal.rmdupBam, mark)
    qcMark = deal.QC(deal.bam, dat.r1, dat.dataDir, bamMark, bed=args.bed)
    if not args.noanalysis:
        chrMark = deal.cut(bamMark, deal.bam, args.pair)
        mark = deal.realign(chrMark)
    else:
        mark = bamMark
else:
    sys.exit(1)


omark = [mark, qcMark]
if not args.noanalysis:
    if args.pair:
        msa = analysis.analy(outdir=args.outdir+'analysis', config=args.config, bed=args.bed)
        msa.purity(args.pair, deal.bam, bamMark)
        msiMark = msa.msi(args.pair, deal.bam, bamMark)
        hlaMark, hlaDir = msa.hla(deal.bam, bamMark)
        sis = snv.somatic(args.pair, deal.realignBam, deal.bed, markd=mark, outdir=args.outdir+'snv-indel', config=args.config)
        antMark = msa.antigen(sis.result, sis.mvcf, hlaDir, mergeMark([sis.mark, hlaMark]))
    
        svs = sv.somatic(args.outdir, config=args.config)
        #svs.run(args.pair, deal.chrBam, chrMark)
        if args.bed == '':
            CP = svs.cnv(deal.bed, args.pair, deal.bam, bamMark)
        else:
            svs.geneCNV(args.bed, args.pair, deal.bam, bamMark)
        mantaMark, mantaDir = svs.manta(args.pair, deal.bam, bamMark)
        strelkaMark = sis.strelka(args.pair, deal.bam, mantaMark, mantaDir)
        ccfMark = msa.ccf(sis.result, svs.segs, CP, mergeMark([sis.mark, svs.cnvmark]))
        #omark = [qcMark, sis.mark, svs.rmark, svs.cnvmark, mantaMark, strelkaMark]
        #omark = [qcMark, sis.mark, svs.cnvmark, mantaMark, ccfMark, msiMark, hlaMark, antMark]
        germMark = sis.acmg(args.pair, args.outdir+"bam/realign", markd=mark)
        omark = [qcMark]
        if 'snv' in args.step:
            omark.append(sis.mark)
        if 'sv' in args.step:
            omark.append(mantaMark)
        if 'cnv' in args.step:
            omark.append(svs.cnvmark)
        if 'ana' in args.step:
            omark.extend([ccfMark, msiMark, hlaMark, antMark])
        if 'germ' in args.step:
            omark.append(germMark)
        if 'snv' in args.step and args.bed != '':
            re = result.over(args.outdir, config=args.config)
            resmark = re.run(args.pair, args.bed, mergeMarkL([sis.mark, svs.cnvmark]))
            omark.append(resmark)

    else:
        sis = snv.somatic(False, deal.realignBam, deal.bed, markd=mark, outdir=args.outdir+'snv-indel', config=args.config)
        svs = sv.somatic(args.outdir, config=args.config)
        mantaMark, mantaDir = svs.mantaGermline(deal.bam, bamMark)
        omark = [qcMark, sis.mark, mantaMark]

        

end = mergeMarkL(omark)

print("end : "+end+"\n")
sys.stdout = stdout
makefile.close 

if os.getcwd()+'/' != args.outdir:
    if os.path.exists(args.outdir+mft):
        os.remove(args.outdir+mft)
    shutil.move(mft, args.outdir)

CMD="make -j %s -f %s -k -s" % (args.numberJob, args.outdir+mft)
#if args.norun:
print(CMD)
#else:
#    stdout = open(args.outdir+obase+'.log', 'w')
#    stderr = open(args.outdir+obase+'.err', 'w')
#    jobo = subprocess.Popen(CMD.split(' '), stdout=stdout, stderr=stderr, cwd=args.outdir)
#    jobo.wait()
#    stdout.close
#    stderr.close



