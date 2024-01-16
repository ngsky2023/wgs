#!/share/public/software/R-3.3.3/bin/R
library('getopt')
.libPaths(c( .libPaths(),"/share/public/software/lib_ssinfo/R/R-3.3.3"))
spec = matrix(c(
        'help' , 'h', 0, "logical",
        'n_sigs' , 'n', 1, "numeric",
	'lst' ,'l' ,1, "character",
	'ref' ,'r' ,1 , "character",
        'outdir' , 'o', 1, "character"
        ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
        cat(getopt(spec, usage=TRUE));
        cat("Usage example: \n")
        cat(" 
R Signatures.R  --n_sigs|-n 2 --lst|-l 1.vcf,2.vcf -ref|-r /home/zhanghk/software/CancerPipe/hg19/hg19.fa -outdir|-o ./
Options: 
--help          -h      NULL            get this help
--n_sigs        -n      numeric		define the mutational signatures. it should be lower than sample numbers. [forced]
--lst		-l	character	the lst contain mutect out result.[forced]
--ref		-r	character	the ref.fa path [forced]
--outdir	-o      character       the prefix for output graph [forced]
#-----------------------------------------------------------------
\n")
        q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$n_sigs) ){ print_usage(spec) }
if ( is.null(opt$lst) )   { print_usage(spec) }
if ( is.null(opt$ref) )   { print_usage(spec) }
if ( is.null(opt$outdir)) { print_usage(spec) }
#################################################################################
library("SomaticSignatures")
library("NMF")
library("ggplot2")
#################################################################################

split.vr <- unlist(strsplit(opt$lst,","))

vr<-NULL
for (i in 1:length(split.vr)){
	vr.temp <- readVcfAsVRanges(split.vr[i],opt$ref)
	if (is.null(vr)){vr<-vr.temp}
	else{vr<-c(vr,vr.temp)}
}

fa_A = FaFile(opt$ref)
sca_motifs = mutationContext(vr,fa_A)
sca_mm = motifMatrix(sca_motifs, group = "sampleNames", normalize = TRUE)
setwd(opt$outdir)

pdf(file="Mutation_spectrum_over_sample.pdf")
plotMutationSpectrum(sca_motifs, "sampleNames")
dev.off()

pdf(file="nmf-heatmap.pdf")
sca_mm=sca_mm[sca_mm[,'NORMAL']>0,]
sigs_nmf = identifySignatures(sca_mm, opt$n_sigs, nmfDecomposition, method="lee")
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()
write.table(sigs_nmf@signatures,"signatures.xls")

pdf(file="nmf-signatures.pdf")
plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")
dev.off()

pdf(file="sample-signatures_contribution.pdf")
#plotSamples(sigs_nmf)
plotSamples(sigs_nmf, normalize=TRUE)
dev.off()

png(file="Mutation_spectrum_over_sample.png")
plotMutationSpectrum(sca_motifs, "sampleNames")
dev.off()

png(file="nmf-heatmap.png")
sigs_nmf = identifySignatures(sca_mm, opt$n_sigs, nmfDecomposition, method="lee")
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()

png(file="nmf-signatures.png")
plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")
dev.off()

png(file="sample-signatures_contribution.png")
#plotSamples(sigs_nmf)
plotSamples(sigs_nmf,normalize=TRUE) ## to fix bar not 100%
dev.off()

png(file="assessNumberSignatures.png")
n_sigs = 2:8
gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)
plotNumberSignatures(gof_nmf)
dev.off()

