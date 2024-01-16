library('deconstructSigs')
args = commandArgs(trailingOnly=TRUE)

mut.list = as.character(args[1])

sigs.input = mut.to.sigs.input(mut.ref = mut.list, sample.id = "Sample", chr = "CHROM",pos = "POS",ref = "REF",alt = "ALT")
signatures = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, contexts.needed = TRUE, tri.counts.method = 'genome')

weights = t(signatures$weights[,signatures$weights['Sample',]>0])
colnames(weights) = c('weights')
write.table(weights, file='weights.signatures.txt', sep="\t", col.names=T, row.names=T, quote=F)

pdf(file='pie.signatures.pdf')
makePie(signatures)
dev.off()

pdf(file='bar.signatures.pdf')
plotSignatures(signatures)
dev.off()

pdf(file='tumor.signatures.pdf')
plotTumor(signatures[['tumor']])
dev.off()
