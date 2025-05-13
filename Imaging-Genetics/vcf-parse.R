library(VariantAnnotation)
# fl = system.file("extdata", "")
xmat = readVcf("./chr1-snp.recode.vcf", package="VariantAnnotation", "hg19")


##Using vcftools to filter out snp on specific genes
vcftools --vcf ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr1.vcf --bed alz-gene.bed --remove-indels --recode --recode-INFO-all --out chr1-snp
    

vcftools --vcf ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr1.vcf --chr 1 --from-bp 207669472 --to-bp 207815110 --remove-indels --recode --recode-INFO-all --out chr1-snp

###########################################################
##Chromsome 1
###########################################################
library(VariantAnnotation)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genesym <- c("AMPD2")

geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL", columns="ENTREZID")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txbygene = transcriptsBy(txdb, "gene")
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL
chrname = unlist(lapply(strsplit(seqlevels(gnrng), split="chr"), function(x){x[2]}))
seqlevels(gnrng) <-chrname

param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT"))

file = "ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr1.vcf.gz"
xmat = readGT(file, param=param)

###########################################################
##Chromsome 2
###########################################################
library(VariantAnnotation)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genesym <- c("AMPD2")

geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL", columns="ENTREZID")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txbygene = transcriptsBy(txdb, "gene")
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL
chrname = unlist(lapply(strsplit(seqlevels(gnrng), split="chr"), function(x){x[2]}))
seqlevels(gnrng) <-chrname

param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT"))

file = "ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr1.vcf.gz"
xmat = readGT(file, param=param)





hdr <- scanVcfHeader(file)
list.files(system.file("extdata", package="VariantAnnotation"))
f1 <-system.file("extdata","chr22.vcf.gz", package="VariantAnnotation")
vcf<-readVcf(f1,"hg19")
vcf_scan=scanVcf("./chr1-snp.recode.vcf")


GT <- readGT(f1)

## ---------------------------------------------------------------------
## Data subsets with ScanVcfParam
## ---------------------------------------------------------------------

## Subset on genome coordinates:
## 'file' must have a Tabix index
rngs <- GRanges("22", IRanges(c(14370, 1110000), c(17330, 1234600)))
names(rngs) <- c("geneA", "geneB")
param <- ScanVcfParam(which=rngs) 
compressVcf <- bgzip(fl, tempfile())
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)
vcf <- readVcf(tab, "hg19", param)

## When data are subset by range ('which' argument in ScanVcfParam),
## the 'paramRangeID' column provides a map back to the original 
## range in 'param'.
rowData(vcf)[,"paramRangeID"]
vcfWhich(param)

## Subset on samples:







