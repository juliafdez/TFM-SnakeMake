#BiocManager::install("ChAmp")
#BiocManager::install("annotatr")
#BiocManager::install("BiocParallel")
library(R.utils)
#gunzip('C:/Users/julia/Documents/MASTER/TFM/data/SRR1272798_1_bismark_bt2_pe.bismark.cov.gz')

#bamFile <- BamFile("C:/Users/julia/Downloads/SRR1272791_1_bismark_bt2_pe.bam")
#seqinfo(bamFile)
#bam=scanBam("C:/Users/julia/Downloads/SRR1272791_1_bismark_bt2_pe.bam")
#data(RRBS)
#RRBS
#design

#BiocManager::install("Rsamtools")

#library(Rsamtools)
#
##head(bam)
library(DSS)
require(bsseq)

#dat1 = read.table("C:/Users/julia/Documents/MASTER/TFM/data/SRR1272791_1_bismark_bt2_pe.bismark.cov")
#dat2 = read.table("C:/Users/julia/Documents/MASTER/TFM/data/SRR1272792_1_bismark_bt2_pe.bismark.cov")
#dat3 = read.table("C:/Users/julia/Documents/MASTER/TFM/data/SRR1272793_1_bismark_bt2_pe.bismark.cov")
#dat4 = read.table("C:/Users/julia/Documents/MASTER/TFM/data/SRR1272794_1_bismark_bt2_pe.bismark.cov")

in1 <- gunzip(snakemake@input[[1]])
in2 <- gunzip(snakemake@input[[2]])
in3 <- gunzip(snakemake@input[[3]])
in4 <- gunzip(snakemake@input[[4]])


cat("AQUI: ", snakemake@input[[1]])

el1 <- paste("/home/julia/github/TFM-SnakeMake/snakemake/",in1, sep='')
el2 <- paste("/home/julia/github/TFM-SnakeMake/snakemake/",in2, sep='')
el3 <- paste("/home/julia/github/TFM-SnakeMake/snakemake/",in3, sep='')
el4 <- paste("/home/julia/github/TFM-SnakeMake/snakemake/",in4, sep='')

print(el1)

dat1 = read.table(el1)
dat2 = read.table(el2)
dat3 = read.table(el3)
dat4 = read.table(el4)

print("HE PASADO")

dat1$V3 <- NULL
dat1$V4 <- NULL
dat1$V7 <- dat1$V5 + dat1$V6
dat1$V6 <- NULL

dat2$V3 <- NULL
dat2$V4 <- NULL
dat2$V7 <- dat2$V5 + dat2$V6
dat2$V6 <- NULL

dat3$V3 <- NULL
dat3$V4 <- NULL
dat3$V7 <- dat3$V5 + dat3$V6
dat3$V6 <- NULL

dat4$V3 <- NULL
dat4$V4 <- NULL
dat4$V7 <- dat4$V5 + dat4$V6
dat4$V6 <- NULL

colnames(dat1) <- c("chr","pos","X","N")
colnames(dat2) <- c("chr","pos","X","N")
colnames(dat3) <- c("chr","pos","X","N")
colnames(dat4) <- c("chr","pos","X","N")

heaad <- head(dat1)

#out1 <- snakemake@output[['out1']]
#
#write.table(heaad, file = out1)

BSobj = makeBSseqData( list(dat1, dat2, dat3, dat4),
                       c("C1","C2", "N1", "N2") )[1:1000,]
BSobj

mooth<-BSmooth(BSobj)

stat <- BSmooth.tstat(mooth,group1 = c("C1", "C2"),
                      group2 = c("N1", "N2"),
                      estimate.var = "group2",
                      local.correct = TRUE)

plot(stat)

dmrs0 <- dmrFinder(stat)

pData <- pData(mooth)
pData$col <- rep(c("red", "blue"), each = 2)
pData(mooth) <- pData

plotRegion(mooth, dmrs0[1,], extend = 5000, addRegions = dmrs0)

out3 <- snakemake@output[['out3']]

# Opening the graphical device
pdf(out3)
plot(stat)
plotRegion(mooth, dmrs0[1,], extend = 5000, addRegions = dmrs0)
dev.off()

print("PLOTS OBTENIDOS")

#Esta funci?n realiza los siguientes pasos:
#(1) estimar los niveles medios de metilaci?n para todos los sitios CpG;
#(2) estimar las dispersiones en cada sitio CpG;
#(3) realizar la prueba de Wald.

dmlTest = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
# iteration: 1234
# iteration: 1234
head(dmlTest)

out1 <- snakemake@output[['out1']]

write.table(dmlTest, file = out1)

#Con suavizado:
dmlTest.sm = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"),
                     smoothing=TRUE)


#Con los resultados de la prueba, uno puede llamar a DML usando la callDMLfunci?n.
#Los LMD de resultados se ordenan por importancia.
#
dmls = callDML(dmlTest, p.threshold=0.001)
head(dmls)

out2 <- snakemake@output[['out2']]

write.table(dmls, file = out2)

dmls2 = callDML(dmlTest, delta=0.1, p.threshold=0.001)
head(dmls2)

#dmrs = callDMR(dmlTest, delta=0, p.threshold=1e-5, minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)
#head(dmrs)

#dmrs2 = callDMR(dmlTest, delta=0.1, p.threshold=0.05)
#head(dmrs2)

#showOneDMR(dmrs[1,], BSobj)
