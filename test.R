getwd()
q()
q()
org=data.frame(row.names = c('Mus musculus','Rattus norvegicus','Homo sapiens'))
#org$code=c('mmu','rno','hsa')
org$db=c("org.Mm.eg.db","org.Rn.eg.db","org.Hs.eg.db")
org=as.matrix(org)
org
library(shiny)
library(shinycssloaders)
runApp('withinInTAD_test.R')
runApp('withinInTAD_test.R')
runApp('Desktop/work/pipelines/inTAD/inTAD-shiny/test_5.R')
runApp('withinInTAD_test.R')
library(shiny); runApp('Desktop/work/pipelines/inTAD/inTAD-shiny/test_6_process.R')
setwd("/Users/venu/Desktop/work/pipelines/inTAD/inTAD-shiny/")
dir()
library(shiny); runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
runApp('test_6_process.R')
source("http://bioconductor.org/biocLite.R")
biocLite("FGNet")
library(topGO)
biocLite("topGO")
library(topGO)
library(ALL)
biocLite("ALL")
data(ALL)
data("ALL")
library(ALL)
data("ALL")
ls()
head(ALL)
data(geneList)
ls90
ls()
topGOdata
data(topGOdata)
library(InTAD)
inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
# filter gene expression
inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding")
# combine signals and genes in TADs
inTadSig <- combineInTAD(inTadSig, tadGR)
# perform correlation anlaysis
corData <- findCorrelation(inTadSig)
corData
head(corData)
ls()
inTadSig
inTadSig@geneCoords
ls()
data("mbSelEnhSignals") # loading enhancer signals
enhSel[1:3,1:3]
data("enhSelCoords")  # loading enhancer signal coordinates 
as.data.frame(enhSelGR[1:3])
data("mbSamplesRPKM")# loading gene expression counts
rpkmCountsSel[1:3,1:3]
data("txsSel") # loading gene coordiantes
as.data.frame(txsSel[1:3]
as.data.frame(txsSel[1:3])
ls()
rpkmCountsSel
gsub("\\..*", "", row.names(rpkmCountsSel))
unique(gsub("\\..*", "", row.names(rpkmCountsSel)))
ensembl=unique(gsub("\\..*", "", row.names(rpkmCountsSel)))
suppressMessages(library(clusterProfiler))
library(org.Hs.eg.db)
up=bitr(ensembl, fromType="ENSEMBL", toType="UNIPROT", OrgDb=org.Hs.eg.db)
head(up)
length(unique(up$UNIPROT))
ego2=enrichGO(gene=up$UNIPROT, keytype='UNIPROT', OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01)
ego2
ego2@result
head(ego2@result)
barplot(ego2, showCategory=8)
barplot(ego2, showCategory=20)
dotplot(ego2)
enrichMap(ego2)
class(ego2@result)
ego.res=(ego2@result)
head(ego.res)
subset(ego.res, ego.res$p.adjust<0.05)
dim(subset(ego.res, ego.res$p.adjust<0.05))
head(subset(ego.res, ego.res$p.adjust<0.05))
dim(ego.res)
dim(subset(ego.res, ego.res$p.adjust<0.05))
barplot
barplot()
?barplot
barplot(ego.res)
barplot(ego)
barplot(ego2)
barplot(ego2, showCategory=5)
barplot(ego2, showCategory=12)
history
history()
savehistory("test.R")
