#Script to plot phylogenetic trees, create distance matrices, calculate regression using ape package in R

#read in alignment
x<-read.dna('SLC45A3_musclealignment.mfa', format='fasta')
#carry out model test on alignment
mtSLC45A3<-modelTest(as.phyDat(x),G=F,I=F)
#calculates DNA distances
d<-dist.dna(x)

#creates trees using the algorithms
tr.upgmaSLC45A3<-upgma(d)
tr.njSLC45A3<-nj(d)
tr.bionjSLC45A3<-bionj(d)
tr.fastmeSLC45A3<-fastme.bal(d,nni=T,spr=T,tbr=T)

#calculates regression
dt.upgmaSLC45A3<-cophenetic(tr.upgmaSLC45A3)
dmatupgma<-as.matrix(d)
nmsupgma<-rownames(dmatupgma)
dt.upgmaSLC45A3<-dt.upgmaSLC45A3[nmsupgma,nmsupgma]
dt.upgmaSLC45A3<-as.dist(dt.upgmaSLC45A3)

dt.njSLC45A3<-cophenetic(tr.njSLC45A3)
dmatnj<-as.matrix(d)
nmsnj<-rownames(dmatnj)
dt.njSLC45A3<-dt.njSLC45A3[nmsnj,nmsnj]
dt.njSLC45A3<-as.dist(dt.njSLC45A3)

dt.bionjSLC45A3<-cophenetic(tr.bionjSLC45A3)
dmatbionj<-as.matrix(d)
nmsbionj<-rownames(dmatbionj)
dt.bionjSLC45A3<-dt.bionjSLC45A3[nmsbionj,nmsbionj]
dt.bionjSLC45A3<-as.dist(dt.bionjSLC45A3)

dt.fastmeSLC45A3<-cophenetic(tr.fastmeSLC45A3)
dmatfastme<-as.matrix(d)
nmsfastme<-rownames(dmatfastme)
dt.fastmeSLC45A3<-dt.fastmeSLC45A3[nmsfastme,nmsfastme]
dt.fastmeSLC45A3<-as.dist(dt.fastmeSLC45A3)

#plots the trees
par(mfrow=c(2,2))
plot(tr.upgmaSLC45A3, main='UPGMA')
plot(tr.fastmeSLC45A3, main='FASTME')
plot(tr.njSLC45A3, main='NJ')
plot(tr.bionjSLC45A3, main='BIONJ')
