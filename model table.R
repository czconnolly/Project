x<-read.dna('ALDH1A2_musclealignment.mfa', format='fasta')
write.table(modelTest(as.phyDat(x),G=F,I=F),'ALDH1A2modeltest.csv',sep='\t', quote=FALSE)