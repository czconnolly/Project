#install.packages('DESeq2')
#source('http://bioconductor.org/biocLite.R')
#biocLite('DESeq2')
require(DESeq2)
require(RCurl)
require(biomaRt)
require(stringr)


dat<-read.delim('/homes/czconnolly/BS32010/Projectgit/RNAseqCounts.txt',
                sep='\t',skip=1,head=T, row.names=1)

nonZeroCounts<-dat[rowSums(dat[,15:21])>0,15:21]

treatments <-factor(c('minus','plus','minus','plus','minus','plus','minus'))
#treatments <- as.factor(str_sub(colnames(nonZeroCounts),-9,-5))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
															as.data.frame(treatments),
															design=~treatments)

dds$treatments <- relevel(dds$treatments,"minus")

dds <- DESeq(dds)

DEresults <- results(dds)
#DEresults<- res[order(DEresults$padj),]
#head(res)

#write.csv(as.data.frame(DEresults), file= 'treatments_treated_results.csv')

#mart<-useMart(biomart='ENSEMBL_MART_ENSEMBL',
#              host='www.ensembl.org',
#              path='/biomart/martservice')

#gg4<-useDataset('ggallus_gene_ensembl',mart=mart)

#annot<-getBM(attributes=c('ensembl_gene_id',
#                          'external_gene_id',
#                          'affy_chicken'),
 #            filter='ensembl_gene_id', 
  #           values=rownames(DEresults), 
 #            mart=gg4)

#annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")

#head(annotResults)
#write.table(annotResults, 'annotResults.txt', sep='\t', quote=FALSE)
#return(annotResults)