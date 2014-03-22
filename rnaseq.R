require(DESeq2)
require(RCurl)
require(biomaRt)
require(stringr)


dat<-read.delim('/homes/czconnolly/BS32010/Projectgit/RNAseqCounts.txt',
                sep='\t',skip=1,head=T)

nonZeroCounts<-dat[rowSums(dat[,6:11])>0,6:11]

treatments <- as.factor(substr(colnames(nonZeroCounts),55,65))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
															as.data.frame(treatments),
															design=~treatments)

dds$treatments <- relevel(dds$treatments,"DRB1_minus.")

dds <- DESeq(dds)

DEresults <- results(dds)

mart<-useMart(biomart='ENSEMBL_MART_ENSEMBL',
              host='www.ensemble.org',
              path='/biomart/martservice')

gg4<-useDataset('ggallus_gene_ensemble',mart=mart)

annot<-getBM(attributes=c('ensemble_gene_id',
                          'external_gene_id',
                          'affy_chicken'),
             filter='ensemble_gene_id', 
             values=rownames(DEresults), 
             mart=gg4)

annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")
