require(DESeq2)
require(RCurl)
require(biomaRt)



dat<-read.delim('/homes/czconnolly/BS32010/Projectgit/RNAseqCounts.txt',
                sep='\t',skip=1,head=T)

#remoteFile <- getURL(paste0("http://www.compbio.dundee.ac.uk/",
#														 "user/pschofield/Teaching/Bioinformatics/",
#														 "Data/RNAseqCounts.txt"))

#geneCounts<-read.delim(textConnection(remoteFile),
#											 head=T,sep="\t",skip=1,row.names=1)

nonZeroCounts<-dat[rowSums(dat[,6:11])>0,6:11]

treatments <- as.factor(substr(colnames(nonZeroCounts),1,1))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
															as.data.frame(treatments),
															design=~treatments)

dds$treatments <- relevel(dds$treatments,"U")

dds <- DESeq(dds)

DEresults <- results(dds)

# get annotations
#mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
#								host="www.ensembl.org", 
#								path="/biomart/martservice")

mart<-useMart(biomart='ENSEMBL_MART_ENSEMBL',
              host='www.ensemble.org',
              path='/biomart/martservice')

#hg19 <- useDataset("hsapiens_gene_ensembl",mart=mart)

gg4<-useDataset('ggallus_gene_ensemble',mart=mart)

#annot <- getBM(attributes=c("ensembl_gene_id",
#														"external_gene_id",
#														"affy_hg_u133b"),
#							 filter="ensembl_gene_id",
#							 values=rownames(DEresults),
#							 mart=hg19)

annot<-getBM(attributes=c('ensemble_gene_id',
                          'external_gene_id',
                          'affy_chicken'),
             filter='ensemble_gene_id', 
             values=rownames(DEresults), 
             mart=gg4)

annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")
