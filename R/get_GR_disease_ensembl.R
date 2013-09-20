##AUTHOR: Olly Burren
##DATE: 25/10/2012
##PURPOSE: Illustrate the use of bioconductor to get lists of genes and their positions
##NOTE: Code also included to illustrate use of GenomicRanges to get overlapping objects

library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#get the details for all genes
gene_details<-getBM(attributes= c('ensembl_gene_id','external_gene_id','chromosome_name','start_position','end_position','strand','gene_biotype'),mart=ensembl)
write.table(file="e67.gene.positions",row.names=F,quote=F,sep="\t")


##here is an example of using R to find overlapping features using GenomicRange objects
##you can ignore this if you have an alternative method - I find it incredibly powerful
##you need to install GenomicRanges and rtracklayer (for the convenience of import bed files via URLs)

#convert gene_details to a GRanges object
library(GenomicRanges)

hs.GRange<-GRanges(seqnames=Rle(c(paste('chr',gene_details$chromosome_name,sep=""))),
    ranges=IRanges(gene_details$start_position, gene_details$end_position,names=gene_details$ensembl_gene_id),
    strand=gene_details$strand,
    external_id=gene_details$external_gene_id,
    biotype=gene_details$gene_biotype)
   

##lets say that I want to find all genes within current T1D susceptibility regions

library(rtracklayer)
regions<-import.bed('http://www.t1dbase.org/webservice/RegionDownloads/=/model/tableBED/species=Human&disease_id=1&type=assoc&build=GRCh37');
ol<-subsetByOverlaps(hs.GRange,regions)
