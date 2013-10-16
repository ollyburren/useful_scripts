library(GenomicRanges)
library(rtracklayer)

## Author: Oliver Burren 
## Function: Demo script of liftOver in R !
## Demonstration of using liftOver from within R to remap a snp.
## note that GRanges can contain many features (in this case snps).


chr<-1;
start<-145079845
snp.name<-'rs4649771'

chain.file.path<-'path to your chain file'

## convert to Genomic Ranges

example.37.gr<-GRanges(
	seqname=Rle(paste("chr",chr,sep="")),
	ranges=IRanges(start=start,end=start),
	snp.name=snp.name)
	

example.37.gr
## GRanges with 1 range and 1 metadata column:
##       seqnames                 ranges strand |    snp.name
##          <Rle>              <IRanges>  <Rle> | <character>
##   [1]     chr1 [145079845, 145079845]      * |   rs4649771
##   ---
##   seqlengths:
##    chr1
##      NA


c<-import.chain(chain.file.path) ## e.g. hg19ToHg18.over.chain
example.36.gr<-unlist(liftOver(example.37.gr,c))  
names(mcols(example.36.gr))<-'snp.name'
example.36.gr
## GRanges with 1 range and 1 metadata column:
##     seqnames                 ranges strand |    snp.name
##        <Rle>              <IRanges>  <Rle> | <character>
##   1     chr1 [143791202, 143791202]      * |   rs4649771
##   ---
##   seqlengths:
##    chr1
##      NA


