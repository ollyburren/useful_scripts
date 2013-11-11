library(GenomicRanges)
library(rtracklayer)
        
## list of target diseases and their id's 
## note in next release these will be accessible by disease code
disease_lu<-data.frame(id=c(14,24,12,2,25,7,20,21,3,13,1,5),
	disease=c('AS','ATD','Celiac','Crohn','JIA','MS','PBC','Psoriasis','RA','SLE','T1D','UC')
	)

## import latest susceptibility regions from immunobase
## note replace build param for build 36
url_template<-'http://www.immunobase.org/webservice/RegionDownloads/=/model/tableBED/species=Human&disease_id=__REPLACE__&type=assoc&build=GRCh37'
        
imb.regions<-GRangesList(lapply(disease_lu$disease,function(x){
	did<-disease_lu[disease_lu$disease==x,]$id
	print(paste("Getting regions for",x,sep=" "))
	durl<-gsub("__REPLACE__",did,url_template)
	import.bed(durl)
}))
names(imb.regions)<-disease_lu$disease
imb.regions<-unlist(imb.regions)

## does one query snp overlap any of these ?

## example snp as a data.frame
## usually these will be obtained from file or pre-existing r object.
example.snp.df<-data.frame(chr='chr19',position=10475652,rsid='rs2304256')

## convert to a GR object   

example.gr<-with(example.snp.df,GRanges(seqname=Rle(chr),ranges=IRanges(start=position,end=position),rsid=rsid))

##because we have only one snp can do this ..
subsetByOverlaps(imb.regions,example.gr)
##GRanges with 10 ranges and 1 metadata column:
##            seqnames               ranges strand |        name
##               <Rle>            <IRanges>  <Rle> | <character>
##         AS    chr19 [10409793, 10631375]      * |     19p13.2
##      Crohn    chr19 [10399000, 10639000]      * |     19p13.2
##        JIA    chr19 [10390709, 10628548]      * |     19p13.2
##         MS    chr19 [10390709, 10634264]      * |     19p13.2
##        PBC    chr19 [10390709, 10628548]      * |     19p13.2
##  Psoriasis    chr19 [10390709, 10628548]      * |     19p13.2
##         RA    chr19 [10427721, 10492274]      * |     19p13.2
##        SLE    chr19 [10472933, 10569000]      * |     19p13.2
##        T1D    chr19 [10395447, 10628468]      * |     19p13.2
##         UC    chr19 [10220000, 10760000]      * |     19p13.2

## When we have more than one query SNP there could be multiple SNPs in a region.

example.snp.df<-data.frame(
	chr=c('chr19','chr1','chr1'),
	position=c(10475652,114377568,114303808),
	rsid=c('rs2304256','rs2476601','rs6679677'))    
example.gr<-with(example.snp.df,GRanges(seqname=Rle(chr),ranges=IRanges(start=position,end=position),rsid=rsid))

ol<-as.matrix(findOverlaps(imb.regions,example.gr))
## turn matrix into list (list name is overlapping index in imb.regions,
## list elements are overlapping indices in example.gr

lol<-split(ol[,2],ol[,1])

## create a report.
do.call('rbind',lapply(seq_along(lol),function(x){
	region<-imb.regions[as.numeric(names(lol)[x]),]
	region<-paste(names(region),region$name,sep=":")
	snps<-paste(example.gr[lol[[x]],]$rsid,sep=",",collapse=",")
	data.frame(region=region,snps=snps)
	}))

                 region                snps
## 1         AS:19p13.2           rs2304256
## 2         ATD:1p13.2 rs2476601,rs6679677
## 3       Crohn:1p13.2 rs2476601,rs6679677
## 4      Crohn:19p13.2           rs2304256
## 5         JIA:1p13.2 rs2476601,rs6679677
## 6        JIA:19p13.2           rs2304256                             
## 7         MS:19p13.2           rs2304256
## 8        PBC:19p13.2           rs2304256
## 9  Psoriasis:19p13.2           rs2304256
## 10         RA:1p13.2 rs2476601,rs6679677
## 11        RA:19p13.2           rs2304256
## 12        SLE:1p13.2 rs2476601,rs6679677                    
## 13       SLE:19p13.2           rs2304256
## 14        T1D:1p13.2 rs2476601,rs6679677
## 15       T1D:19p13.2           rs2304256
## 16        UC:19p13.2           rs2304256

              
