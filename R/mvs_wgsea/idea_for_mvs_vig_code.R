library(snpStats)
library(wgsea)

## Here I attempt to repeat the example in the current wgsea vignette but using mvs


## The mvs method is sensitive to population, if we include multiple populations with
## heterogeneous LD structure I think the test becomes anticonservative. 

data(for.exercise,package="snpStats")
snpsum<-col.summary(snps.10)
snp.indicator<-abs(snpsum$z.HWE) > 1.96

## Notice that I subset on stratum here otherwise the resultant test is anti-conservative !!

case<-snps.10[subject.support$cc==1 & subject.support$stratum=='CEU',!is.na(snp.indicator)]
control<-snps.10[subject.support$cc==0& subject.support$stratum=='CEU' ,!is.na(snp.indicator)]
snp.support.filt<-snp.support[!is.na(snp.indicator),]
snp.indicator<-snp.indicator[!is.na(snp.indicator)]


maf<-col.summary(control)[,"MAF"]
p<-pairtest(case,control)
## load recomb hotspots from hapmap
hs.10.url<-'http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/genetic_map_chr10_b36.txt'

library(GenomicRanges)


## a function for splitting a chromosome into regions based on 0.1cM genetic distance 
## returns a GRange object ..
gr.recomb<-function(url,chr,recomb.thresh=0.1){
	r.df<-read.table(file=url,header=TRUE,sep=" ")
	names(r.df)<-c('position','recombination_rate','genetic_map_position')
	r.df<-r.df[order(r.df$position),]
	gmap<-r.df$genetic_map_position
	r.df$interval<-trunc(gmap/recomb.thresh)+1
	r.df.split<-split(r.df,r.df$interval)
	#here we adjust intervals slightly to use max recomb rate as within an interval this makes sense
	interval.position<-sapply(r.df.split,function(x) x[which.max(x$recombination_rate),]$position)
	idx<-which(r.df$position %in% interval.position)
	r.df<-data.frame(start=r.df[head(idx,n=length(idx)-1),]$position+1,end=r.df[idx[-1],]$position)
	## this to capture the beginning and end of chromosome (999999999 is arbitrarily long)
	r.df<-rbind(r.df,data.frame(start=c(1,max(r.df$end)+1),end=c(min(r.df$start)-1,999999999)))	
	r.df<-r.df[order(r.df$start),]
	with(r.df,
		GRanges(seqnames=Rle(chr),
						IRanges(start=start,end=end,name=1:nrow(r.df)))
						)
}


## regions for chr10 using hapmap data
chr10.recomb.regions<-gr.recomb(hs.10.url,'chr10')

## create lists of SNPs for each recombination region
snp.support.gr<-with(snp.support.filt,GRanges(seqnames=Rle(paste('chr',chromosome,sep="")),
	ranges=IRanges(start=position,end=position),name=rownames(snp.support.filt)))

## which snp goes with which region
ol<-as.matrix(findOverlaps(chr10.recomb.regions,snp.support.gr))

## compute sigma for each region based on control genotypes

sigma.10<-lapply(split(ol[,2],ol[,1]),function(x){
	## allow for a genetic region not contain snps 
	if(length(x)==0)
		return(NA)
	if(length(x)==1){
		snp.name<-snp.support.gr[x,]$name
		return(Matrix(1,dimnames = list(snp.name,snp.name)))
	}
	mvs.sigma.ld(control[,x])
})

## next we compute 1000 permutations based on multivariate sampling

p.mvs.perm<-compute.mvs.perms_slow(sigma.10,1000,snp.support.gr)

## next do wilcoxon

W<-wilcoxon(p,snps.in=which(snp.indicator==1))
Wstar<-wilcoxon(p.mvs.perm,snps.in=which(snp.indicator==1))
Z.value(W=W,Wstar=Wstar,n.in=sum(snp.indicator==1),n.out=sum(snp.indicator==0))

