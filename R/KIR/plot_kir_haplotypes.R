library(ggbio)
library(GenomicRanges)

## quick GFF parser

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     with(gff,GRanges(seqname=Rle(gff$seqname),ranges=IRanges(start=start,end=end),feature=feature,att=attributes))
     #return(gff)
}


kir.gff.file<-'path to gff file create by create_kir_gff.pl'

##PLOT SELECTED

select.snps<-c(
'rs35987710',
'seq-rs2240936',
'seq-rs270790',
'seq-rs34876880',
'seq-rs35093205',
'seq-rs62132666',
'seq-t1d-19-59936350-A-G',
'rs12975559',
'seq-rs11669355',
'rs2569678',
'rs3865507',
'seq-rs10500318',
'seq-rs17173106',
'seq-rs35656676',
'seq-rs3865510',
'seq-rs55761930',
'seq-rs592645',
'seq-rs597598',
'seq-rs598452',
'seq-rs604077',
'seq-rs604999',
'seq-rs648689',
'seq-rs649216',
'seq-t1d-19-60007809-C-G',
'seq-t1d-19-60014013-A-C',
'rs4806585',
'seq-rs4806568',
'seq-t1d-19-60054973-T-C',
'seq-t1d-19-60056605-A-T',
'seq-t1d-19-60056721-C-T')

## we colour this SNP red.
hilight.snp<-"seq-rs597598"


kir.gr<-gffRead(kir.gff.file)
kir.grl<-split(kir.gr,seqnames(kir.gr))

## PLOT ALL

track.list<-lapply(kir.grl,function(x){
	probe.gr<-x[grep("probe",x$feature),]
	probe.gr$att<-gsub("id=","",probe.gr$att)
	probe.gr<-probe.gr[which(probe.gr$att %in% select.snps),]
	probe.gr$color<-'black'
	hilight.index<-grep(hilight.snp,probe.gr$att)
	probe.gr[grep(hilight.snp,probe.gr$att),]$color<-'red'
	gene.gr<-x[x$feature=="gene",]
	
	## annoyingly we need to use 'identity' here as 'stepping' throws us of.
	## drop to defining rectangles manually 
	gene.gr$ymax<-length(probe.gr) ## depends on number of probes to plot
	gene.gr$ymin<-0
	gene.gr$xmin<-start(gene.gr)
	gene.gr$xmax<-end(gene.gr)
	
	## get allele name if available ..
	gene.gr$att<-gsub('.*allele [\\"]*([^ \\"]+)[\\"]*[\\"]*.*',"\\1",gene.gr$att)
	## if not drop to allele.
	gene.gr$att<-gsub('.*gene [\\"]*([^ \\"]+)[\\"]*[\\"]*.*',"\\1",gene.gr$att)
	## add allele / gene names to plot.
	gene.names<- annotate("text", label = gene.gr$att ,x=start(gene.gr) + width(gene.gr)/2,y=length(probe.gr)/2,size = 5, colour = "white",angle=-90,alpha=0.8)
	
	## plot gene locations
	gene.part<- geom_rect(gene.gr,stat="identity",fill="grey",color="white",aes(ymax=ymax,ymin=ymin,xmin=xmin,xmax=xmax),alpha=0.8)
	
	## plot probe locations
	## use the one below for hilighting a snp
	# probe.part<- geom_rect(probe.gr,fill=probe.gr$color,color=probe.gr$color,aes(group=att,ylab=att),group.selfish=TRUE)
	## use this if we want to colour probeA and B differently
	probe.part<- geom_rect(probe.gr,aes(group=att,ylab=att,fill=feature,color=feature),group.selfish=TRUE)
	if(length(probe.gr) >0){
		ggplot() + gene.part + probe.part + gene.names + theme_bw() + scale_x_sequnit("kb")
	}else{
		## allow for no mapping probes
		gene.part  + gene.names +theme_bw() +   scale_x_sequnit("kb")
	}
})


tracks(track.list)

## ZOOM IN AND PLOT A SPECIFIC GENE

zoom.gene<-"KIR2DL4" ## an example.

zoom.track.list<-lapply(kir.grl,function(x){
	gene.gr<-x[grep(zoom.gene,x$att),]
	#start(gene.gr)<-start(gene.gr)-1000
	#end(gene.gr)<-end(gene.gr)+1000
	
	probe.gr<-x[grep("probe",x$feature),]
	probe.gr$att<-gsub("id=","",probe.gr$att)
	probe.gr<-probe.gr[which(probe.gr$att %in% select.snps),]
	probe.gr$color<-'black'
	hilight.index<-grep(hilight.snp,probe.gr$att)
	probe.gr[grep(hilight.snp,probe.gr$att),]$color<-'red'
	probe.gr<-subsetByOverlaps(probe.gr,gene.gr)
	probe.gr<-shift(probe.gr,start(gene.gr)*-1)
	gene.gr<-shift(gene.gr,start(gene.gr)*-1)
	
	## annoyingly we need to use 'identity' here as 'stepping' throws us of.
	## drop to defining rectangles manually 
	#gene.gr<-x[x$feature=="gene",]
	#gene.gr<-x[x$feature=="gene",]
	gene.gr$ymax<-length(probe.gr)
	gene.gr$ymin<-0
	gene.gr$xmin<-start(gene.gr)
	gene.gr$xmax<-end(gene.gr)
	
	## get allele name if available ..
	gene.gr$att<-gsub('.*allele [\\"]*([^ \\"]+)[\\"]*[\\"]*.*',"\\1",gene.gr$att)
	## if not drop to allele.
	gene.gr$att<-gsub('.*gene [\\"]*([^ \\"]+)[\\"]*[\\"]*.*',"\\1",gene.gr$att)
	gene.names<- annotate("text", label = gene.gr$att ,x=start(gene.gr) + width(gene.gr)/2,y=length(probe.gr)/2,size = 5, colour = "white",angle=-90,alpha=0.8)
	
	## plot gene locations
	gene.part<- geom_rect(gene.gr,stat="identity",fill="grey",color="white",aes(ymax=ymax,ymin=ymin,xmin=xmin,xmax=xmax),alpha=0.8)
	
	## plot probe locations
	#probe.part<- geom_rect(probe.gr,fill=probe.gr$color,color=probe.gr$color,aes(group=att,ylab=att),group.selfish=TRUE)
	probe.part<- geom_rect(probe.gr,aes(group=att,ylab=att,fill=feature,color=feature),group.selfish=TRUE)
	if(length(probe.gr) >0){
		ggplot() + gene.part + probe.part + gene.names + theme_bw() + scale_x_sequnit("kb")
	}else{
		## allow for no mapping probes
		gene.part  + gene.names +theme_bw() +   scale_x_sequnit("kb")
	}
})

tracks(zoom.track.list)



