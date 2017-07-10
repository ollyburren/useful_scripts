library(data.table)


#download data using wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
#untar to a directory

RECOMB.DIR<-'FILL ME IN - LOCATION OF FILES THAT YOU DOWNLOADED'


computeInterval<-function(file,req.cm=0.1){
    chr<-sub(".*(chr[^.]+)\\.txt","\\1",basename(file))
    message(sprintf("Computing %.1f cM blocks for %s",req.cm,chr))
    dat<-fread(file)
    setnames(dat,c('chr','end','rate','cm'))
    dat$start<-head(c(1,dat$end+1),-1)
    res<-do.call('rbind',lapply(split(dat,floor(dat$cm/req.cm)),function(z){
        c(min(z$start),max(z$end))
    }))
    if(max(res[,2]) < max(dat$end))
        res<-rbind(res,c(max(res[,2])+1,max(dat$end)))
    return(data.table(chr=chr,start=res[,1],end=res[,2]))
    
}

rfiles<-list.files(path=RECOMB.DIR,pattern='.*chr[^.]+\\.txt',full.names=TRUE)
results<-lapply(rfiles,computeInterval)

## checking for region size distro.
lapply(results,function(r){
    i<-IRanges(start=r$start,end=r$end)
    summary(width(i))
})



