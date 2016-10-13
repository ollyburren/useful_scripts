## code to draw a graph of trait-loci by year to show the explosion of gwas

library(data.table)

#gw<-fread("~/Downloads/gwas_catalog_v1.0-associations_e85_r2016-09-27.tsv")
gw.live<-fread("https://www.ebi.ac.uk/gwas/api/search/downloads/full")

gw.f<-gw[,.(REGION,DATE,`P-VALUE`,`DISEASE/TRAIT`)]
setnames(gw.f,c('region','date','pval','trait'))

gw.f<-subset(gw.f,pval<5e-8)

setkeyv(gw.f,c('trait','region'))
gw.f<-unique(gw.f)

gw.f$year<-sub("^([^\\-]+)\\-.*","\\1",gw.f$date)
gw.f[,list(yc=nrow(.SD)),by=year]
pm<-gw.f[,list(yc=nrow(.SD)),by=year]
pm<-pm[order(pm$year),]

library(ggplot2)
library(cowplot)
ggplot(data = pm, aes(x = year, y = cumsum(yc),group=1)) + geom_line() +
    geom_point()  +
    scale_x_discrete(labels = pm$year) + xlab("Year") + ylab("Trait-Loci") + theme_bw(base_size=22) + theme(axis.text.x = element_text(angle=90, hjust = 1))

+ theme(axis.text.x = element_text(angle=90, hjust = 1)) + theme_bw()
