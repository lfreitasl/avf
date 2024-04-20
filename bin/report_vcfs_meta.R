#! micromamba run -n base Rscript

##Script for producing metadata to generated vcfs from dartR
library(dartR)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)
meta<-args[1]
pref<-args[2]
meta<-read.csv(meta)



filenames<-list.files("./",pattern = "\\.vcf$")
files<-list.files("./",pattern = "\\.vcf$",full.names = T)
vcfs<-lapply(files, read.vcfR)
gls<-lapply(vcfs, vcfR2genlight)
n_locs<-sapply(gls, function(x) x$n.loc)
n_inds<-sapply(gls, function(x) length(x$ind.names))

#Making dataframe to export
vcf_meta<-data.frame(filenames,n_inds,n_locs)

#Exporting information
write.csv(vcf_meta,"./vcfs_info.csv",row.names = F,quote = F)

myg<-gl.read.vcf(files)
meta<-meta[meta$samples%in%myg$ind.names,]
meta<-meta[match(myg$ind.names, meta$samples),]
write.csv(meta, paste(pref,"sorted_meta.csv", sep="_"),row.names = F,quote = F)