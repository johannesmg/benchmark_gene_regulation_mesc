library(tidyverse)
library(parmigene)
library(parcor)
library(parallel)

load("expression_set_gpl_6246_full.RData")
available.transcription.factors <- deframe(read_csv("available_transcription_factors.csv"))
eset <- eset[available.transcription.factors,]


sample.cols <- lapply(2^(4:10),function(x){replicate(5,sample(1:ncol(eset),x),simplify=FALSE)})
sample.cols[[length(sample.cols)+1]] <- replicate(1,sample(1:ncol(eset),ncol(eset)),simplify=FALSE)
sample.cols <- do.call("c",sample.cols)
save(sample.cols,file="eset_random_samples.RData")

pls.sample <- function(sample.col.in,eset.in){
    pls.out <- list()
    pls.out[["columns"]] <- sample.col.in
    pls.out[["pls"]] <- pls.net(X=t(eset.in[,sample.col.in]),k=3,verbose=FALSE)
    save(pls.out,file=paste("pls_net_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
}

pearson.sample <- function(sample.col.in,eset.in){
    out.object <- list()
    out.object[["columns"]] <- sample.col.in
    out.object[["mat"]] <- cor(t(eset.in[,sample.col.in]))
    save(out.object,file=paste("pearson_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
}

clr.pearson.sample <- function(sample.col.in){
    load(file=paste("pearson_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
    clr.out <- list()
    clr.out[["columns"]] <- sample.col.in
    clr.out[["clr"]] <- parmigene::clr(out.object[["mat"]])
    print("saving")
    save(clr.out,file=paste("clr_pearson_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
}
    
mi.sample <- function(sample.col.in,eset.in){
    mi.out <- list()
    clr.out <- list()
    mi.out[["columns"]] <- sample.col.in
    clr.out[["columns"]] <- sample.col.in
    print("computing mi")
    mi.out[["mi"]] <- knnmi.all(eset.in[,sample.col.in],k=9)
    print("computing clr")
    clr.out[["clr"]] <- parmigene::clr(mi.out[["mi"]])
    print("saving")
    save(mi.out,file=paste("mi_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
    save(clr.out,file=paste("clr_6246_fixed_sampled_",paste(sample.col.in[1:10],collapse="_"),".RData",sep=""))
}

#aaa <- load(file="eset_sample_reconstruction.RData")
load(file="eset_random_samples.RData")


foo <- mclapply(sample.cols,pearson.sample,eset.in=eset[,],mc.cores=4)
foo <- mclapply(sample.cols,pls.sample,eset.in=eset[,],mc.cores=4)
foo <- mclapply(sample.cols,mi.sample,eset.in=eset[,],mc.cores=1)
foo <- mclapply(sample.cols,clr.pearson.sample,mc.cores=1)
