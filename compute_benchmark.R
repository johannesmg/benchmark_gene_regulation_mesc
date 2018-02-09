prepare.score.gs <- function(score.object,gs.table){
    prepared.score <- melt(score.object[["mat"]])
    prepared.score$size <- score.object[["size"]]
    prepared.score$colidx <- score.object[["colidx"]]
    prepared.score[,1] <- as.character(prepared.score[,1])
    prepared.score[,2] <- as.character(prepared.score[,2])
    print("concatenating TF-target strings")
    prepared.score$pair <- apply(prepared.score[,c(2,1)],1,str_c,collapse="_")
    gs.table$pair <- apply(gs.table[,c("TF","target")],1,str_c,collapse="_")
    prepared.score <- left_join(prepared.score[,c("pair","value","size","colidx")],gs.table,by="pair")
    prepared.score <- gather(prepared.score,gs,is.target,ChIP_seq:Overexpression)
    prepared.score <- prepared.score[,c("pair","value","size","gs","colidx","is.target")]
    prepared.score <- separate(prepared.score,pair,c("TF","target"),sep="_")
    prepared.score[is.na(prepared.score)] <- 0
    prepared.score
}


benchmark.size.dependent <- function(score.object,gs.table,prec.stop=0.05){
    prepared.score <- prepare.score.gs(score.object,gs.table)
    print("computing benchmark")
    prepared.score %>%
        group_by(gs,TF,size) %>%
        mutate(ss=sum(is.target)) %>%
        ungroup %>%
        group_by(gs,size,colidx) %>%
                                        #we have to restrict to TFs present in a given gold standard
        mutate(no.targets=sum(is.target[ss>0]),potential.targets=length(is.target[ss>0])) %>%
        summarise(aupr=performance(prediction(predictions=value[ss>0],labels=is.target[ss>0]),measure="aucpr",prec.stop=prec.stop)@y.values[[1]],no.targets=no.targets[1],potential.targets=potential.targets[1],raupr=no.targets[1]/potential.targets[1]*prec.stop,aupri=aupr/raupr)->
        aupri.network

    aupri.network
}



benchmark.size.dependent.tf <- function(score.object,gs.table,prec.stop=0.05){
    prepared.score <- prepare.score.gs(score.object,gs.table)
    print("computing benchmark")
    prepared.score %>%
        group_by(TF,gs,size) %>%
        mutate(ss=sum(is.target)) %>%
#        ungroup %>%
        filter(TF %in% TF[ss>=100]) %>%
        group_by(gs,TF,size) %>%
        summarise(raupr=sum(is.target)/length(is.target)*prec.stop,aupr=performance(prediction(predictions=value,labels=is.target),measure="aucpr",prec.stop=prec.stop)@y.values[[1]],aupri=aupr/raupr) ->
        by.tf.kd

    by.tf.kd
}


set.lower.tri.zero <- function(score.matrix,ensg.ids){
   foo <- score.matrix[ensg.ids,ensg.ids]
   foo[lower.tri(foo)] <- 0
   diag(foo) <- 0
   score.matrix[ensg.ids,ensg.ids] <- foo
   score.matrix
}


reformat.score.matrix <- function(eset,score.matrix,factors,set.lower.tri.zero=TRUE){
    rownames(score.matrix) <- rownames(eset)
    colnames(score.matrix) <- rownames(eset)
    ret.score <- score.matrix[,factors]
    if (set.lower.tri.zero){
        ret.score <- set.lower.tri.zero(ret.score,factors)
    }
    colnames(ret.score) <- names(factors)[match(factors,colnames(ret.score))]
    ret.score
}

special.plotting.fn <- function(gs.in){
    filter(benchmark.tf,gs==gs.in,L1 %in% c("ARACNe_50","CLR_MI","Pcor_pls")) %>%
    dplyr::select(one_of(c("TF","aupri","L1"))) %>%
    mutate(aupri=log2(aupri)) %>%
    spread(L1,aupri) %>%
        ggpairs(columns=2:4,lower = list(continuous = wrap("points", size = 0.5)))->p1
    p1
}

plotting.modifier <- function(po){
    for (i in 2:3){
    for (j in (1:(i-1))){
        po[i,j]+
            geom_hline(yintercept=0,colour="darkgray",linetype="dotted")+
            geom_vline(xintercept=0,colour="darkgray",linetype="dotted")+
            geom_abline(intercept=0,slope=1,colour="gray40")+
            theme(plot.margin=unit(rep(5,4),"cm"))->
            po[i,j]
    }
}
po
}


factors.ordered.hclust <- function(mat,reorder.col=TRUE,reorder.row=TRUE){
    mat.in <- as.matrix(as.data.frame(mat)[,-1])
    rownames(mat.in) <- as.data.frame(mat)[,1]
    dd.col <- as.dendrogram(hclust(dist(t(mat.in))))
    col.ord <- order.dendrogram(dd.col)
    dd.row <- as.dendrogram(hclust(dist((mat.in))))
    row.ord <- order.dendrogram(dd.row)
    if (reorder.col){
        mat.in <- mat.in[,col.ord]
    }        
    if (reorder.row){
        mat.in <- mat.in[row.ord,]
    }        
    xx_names <- attr(mat.in, "dimnames")
    df <- as.data.frame(mat.in)
    colnames(df) <- xx_names[[2]]
    df$row <- xx_names[[1]]
    df$row <- with(df, factor(row, levels=row, ordered=TRUE))
    df <- melt(df)
#    df$col <- xx_names[[2]]
    df$variable <- with(df, factor(variable, levels=colnames(mat.in), ordered=TRUE))
    df
}

library(tidyverse)
library(parallel)
library(reshape2)
library(stringr)
library(ROCR)
library(ggplot2)
library(cowplot)
library(GGally)

source("./raucpr/precision_recall.r")
source("aucpr_rocr_interface.R")


#################################################################
## BENCHMARKING COMPUTATIONS
#################################################################

max.prec=0.05
gs_table <- read_csv(file="./gold_standard_table.csv")
available.transcription.factors <- deframe(read_csv("available_transcription_factors.csv"))

load("./expression_set_gpl_6246_full.RData")
eset <- eset[available.transcription.factors,]
eset <- eset[,]

interaction.score.predict <- list()
a <- load(file="spearman_score_6246.RData")
interaction.score.predict[["Spearman"]] <- reformat.score.matrix(eset=eset,score.matrix=spearman.score,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="pearson_score_6246.RData")
interaction.score.predict[["Pearson"]] <- reformat.score.matrix(eset=eset,score.matrix=pearson.score,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="pcor_lasso_k_3_adalasso_gpl_6246_sol.RData")
interaction.score.predict[["Pcor_lasso"]] <- reformat.score.matrix(eset=eset,score.matrix=lasso.adalasso.sol[["pcor.lasso"]],factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="pcor_pls_k_3_gpl_6246.RData")
interaction.score.predict[["Pcor_pls"]] <- reformat.score.matrix(eset=eset,score.matrix=pls.sol[["pcor"]],factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="aracne_gpl_6246_tau_0.15.RData")
interaction.score.predict[["ARACNe_15"]] <- reformat.score.matrix(eset=eset,score.matrix=save.object[["aracne"]],factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="aracne_gpl_6246_tau_0.5.RData")
interaction.score.predict[["ARACNe_50"]] <- reformat.score.matrix(eset=eset,score.matrix=save.object[["aracne"]],factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="mi_parmigene_gpl_6246.RData")
interaction.score.predict[["MI"]] <- reformat.score.matrix(eset=eset,score.matrix=all.mi,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="mrnet_parmigene_gpl_6246.RData")
interaction.score.predict[["MRNET"]] <- reformat.score.matrix(eset=eset,score.matrix=mrnet.mat,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="clr_pearson_gpl_6246.RData")
interaction.score.predict[["CLR_Pearson"]] <- reformat.score.matrix(eset=eset,score.matrix=clr.mat,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
a <- load(file="clr_parmigene_gpl_6246.RData")
interaction.score.predict[["CLR_MI"]] <- reformat.score.matrix(eset=eset,score.matrix=clr.mat,factors=available.transcription.factors)
eval(parse(text=paste("rm(",a,")",sep="")))
gc()


matrix.list <- mclapply(interaction.score.predict,function(x){rv <- list();rv[["mat"]] <- x;rv$colidx=1;rv$size=ncol(eset);rv},mc.cores=length(interaction.score.predict))
top.size.benchmark <- mclapply(matrix.list,benchmark.size.dependent,gs.table=gs_table,prec.stop=max.prec,mc.cores=length(matrix.list))
top.size.benchmark <- do.call("bind_rows",c(top.size.benchmark,.id="score"))
save(top.size.benchmark,file="top_size_benchmark.RData")

matrix.list <- mclapply(interaction.score.predict,function(x){rv <- list(); rv[["mat"]] <- x;rv$colidx=1;rv$size=1;rv},mc.cores=length(interaction.score.predict))
per.tf.benchmark <- mclapply(matrix.list,benchmark.size.dependent.tf,gs.table=gs_table,prec.stop=max.prec,mc.cores=length(matrix.list))
save(per.tf.benchmark,file="per_tf_benchmark.RData")



ff <- list.files(path="./",pattern="*.RData",full.names=TRUE)
ff <- ff[grep("clr_6246_fixed",ff)]
matrix.list <- mclapply(ff,function(x){load(file=x);rv <- list(); rv[["mat"]] <- reformat.score.matrix(eset=eset,score.matrix=clr.out[["clr"]],factors=available.transcription.factors);rv$colidx=sum(clr.out[["columns"]]);rv$size=length(clr.out[["columns"]]);rv},mc.cores=length(ff))
clr.network.benchmark <- mclapply(matrix.list,benchmark.size.dependent,gs.table=gs_table,prec.stop=max.prec,mc.cores=length(matrix.list))
save(clr.network.benchmark,file="clr_size_benchmark.RData")

ff <- list.files(path="./",pattern="*.RData",full.names=TRUE)
ff <- ff[grep("clr_pearson_6246_fixed",ff)]
matrix.list <- mclapply(ff,function(x){load(file=x);rv <- list(); rv[["mat"]] <- reformat.score.matrix(eset=eset,score.matrix=clr.out[["clr"]],factors=available.transcription.factors);rv$colidx=sum(clr.out[["columns"]]);rv$size=length(clr.out[["columns"]]);rv},mc.cores=length(ff))
clr.pearson.network.benchmark <- mclapply(matrix.list,benchmark.size.dependent,gs.table=gs_table,prec.stop=max.prec,mc.cores=length(matrix.list))
save(clr.pearson.network.benchmark,file="clr_pearson_size_benchmark.RData")

ff <- list.files(path="./",pattern="*.RData",full.names=TRUE)
ff <- ff[grep("pls_net_6246_fixed_sampled",ff)]
matrix.list <- mclapply(ff,function(x){load(file=x);rv <- list(); rv[["mat"]] <- reformat.score.matrix(eset=eset,score.matrix=pls.out[["pls"]][["pcor"]],factors=available.transcription.factors);rv$colidx=sum(pls.out[["columns"]]);rv$size=length(pls.out[["columns"]]);rv},mc.cores=26)
pls.network.benchmark <- mclapply(matrix.list,benchmark.size.dependent,gs.table=gs_table,prec.stop=max.prec,mc.cores=length(matrix.list))
save(pls.network.benchmark,file="pls_net_size_benchmark.RData")



#################################################################
## BENCHMARKING FIGURES
#################################################################


load(file="top_size_benchmark.RData")

top.size.benchmark %>%
    ungroup %>%
    dplyr::select(score,gs,aupr,no.targets,potential.targets,raupr,aupri) %>%
    filter(!score %in% c("pcor_mean","random")) %>%
    dplyr::select(one_of(c("score","gs","aupri"))) %>%
    spread(gs,aupri) ->
    out.mat

factors.ordered.hclust(out.mat) %>%
    group_by(variable) %>%
    mutate(ff=round(value,1),rk=max(rank(value))-rank(value),value=(value-1)/max(value)) %>%
    ggplot(aes(x=variable,y=row,fill=value,label=ff))+
    geom_tile()+
    geom_text(size=2)+
    scale_fill_gradient2("Performance",breaks=c(-0.25,0,0.99),labels=c("min","random","max"),low="blue",high="red")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position = "top",legend.text = element_text(angle=90))+
    scale_x_discrete("")+
    scale_y_discrete("") ->
    p.bench.hm

group_by(top.size.benchmark, gs) %>%
    mutate(rk=max(rank(aupri))-rank(aupri)) %>%
    group_by(score) %>%
    summarise(mrk=mean(rk)) %>%
    ungroup %>%
    arrange(mrk)

save_plot(file="benchmark_top_size.pdf",plot=p.bench.hm,ncol=1,nrow=1,base_height=4,device=cairo_pdf,base_aspect_ratio=1.3)

#################################################################
#################################################################
#################################################################
#################################################################



load(file="pls_net_size_benchmark.RData")
load(file="clr_size_benchmark.RData")
load(file="clr_pearson_size_benchmark.RData")

network.benchmark.all <- bind_rows(do.call("rbind",clr.network.benchmark),do.call("rbind",pls.network.benchmark),do.call("rbind",clr.pearson.network.benchmark),.id="score")
network.benchmark.all$score <- c("clr_mi","pcor_pls","clr_pearson")[as.numeric(network.benchmark.all$score)]

network.benchmark.all %>%
    ungroup %>%
    mutate(gs=ifelse(gs=="Overexpression","Overexpr.",gs)) ->
    network.benchmark.all

group_by(network.benchmark.all,gs,size,score) %>%
    mutate(aupri=log2(aupri)) %>%
    summarise(ymax=median(aupri)+sd(aupri),ymin=median(aupri)-sd(aupri),aupri=median(aupri)) %>%
    ggplot(aes(x=size,y=aupri,group=gs,colour=gs,ymax=ymax,ymin=ymin))+
    geom_line(size=0.4)+
    geom_point(size=0.5)+
    geom_errorbar(size=0.3)+
    scale_y_continuous("Log2 Performance")+
    scale_x_continuous("# Arrays")+
    scale_colour_brewer("",palette="Set1")+
    facet_wrap(~score,ncol=4)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position="top",axis.line=element_line()) ->
    p.bench.all

save_plot(file="benchmark_sampling.pdf",plot=p.bench.all,ncol=1,nrow=1,base_height=4,device=cairo_pdf,base_aspect_ratio=1.3)


#################################################################
#################################################################
#################################################################
#################################################################


load(file="per_tf_benchmark.RData")

benchmark.tf <- melt(per.tf.benchmark,id.vars=1:6)
p.tf.ind <- lapply(unique(benchmark.tf$gs),special.plotting.fn)
names(p.tf.ind) <- unique(benchmark.tf$gs)
p.tf.ind <- lapply(p.tf.ind,plotting.modifier)
p.bench.ind <- plot_grid(plot_grid(p.tf.ind[[1]][2,1], p.tf.ind[[1]][3,1]+labs(title=names(p.tf.ind)[1])+ theme(plot.title = element_text(hjust = 0.5),axis.line=element_line(size=0.5)), p.tf.ind[[1]][3,2],ncol=3,align="h"),plot_grid(p.tf.ind[[2]][2,1], p.tf.ind[[2]][3,1]+labs(title=names(p.tf.ind)[2])+ theme(plot.title = element_text(hjust = 0.5),axis.line=element_line()), p.tf.ind[[2]][3,2],ncol=3,align="h"),plot_grid(p.tf.ind[[3]][2,1], p.tf.ind[[3]][3,1]+labs(title=names(p.tf.ind)[3])+ theme(plot.title = element_text(hjust = 0.5),axis.line=element_line()), p.tf.ind[[3]][3,2],ncol=3,align="h"),plot_grid(p.tf.ind[[4]][2,1], p.tf.ind[[4]][3,1]+labs(title=names(p.tf.ind)[4])+ theme(plot.title = element_text(hjust = 0.5),axis.line=element_line()), p.tf.ind[[4]][3,2],ncol=3,align="h"),align="vh",nrow=4,ncol=1,scale=0.9,vjust=-1)

save_plot(file="benchmark_individual_tfs.pdf",plot=p.bench.ind,ncol=1,nrow=1,base_height=4,device=cairo_pdf,base_aspect_ratio=1.3)

