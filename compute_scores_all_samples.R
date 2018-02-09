library(tidyverse)
library(parmigene)
library(parcor)


available.transcription.factors <- deframe(read_csv("available_transcription_factors.csv"))

load("expression_set_gpl_6246_full.RData")
eset <- eset[available.transcription.factors,]

pearson.score <- cor(t(eset),method="pearson")
save(pearson.score,file="pearson_score_6246.RData")

spearman.score <- cor(t(eset),method="spearman")
save(spearman.score,file="spearman_score_6246.RData")



all.mi <- knnmi.all(eset,k=9)
save(all.mi,file="mi_parmigene_gpl_6246.RData")
mrnet.mat <- parmigene::mrnet(all.mi)
save(mrnet.mat,file="mrnet_parmigene_gpl_6246.RData")
clr.mat <- parmigene::clr(all.mi)
save(clr.mat,file="clr_parmigene_gpl_6246.RData")



tau.sequence <- c(0.15,0.5)
for (tau in tau.sequence){
    print(tau)
    arac.mat <- aracne.m(all.mi,tau=tau)
    save.object <- list(aracne=arac.mat,tau=tau)
    save(save.object,file=paste("aracne_gpl_6246_tau_",as.character(tau),".RData",sep=""))
}


aa <- load(file="pearson_score_6246.RData")
clr.mat <- parmigene::clr(pearson.score)
save(clr.mat,file="clr_pearson_gpl_6246.RData")


print("start pls")
print(timestamp())
pls.sol <- pls.net(X=t(eset),k=3,verbose=TRUE)
save(pls.sol,file="pcor_pls_k_3_gpl_6246.RData")
lasso.adalasso.sol <- adalasso.net(X=t(eset),k=3,both=TRUE,verbose=TRUE)
save(lasso.adalasso.sol,file="pcor_lasso_k_3_adalasso_gpl_6246_sol.RData")


