library(readr)
library(readxl)
library(data.table)
require(glmnet)
require(plyr)
require(dplyr)
require(ggplot2)
require(samr)
setwd('~/Dropbox/CSF metabolites/')
serum1 = read.csv('SERUM_Log2_S2Nfiltered.tsv', sep  = "\t")
csf1 = read_delim('CSF_Log2_S2Nfiltered.tsv', delim = "\t")
revised_groups = fread('~/Dropbox/CSF metabolites/narcosharadrevised_EM_Dec12.csv', header = T)
revised_groups=revised_groups[, -1, with =F]
serumIDS = read_csv('~/Documents/CSF metabolites/Serum Id_dec2_2016.csv')
neHCRTCSF=read_csv('~/Dropbox/CSF metabolites/CSF data from Shahrad with new HCRT and percentage result.csv')
neHCRTCSF$`A:%B(F)`=100-as.numeric(gsub("^\\.", NA, neHCRTCSF$`A:%B(F)`))
neHCRTCSF$`HCRT-11-4-16` = as.numeric(gsub("^\\.", NA, neHCRTCSF$`HCRT-11-4-16`))
HCRT_CSF_percent= select(neHCRTCSF, ExtIdentifier, SampleDescription, `A:%B(F)`)
colnames(HCRT_CSF_percent)[3] = c("HCRT_Percent")
######################
ids_ling = read_excel('~/Documents/CSF metabolites/Shahrad sera-CSF list with infor1-groups to EM5.xlsx', sheet = 1, col_names = T)
idstomerge = select(ids_ling, `ID to Sharhad (ST)`, groups)
colnames(idstomerge) = c("EMID", "Group")
#######################
serum_proteins =serum1[30:ncol(serum1)] 
serum_proteins = cbind.data.frame(serum1$Target, serum_proteins)
csf_proteins = csf1[30:ncol(csf1)]
csf_proteins = cbind.data.frame(csf1$Target, csf_proteins)
serum_proteins$MedianSignalToNoisePass2InAnyGroup = NULL
csf_proteins$MedianSignalToNoisePass2InAnyGroup = NULL
#transpose data frame
require(data.table)
csf_t = setNames(data.frame(t(csf_proteins[,-1])), csf_proteins[,1])
serum_t = setNames(data.frame(t(serum_proteins[,-1])), serum_proteins[,1])
##select only proteins that are present in CSF
#merge the data frames 
csf_t$EXTID = row.names(csf_t)
serum_t$EXTID = row.names(serum_t)

#######add in EMIDS into the data frame to match 
serum_IDmatched=merge.data.frame(serum_t, serumIDS, by.x = "EXTID", by.y = "ExtIdentifier")
colnames(serum_IDmatched)=paste0(colnames(serum_IDmatched), "_SERUM")
##########get IDS from the CSF file
require(plyr)
require(dplyr)
CSFIDS=select(CSF, SampleDescription, ExtIdentifier)
csf_IDmatched = merge.data.frame(csf_t, CSFIDS, by.x = "EXTID", by.y = "ExtIdentifier")
colnames(csf_IDmatched)=paste0(colnames(csf_IDmatched), "_CSF")
#######################################################################################
merg1 = merge.data.frame(csf_IDmatched, serum_IDmatched, by.x = "SampleDescription_CSF", by.y = "SampleId_SERUM")
merg1$EXTID_CSF = NULL
merg1$EXTID_SERUM = NULL
colnames(merg1)=gsub(" ", "_", colnames(merg1)) ##remove spaces from colnames
###perform cross-correlations between csf and serum app####################################### 
##CSF_APP/Serum_APP with all other proteins in the panel
which(colnames(merg1) == "Apo_E4_CSF") ####get column index
metnames = merg1[-c(1, 167)]
amyloid = merg1$Apo_E4_CSF
df = data.frame(label = paste(names(metnames)), Estimate=numeric(length(metnames)), P.value=numeric(length(metnames)))
for (i in 1:length(metnames)){
  test = cor.test(amyloid, metnames[[i]], method = "p")
  df$Estimate[i]=test$estimate
  df$P.value[i]=test$p.value
}
df=df[order(df$P.value),]
df=df[order(df$P.value),]
adjustedp1=p.adjust(df$P.value, method = "hochberg")
df=cbind(df, adjustedp1)
df
sigam_csf=filter(df, df$adjustedp1 <= 0.05)
write_csv(sigam, path = '~/Dropbox/CSF metabolites/Significant_correlation_amyloidbeta_Dec7.csv')
###############CSF_APP with other serum proteins##################################################
which(colnames(merg1) == "amyloid_precursor_protein_SERUM") ####get column index
require(stringr)
metnames = merg1[str_detect(colnames(merg1), "_SERUM")]
metnames = metnames[-1]
amyloid = merg1$Apo_E4_SERUM
df = data.frame(label = paste(names(metnames)), Estimate=numeric(length(metnames)), P.value=numeric(length(metnames)))
for (i in 1:length(metnames)){
  test = cor.test(amyloid, metnames[[i]], method = "p")
  df$Estimate[i]=test$estimate
  df$P.value[i]=test$p.value
}
df=df[order(df$P.value),]
adjustedp1=p.adjust(df$P.value, method = "hochberg")
df=cbind(df, adjustedp1)
df
require(plyr);require(dplyr)
sigam_csf=filter(df, df$adjustedp1 <= 0.05)
write_csv(sigam, path = '~/Dropbox/CSF metabolites/Significant_correlation_amyloidbeta_Dec7.csv')

#########################K-means clustering of all the 1725 proteins################################

total_matrix = data.matrix(merg1[-1])
total_matrix1_cor = cor(total_matrix)
total1_dist = dist(total_matrix1_cor, method = "e")
total1_cluster = hclust(total1_dist, method = "complete")
plot(total1_cluster)

clusters1=cutree(total1_cluster, h=11)
clusters2 = as.data.frame(clusters1)
clusters2$TargetID = row.names(clusters2)
clusters2
write.csv(clusters2, file = '~/Dropbox/total_serum_CSF_clusters_jan13.csv')
d <- cut(as.dendrogram(total1_cluster), h=10)

plot(d$lower[[1]])

#####################################################################################################
merg2= setNames(data.frame(t(merg1[,-1])), merg1[,1]) ##transpose the data frame
colnames(merg1) = gsub("\\.x", "_CSF", colnames(merg1))
colnames(merg1) = gsub("\\.y", "_SERUM", colnames(merg1))
colnames(merg1)
merg2= setNames(data.frame(t(merg1[,-1])), merg1[,1])
merg2$TargetIDS = row.names(merg2)
merg2$source=ifelse(grepl("\\.x", merg2$TargetIDS), "CSF", "SERUM")
merg2$TargetIDS = gsub("\\.x|\\.y", "", merg2$TargetIDS)
setDT(merg2)
merg2[, 1:10, with = FALSE]
#######################################################################################
merg2 = as.data.frame(merg2)
merg2_list=split.data.frame(merg2, merg2$source)
merg2_list$CSF$TargetIDS[!merg2_list$CSF$TargetIDS %in% merg2_list$SERUM$TargetIDS]
merg2_list$SERUM$TargetIDS[!merg2_list$SERUM$TargetIDS %in% merg2_list$CSF$TargetIDS]

#############################
library(limma)
require(reshape2)
require(data.table)
setDT(merg2)
mod1=dcast(merg2, TargetIDS~source,fun.aggregate = length, value.var = values)
mod1$TargetIDS = NULL
vennCounts(mod1) %>% vennDiagram()
########################################################################################
x=melt(merg2, id.vars = c("TargetIDS", "source"))
y=x[, list(value=mean(value)), by = list(TargetIDS, source)]
z=dcast(y, TargetIDS~source, value.var = "value")
z$CSF_ratio = z$CSF/z$SERUM
z$serum_ratio = z$SERUM/z$CSF
##melt the DF


long_merg2$Ratio[grep("Cystatin C", long_merg2$TargetIDS)] %>% hist(main = "Cystatin C - CSF/Serum ratio")
long_merg2$Ratio[grep("Cystatin M", long_merg2$TargetIDS)] %>% hist(main = "Cystatin C - CSF/Serum ratio")

###########CSF clusters###############################################################

csf_matrix = data.matrix(csf_IDmatched[2:616])
csf_matrix1_cor = cor(csf_matrix)
csf1_dist = dist(csf_matrix1_cor, method = "e")
csf_matrix1_cor
csf1_dist = dist(csf_matrix1_cor)
csf1_cluster = hclust(csf1_dist, method = "complete")
plot(csf1_cluster)
cutree(csf1_cluster, h=max(csf1_cluster$height/4))
clusters1=cutree(csf1_cluster, h=11)
clusters2 = as.data.frame(clusters1)
clusters2$TargetID = row.names(clusters2)
clusters2
########################################################################################
serum_matrix = data.matrix(serum_IDmatched[2:1110])
serum_matrix1_cor = cor(serum_matrix)
serum1_dist = dist(serum_matrix1_cor, method = "e")
serum1_cluster = hclust(serum1_dist, method = "complete")
plot(serum1_cluster)
cutree(serum1_cluster, h=max(csf1_cluster$height/4))
serum_clusters1=cutree(serum1_cluster, h=11)
serum_clusters2 = as.data.frame(serum_clusters1)
serum_clusters2$TargetID = row.names(serum_clusters2)
serum_clusters2
############################################get ratios from long_merg2 and match them with csf clusters
###aggregate long merg2 by mean ratios
setDT(long_merg2)
long_merg3 = long_merg2
Longmerg_ratios=long_merg3[, list(Ratio = mean(Ratio, na.rm = T)), by =list(TargetIDS)]
Longmerg_ratios = as.data.frame(Longmerg_ratios)
####################match ratios to clusters#######################################
clus_ratios_csf = merge.data.frame(z, clusters2,  by.x = "TargetIDS", by.y = "TargetID", all.x = T)
clus_ratios_serum = merge.data.frame(z, serum_clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)
arrange(clus_ratios_csf, desc(ratio))
arrange(clus_ratios_serum, desc(serum_ratio))
total_clusters = merge.data.frame(clus_ratios_serum, clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)
arrange(total_clusters, desc(CSF_ratio))
write.csv(total_clusters, file ='~/Desktop/Total_Clusters_Dec5.csv')
#######################subgroup###############################################analysis
csf_groups = merge.data.frame(csf_IDmatched, revised_groups, by.x = "SampleDescription", by.y = "EMID") ##get diagnosis information by merge
csf_groups$totalAPO = (csf_groups$`Apo E2`+csf_groups$`Apo E3`+csf_groups$`Apo E4`)/3
serum_groups = merge.data.frame(serum_IDmatched, revised_groups, by.x = "SampleId", by.y = "EMID") 
serum_groups$totalAPO = (serum_groups$`Apo E2`+serum_groups$`Apo E3`+serum_groups$`Apo E4`)/3

############ NARC -CSF ############################
csf_groups[1:10]
csf_narc = filter(csf_groups, Group == "1")
csf_narc$SampleDescription = NULL
csf_narc$EXTID = NULL
csf_narc$Group = NULL

csf_narc_matrix = data.matrix(csf_narc)
csf_narc_matrix1_cor = cor(csf_narc_matrix)
csf_narc1_dist = dist(csf_narc_matrix1_cor, method = "e")
csf_narc1_cluster = hclust(csf_narc1_dist, method = "complete")
plot(csf_narc1_cluster)
cutree(csf_narc1_cluster, h=max(csf1_cluster$height/5))
csf_narc_clusters1=cutree(csf_narc1_cluster, h=11)
csf_narc_clusters2 = as.data.frame(csf_narc_clusters1)
csf_narc_clusters2$TargetID = row.names(csf_narc_clusters2)
csf_narc_clusters2
############ CTRL -CSF ############################
csf_ctrl = filter(csf_groups, Group != "1")
###merge the HCRT data into the ctrls file for cluster analysis
csf_ctrl = merge.data.frame(csf_ctrl, HCRT_CSF_percent, by.x = "SampleDescription", by.y = "SampleDescription", all.x = T)
csf_ctrl$SampleDescription = NULL
csf_ctrl$EXTID = NULL
csf_ctrl$Group = NULL
csf_ctrl$HCRT_Percent = NULL
csf_ctrl_matrix = data.matrix(csf_ctrl)
csf_ctrl_matrix1_cor = cor(csf_ctrl_matrix)
csf_ctrl1_dist = dist(csf_ctrl_matrix1_cor, method = "e")
csf_ctrl1_cluster = hclust(csf_ctrl1_dist, method = "complete")
plot(csf_ctrl1_cluster)
cutree(csf_ctrl1_cluster, h=max(csf1_cluster$height/5))
csf_ctrl_clusters1=cutree(csf_ctrl1_cluster, h=10)
csf_ctrl_clusters2 = as.data.frame(csf_ctrl_clusters1)
csf_ctrl_clusters2$TargetID = row.names(csf_ctrl_clusters2)
csf_ctrl_clusters2
############ NARC -SERUM ############################
serum_narc = filter(serum_groups, Group == "1")
serum_narc$SampleId = NULL
serum_narc$x = NULL
serum_narc$Group = NULL

serum_narc_matrix = data.matrix(serum_narc)
serum_narc_matrix1_cor = cor(serum_narc_matrix)
serum_narc1_dist = dist(serum_narc_matrix1_cor, method = "e")
serum_narc1_cluster = hclust(serum_narc1_dist, method = "complete")
plot(serum_narc1_cluster)
#cutree(serum_narc1_cluster, h=max(serum1_cluster$height/5))
serum_narc_clusters1=cutree(serum_narc1_cluster, h=11)
serum_narc_clusters2 = as.data.frame(serum_narc_clusters1)
serum_narc_clusters2$TargetID = row.names(serum_narc_clusters2)
serum_narc_clusters2
############ CTRL -SERUM ############################
serum_ctrl = filter(serum_groups, Group != "1")
serum_ctrl$SampleId = NULL
serum_ctrl$EXTID = NULL
serum_ctrl$Group = NULL

serum_ctrl_matrix = data.matrix(serum_ctrl)
serum_ctrl_matrix1_cor = cor(serum_ctrl_matrix)
serum_ctrl1_dist = dist(serum_ctrl_matrix1_cor, method = "e")
serum_ctrl1_cluster = hclust(serum_ctrl1_dist, method = "complete")
plot(serum_ctrl1_cluster)
cutree(serum_ctrl1_cluster, h=max(serum1_cluster$height/5))
serum_ctrl_clusters1=cutree(serum_ctrl1_cluster, h=11)
serum_ctrl_clusters2 = as.data.frame(serum_ctrl_clusters1)
serum_ctrl_clusters2$TargetID = row.names(serum_ctrl_clusters2)
serum_ctrl_clusters2

#######################merge all the clusters to have one master sheet##############################
total_clusters$ratio = NULL
total_clusters = merge.data.frame(total_clusters, csf_narc_clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)
total_clusters = merge.data.frame(total_clusters, csf_ctrl_clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)
total_clusters = merge.data.frame(total_clusters, serum_narc_clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)
total_clusters = merge.data.frame(total_clusters, serum_ctrl_clusters2, by.x = "TargetIDS", by.y = "TargetID", all.x = T)

#####################################Kmeans################################# CSF ######################
csf_totalcor = as.matrix(cor(csf_groups[3:617]))
csf_total_allgrp = kmeans(csf_totalcor, centers = 10)
csf_narc_k = kmeans(csf_narc_matrix1_cor, centers = 10)
csf_ctrl_k = kmeans(csf_ctrl_matrix1_cor, centers = 10)
csf_total_allgrp_1 = as.data.frame(csf_total_allgrp$cluster)
csf_narc_K_1 = as.data.frame(csf_narc_k$cluster)
csf_ctrl_K_1 = as.data.frame(csf_ctrl_k$cluster)
csf_total_allgrp_1$TargetIDs = rownames(csf_total_allgrp_1)
csf_narc_K_1$TargetIDs = rownames(csf_narc_K_1)
csf_ctrl_K_1$TargetIDs = rownames(csf_ctrl_K_1)
colnames(csf_narc_K_1)[1] = c("CSF_NARCOLEPSY")
colnames(csf_ctrl_K_1)[1] = c("CSF_CTRLS")
colnames(csf_total_allgrp_1)[1] = c("TOTAL_BOTH")
#####################################Kmeans################################ SERUM ##################


serum_totalcor = as.matrix(cor(serum_groups[3:617]))
serum_total_allgrp = kmeans(serum_totalcor, centers = 10)
serum_narc_k = kmeans(serum_narc_matrix1_cor, centers = 10)
serum_ctrl_k = kmeans(serum_ctrl_matrix1_cor, centers = 10)
serum_total_allgrp_1 = as.data.frame(serum_total_allgrp$cluster)
serum_narc_K_1 = as.data.frame(serum_narc_k$cluster)
serum_ctrl_K_1 = as.data.frame(serum_ctrl_k$cluster)
serum_total_allgrp_1$TargetIDs = rownames(serum_total_allgrp_1)
serum_narc_K_1$TargetIDs = rownames(serum_narc_K_1)
serum_ctrl_K_1$TargetIDs = rownames(serum_ctrl_K_1)
colnames(serum_narc_K_1)[1] = c("serum_NARCOLEPSY")
colnames(serum_ctrl_K_1)[1] = c("serum_CTRLS")
colnames(serum_total_allgrp_1)[1] = c("TOTAL_BOTH")

n_ser=serum_narc_K_1[grep(7, serum_narc_K_1$serum_NARCOLEPSY),]
c_ser=serum_ctrl_K_1[grep(7, serum_ctrl_K_1$serum_CTRLS),]


################################################ samr ###############################################
require(samr) ##csf
y = ifelse(csf_groups$Group == "1", 1, 2)
x = t(csf_groups[3:617])
d = list(x=x, y=y, geneid=rownames(x), logged2 =F)

samr.obj<-samr(d,  resp.type="Two class unpaired", assay.type = "array", center.arrays = T, nperms = 1000)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=1.5
samr.plot(samr.obj,delta)


### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
ggplot(csf_groups, aes(factor(y), csf_groups$`sE-Selectin`, color = factor(Group)))+geom_point(size = 4) + stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")

##get all significant proteins from both groups ####
csf_upctrls=siggenes.table$genes.up %>% as.data.frame()
csf_upcases=siggenes.table$genes.lo %>% as.data.frame()
csf_samr=rbind.data.frame(csf_upcases, csf_upctrls)
csf_samr$`Score(d)`=as.numeric(as.character(csf_samr$`Score(d)`))
ggplot(csf_samr, aes(reorder(csf_samr$`Gene Name`, csf_samr$`Score(d)`), csf_samr$`Score(d)`)) + geom_bar(stat = "identity") +coord_flip()+theme(axis.text.y = element_text(size = 12))

####make a heatmap 
csf_samr_sub=csf_groups1[colnames(csf_groups1) %in% csf_samr$`Gene Name`]
csf_samr_mat = as.matrix(t(csf_samr_sub))
csf_samr_cor = dist(cor(csf_samr_mat), method = "e")
csf_samr_clust = hclust(csf_samr_cor, method = "ward.D2")
heatmap.2(csf_samr_mat, trace="none",col=greenred(10), Colv= as.dendrogram(csf_samr_clust), Rowv = NA, dendrogram = "col",density.info = "none", key = F)

################################################ SAMR SERUM ########################################

y = ifelse(serum_groups$Group == "1", 1, 2)
x = t(serum_groups[3:1111])
d = list(x=x, y=y, geneid=rownames(x), logged2 =F)

samr.obj<-samr(d,  resp.type="Two class unpaired", assay.type = "array", center.arrays = T, nperms = 1000)
delta=1.52
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)

### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
serum_upctrls=siggenes.table$genes.up %>% as.data.frame()
serum_upcases=siggenes.table$genes.lo %>% as.data.frame()
serum_samr=rbind.data.frame(serum_upcases, serum_upctrls)
serum_samr$`Score(d)`=as.numeric(as.character(serum_samr$`Score(d)`))
serum_upcases$`Score(d)`=as.numeric(as.character(serum_upcases$`Score(d)`))
serum_upctrls$`Score(d)`=as.numeric(as.character(serum_upctrls$`Score(d)`))
ggplot(serum_samr, aes(reorder(serum_samr$`Gene Name`, serum_samr$`Score(d)`), serum_samr$`Score(d)`)) + geom_bar(stat = "identity") +coord_flip()+theme(axis.text.y = element_text(size = 9))
ggplot(serum_upcases, aes(reorder(serum_upcases$`Gene Name`, serum_upcases$`Score(d)`), serum_upcases$`Score(d)`)) + geom_bar(stat = "identity") +coord_flip()+theme(axis.text.y = element_text(size = 12))
ggplot(serum_upctrls, aes(reorder(serum_upctrls$`Gene Name`, serum_upctrls$`Score(d)`), serum_upctrls$`Score(d)`)) + geom_bar(stat = "identity") +coord_flip()+theme(axis.text.y = element_text(size = 12))

####make a heatmap 
serum_samr_sub=serum_groups[colnames(serum_groups) %in% serum_samr$`Gene Name`]
serum_samr_mat = as.matrix(t(serum_samr_sub))
serum_samr_cor = dist(cor(serum_samr_mat), method = "e")
serum_samr_clust = hclust(serum_samr_cor, method = "ward.D2")
plot(serum_samr_clust)
heatmap.2(serum_samr_mat, trace="none",col=greenred(10), Rowv = NA, dendrogram = "col", Colv = as.dendrogram(serum_samr_clust), density.info = "none", key = F)





ggplot(serum_groups, aes(factor(y), serum_groups$Aggrecan, color = factor(Group)))+geom_point(size = 4) + stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
csf_total_K = merge.data.frame(csf_narc_K_1, csf_ctrl_K_1, by.x = "TargetIDs", by.y = "TargetIDs")
csf_total_K = merge.data.frame(csf_total_K, csf_total_allgrp_1,  by.x = "TargetIDs", by.y = "TargetIDs")
write.csv(csf_total_K, file ='~/Desktop/CSF_TOTAL_CLUSTERS_DEC6.csv')
###########test
require(vegan)
fit <- cascadeKM(csf_ctrl_matrix1_cor, 1, 10, iter = 50)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
################################################# HEATMAPS for cluster analysis ####################################################
csf_groups1 = filter(csf_groups, csf_groups$newgroups == 2 | csf_groups$newgroups == 3)
csf_groups1$SampleDescription = NULL
csf_groups1$EXTID = NULL
csf_groups1$Group = NULL
csf_groups1$SampleDescription = NULL

which(colnames(csf_groups1) == "amyloid precursor protein") ####get column index
metnames = csf_groups1[-c(616:621, 219)]
amyloid = csf_groups1$`amyloid precursor protein`
df = data.frame(label = paste(names(metnames)), Estimate=numeric(length(metnames)), P.value=numeric(length(metnames)))
for (i in 1:length(metnames)){
  test = cor.test(amyloid, metnames[[i]], method = "p")
  df$Estimate[i]=test$estimate
  df$P.value[i]=test$p.value
}
df=df[order(df$P.value),]
adjustedp1=p.adjust(df$P.value, method = "b")
df=cbind(df, adjustedp1)
df
sigam_csf=df#filter(df, df$adjustedp1 <= 0.05)
write_csv(sigam, path = '~/Dropbox/CSF metabolites/Significant_correlation_amyloidbeta_Dec7.csv')

sig_amyloid_all=csf_groups1[colnames(csf_groups1)%in% sigam$label]
sig_amyloid_all$`amyloid precursor protein`= csf_groups1$`amyloid precursor protein`
sig_amy_cor=as.matrix(cor(sig_amyloid_all))
sig_amy_dist = dist(sig_amy_cor, method = "e")
sig_amy_clus = hclust(sig_amy_dist, method = "complete")
sig_amy_clus1 = cutree(sig_amy_clus, h=8)
plot(sig_amy_clus)
abline(h=8)
col_breaks = c(seq(-1,0,length=100), # for red
               seq(0,0.8,length=100),  # for yellow
               seq(0.81,1,length=100)) # 

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
require(gplots)
clusterCols <- rainbow(length(unique(sig_amy_clus1))) #assign rainbow colors to each of the cluster
myClusterSideBar <- clusterCols[sig_amy_clus1] #assign it to heatmap clusters below
myheatcol <- rev(redgreen(75)) #assign red/green colors to the heatmap
#plot the heatmap
heatmap.2(sig_amy_cor, density.info = "none", trace = "none", col = myheatcol,  Rowv = as.dendrogram(sig_amy_clus), key = F, dendrogram = "row", scale = "none", RowSideColors= myClusterSideBar, margins = c(10, 15))

##########################################################################################################################################
metnames1 = csf_ctrl[-c(219, 615)]
amyloid1 = csf_ctrl$`amyloid precursor protein`
df2 = data.frame(label = paste(names(csf_ctrl[-c(219, 615)])), Estimate=numeric(614), P.value=numeric(614))
for (i in 1:length(metnames)){
  test = cor.test(amyloid1, metnames1[[i]], method = "p")
  df2$Estimate[i]=test$estimate
  df2$P.value[i]=test$p.value
}
df2=df2[order(df2$P.value),]
df2=df2[order(df2$P.value),]
adjustedp2=p.adjust(df2$P.value, method = "b")
df2=cbind(df2, adjustedp2)
sigam_ctrl=filter(df2, df$adjustedp2 <= 0.05)
write_csv(sigam, path = '~/Dropbox/CSF metabolites/Significant_correlation_amyloidbeta_Dec7.csv')

require(gplots)
clusterCols <- rainbow(length(unique(csf_ctrl_clusters1))) #assign rainbow colors to each of the cluster
myClusterSideBar <- clusterCols[csf_ctrl_clusters1] #assign it to heatmap clusters below
myheatcol <- rev(redgreen(75)) #assign red/green colors to the heatmap
#plot the heatmap
heatmap.2(csf_ctrl_matrix1_cor, density.info = "none", trace = "none", col = myheatcol,  Rowv = as.dendrogram(csf_ctrl1_cluster), key = F, dendrogram = "row", RowSideColors= myClusterSideBar, margins = c(10, 15))

##########################################################################################################################################

heatmap.2(csf_ctrl_matrix1_cor, density.info = "none", trace = "none", col = myheatcol, key = T, dendrogram = "row", RowSideColors= myClusterSideBar, margins = c(10, 15))



plot(csf_ctrl1_cluster)


d <- cut(as.dendrogram(csf_ctrl1_cluster), h=10)

plot(d$lower[[1]])

save.image('~/Dropbox/CSF metabolites/Analysis_CSF_SERUM_Dec6_6pm.Rdata')

mat=as.matrix(sig_amyloid_all)
heatmap.2(mat, trace="none",col=greenred(10), Rowv = NA, dendrogram = "col", Colv = as.dendrogram(sig_amy_clus),density.info = "none", key = T, scale = "none")

###only ctrl groups#####

heatmap.2(csf_ctrl_matrix, trace="none",col=greenred(10), Rowv = NA, dendrogram = "col", Colv = as.dendrogram(csf_ctrl1_cluster),density.info = "none", key = F)
##############################################################serum##############################################################################################
serum_groups1 =serum_groups# filter(serum_groups, serum_groups$newgroups == 2 | serum_groups$newgroups == 3)
serum_groups1$SampleDescription = NULL
serum_groups1$EXTID = NULL
serum_groups1$Group = NULL
serum_groups1$SampleDescription = NULL

which(colnames(serum_groups1) == "amyloid precursor protein") ####get column index
metnames = serum_groups1[-c(1, 1111:1116, 349)]
amyloid = serum_groups1$`amyloid precursor protein`
df = data.frame(label = paste(names(metnames)), Estimate=numeric(1108), P.value=numeric(1108))
for (i in 1:length(metnames)){
  test = cor.test(amyloid, metnames[[i]], method = "p")
  df$Estimate[i]=test$estimate
  df$P.value[i]=test$p.value
}
df_serum=df[order(df$P.value),]
adjustedp1_serum=p.adjust(df_serum$P.value, method = "hochberg")
df=cbind(df_serum, adjustedp1_serum)
df
sigam_serum=filter(df, df$adjustedp1 <= 0.05)
write_csv(sigam_serum, path = '~/Dropbox/serum metabolites/Significant_correlation_amyloidbeta_serum_Jan11.csv')
###############################################compare correlations between csf serum ayloid and serum amyloid################################
sigam_all=merge.data.frame(sigam_serum, sigam_csf, by.x = "label", by.y = "label", all = T)
colnames(sigam_all) = gsub(".x", "_SERUM", colnames(sigam_all))
colnames(sigam_all) = gsub(".y", "_CSF", colnames(sigam_all))
sigam_all
arrange(sigam_all, adjustedp1_CSF, adjustedp1_SERUM)
############################limma ####################
require(Biobase)
require(limma)
require(plyr)
require(dplyr)
require(data.table)
object = new("ExpressionSet", exprs = as.matrix(t(csf_groups[3:617])))
object
design1 = model.matrix(~y)
fit = eBayes(lmFit(object, design1))
topTable(fit, coef = "y", adjust = "BH", p = 40)

require(plyr)
require(dplyr)
###retreive clusters for ctrl csf and narc csf 
csf_narc_K_1[grep("amy", csf_narc_K_1$TargetIDs),]
csf_narc_K_1[grep(3, csf_narc_K_1$CSF_NARCOLEPSY),] %>% write.csv(file = '~/Desktop/NARC_cluster.csv')

csf_ctrl_K_1[grep("amy", csf_ctrl_K_1$TargetIDs),]
csf_ctrl_K_1[grep(3, csf_ctrl_K_1$CSF_CTRLS),] %>% write.csv(file = '~/Desktop/CTRL_cluster.csv')


##get all serum and csf protein names from the tow files
serumnames=select(serum1, TargetFullName, Target)
csfnames = select(csf1, TargetFullName, Target)

##generate boxplots for checking normalilty 
serum_melt=melt(serum_groups1, id.vars = c("SampleId","DbID","Dx","0602","HCRT","oldgroups","newgroups"))
ggplot(serum_melt, aes(factor(variable), value))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
csf_melt=melt(csf_groups, id.vars = c("SampleDescription", "EXTID","DbID","Dx","0602","HCRT","oldgroups","newgroups"))
ggplot(csf_melt, aes(factor(variable), value))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))

###back convert log2 to raw values for CSF_group2
csf_rawEM = cbind.data.frame(csf_groups1[c(616:621)], csf_groups1[-c(616:621)])
csfraw_melt=melt(csf_rawEM, id.vars = c(names(csf_groups1[c(616:621)])))
csf_dilEM=select(csf1, Dilution, Target)  
csfraw_melt = as_data_frame(csfraw_melt)
csfraw_dilEM_melt=merge.data.frame(csfraw_melt, csf_dilEM, by.x = "variable", by.y = "Target", all.x = T)
csfraw_dilEM_melt$Dilution[is.na(csfraw_dilEM_melt$Dilution)] = 15
setDT(csfraw_dilEM_melt)
csf_agg=csfraw_dilEM_melt[, list(value = median(value)), by =list(variable, Dilution)]
##read in the japanese soma data
require(readr)
jap=read_delim('SOMA_data.txt', delim = "\t", skip = 6)

###get dilutions for each csf protein from the japanese cohort
jap_dil=read_delim('SOMA_data.txt', delim = "\t", skip = 3)
dil_CSFjap =jap_dil[1:3,]
dil_CSFjap2=data.frame(t(dil_CSFjap[-c(1:3)]))
colnames(dil_CSFjap2)=c("Dilution", "Units", "Target")
dil_CSFjap2=dil_CSFjap2[-1,]
dil_CSFjap2$Dilution=gsub("%", "", dil_CSFjap2$Dilution) %>% as.numeric() %>% round(digits =3)  
setDT(jap)
jap_melt=melt(jap, id.vars = c("Sample", "1=men; 2=women", "years", "Target")) 
jap_dil_melt = merge.data.frame(jap_melt, dil_CSFjap2, by.x = "variable", by.y = "Target", all.x = T)
jap_dil_melt$Dilution[is.na(jap_dil_melt$Dilution)] = 0.05
setDT(jap_dil_melt)
jap_csf_agg=jap_dil_melt[, list(value = median(value)), by = list(variable, Dilution)]
############################################################################################ merge the aggregated data between japanese and our cohort and make plots######################################33
csf_agg_jap_ours=merge.data.frame(jap_csf_agg, csf_agg, by.x = "variable", by.y = "variable", suffixes = c("_JAP", "_OURS"))
##merge in the kmeans clusters 
csf_agg_jap_oursk=merge.data.frame(csf_agg_jap_ours, csf_total_K, by.x = "variable", by.y = "TargetIDs")
require(ggplot2)
ggplot(csf_agg_jap_oursk, aes(log2(value_JAP), value_OURS, size = factor(Dilution_JAP), label = variable))+geom_point(alpha = 0.3)+geom_text_repel(size = 2.5)+facet_wrap(~TOTAL_BOTH, scales = "free")
ggplot(csf_agg_jap_oursk, aes(log2(value_JAP), value_OURS, size = factor(Dilution_JAP), label = variable, color = factor(TOTAL_BOTH)))+geom_point(alpha = 0.3)+facet_wrap(~TOTAL_BOTH, scales = "free")

ggplot(csf_agg_jap_ours, aes(value_JAP, value_OURS, size = Dilution_JAP, label = variable))+geom_text_repel()+facet_wrap(~TOTAL_BOTH, scales = "free")

############################################################################################
jap$totalAPO = (jap$`Apo E2`+jap$`Apo E3`+jap$`Apo E4`)/3
metnames = log2(jap[-c(1:4, 330)])
amyloid =log2(jap$`amyloid precursor protein`)
df = data.frame(label = paste(names(metnames)), Estimate=numeric(length(metnames)), P.value=numeric(length(metnames)))
for (i in 1:length(metnames)){
  test = cor.test(amyloid, metnames[[i]], method = "p")
  df$Estimate[i]=test$estimate
  df$P.value[i]=test$p.value
}
df_jap=df[order(df$P.value),]
adjustedp1_jap=p.adjust(df_jap$P.value, method = "hochberg")
df=cbind(df_jap, adjustedp1_jap)
df
require(plyr);require(dplyr)
sigam_jap=df#filter(df, df$adjustedp1_jap <= 0.05)
sigam_jap

#############correlation matrix of japanese data###############################################
jap_cor = as.matrix(cor(jap[-c(1:4)]))
jap_dist = dist(jap_cor, method = "e")
jap_clust = hclust(jap_dist, method = 'complete')
plot(jap_clust)
d_jap <- cut(as.dendrogram(jap_clust), h=10)
plot(d_jap$lower[[9]])

require(HSAUR3)
km    <- kmeans(scale(jap[-c(1:4)]), centers = 5)
dissE <- daisy(jap[-c(1:4)]) 
dE2   <- dissE^2
sk2   <- silhouette(km$cl, dE2)
plot(sk2)

############################################SAMR analysis japanese data#####################################
require(samr)
#y = ifelse(jap$`1=men; 2=women` == "1", 1, 2)
y=jap$years
x = t(scale(jap[-c(1:4)]))
d = list(x=x, y=y, geneid=rownames(x), logged2 =T)

samr.obj<-samr(d,  resp.type="Quantitative", assay.type = "array", center.arrays = T, nperms = 1000)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.8
samr.plot(samr.obj,delta)

### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
################################################Limma analysis ###############################################


object = new("ExpressionSet", exprs = as.matrix(t(log2(jap[-c(1:4)]))))
object
design1 = model.matrix(~`1=men; 2=women`+years, data = jap)
fit = eBayes(lmFit(object, design1))
topTable(fit, coef = "`1=men; 2=women`", adjust = "BH", number = 50)
topTable(fit, coef = "years", adjust = "BH", number = 50)
volcanoplot(fit, coef = "`1=men; 2=women`",highlight=20, pch=20)






############################################# merge the corre estimates between japan and ours ###
csf_merge_JO = merge.data.frame(sigam_jap, sigam_csf, by.x = "label", by.y = "label", suffixes = c("_JAP", "_OURS"))
################add in the dilution factors for the japanese
dil_jP=select(csf_agg_jap_oursk, Dilution_JAP, variable)
########################merge the dilutions for each protein########################################
csf_merge_JO2 = merge.data.frame(csf_merge_JO, dil_jP, by.x = "label", by.y = "variable")

ggplot(csf_merge_JO, aes(Estimate_JAP, Estimate_OURS, label = label))+geom_text(size = 3.5, angle = 90, check_overlap = T)#+geom_point(size = 5)+geom_text(label = csf_merge_JO$label,check_overlap = T)
ggplot(csf_merge_JO2, aes(Estimate_JAP, Estimate_OURS, label = label, size = factor(Dilution_JAP)))+geom_point(alpha = 0.5)+geom_text(size = 3.5, angle = 45, check_overlap = T)
################################################################################################################################
#write out the final files for fresh analsysis 
###from IDStoling get all the clinical data
clinical_data=select(ids_ling, DbID, `ID to Sharhad (ST)`, Gender, Age, Race, BMI, `Interval from sleepiness(year)`, `Interval from cataplexy(year)`)
clinical_data = merge.data.frame(revised_groups[-c(1, 3:5)], clinical_data, by.x = "EMID", by.y = "ID to Sharhad (ST)")
###make final files #### for further analysis
csf_final = merge.data.frame(clinical_data, csf_groups[-c(2, 618:623)], by.x = "EMID", by.y = "SampleDescription")
write_csv(csf_final, path = '~/Dropbox/CSF metabolites/CSFProteinAptamers/CSF_FINAL_Jan18.csv')
serum_final = merge.data.frame(clinical_data, serum_groups[-c(2, 1112:1117)], by.x = "EMID", by.y = "SampleId")
write_csv(serum_final, path = '~/Dropbox/CSF metabolites/CSFProteinAptamers/SERUM_FINAL_Jan18.csv')
write_csv(jap, path = '~/Dropbox/CSF metabolites/CSFProteinAptamers/JAPCSF_FINAL_Jan18.csv')
##########merge final 
Combined_final = merge.data.frame(clinical_data, merg1, by.x = "EMID", by.y = "SampleDescription_CSF")
require(readr)
write_csv(Combined_final, path = '~/Dropbox/CSF metabolites/CSFProteinAptamers/combined_FINAL_Jan18.csv')


















