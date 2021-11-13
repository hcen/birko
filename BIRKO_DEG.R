library(tidyverse)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
# package for changing EsembleID to gene name
#BiocManager::install("EnsDb.Mmusculus.v79")
#library(EnsDb.Mmusculus.v79)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set work dir to current folder
getwd()

raw.counts = read.table(file="input/kallisto_gene_counts.counts.matrix", row.names=1,check.names=FALSE)
view(raw.counts)
raw.counts <- round(raw.counts,0)


sampleID <- colnames(raw.counts)
sampleID
genotype <- c("HET","HET","HET","HET","HET","HET","HET","KO","KO","KO","KO","KO","KO","KO","KO","KO","WT","WT","WT","WT","WT","WT",
              "HET","HET","HET","HET","HET","KO","KO","KO","KO","KO","WT","WT","WT","WT","WT","WT")
sex <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F",
         "M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M")
sex.genotype <- paste0(sex,sep=".",genotype)
meta.data <- data.frame(sampleID,genotype,sex,sex.genotype)
rownames(meta.data) <- meta.data$sampleID
view(meta.data)
write.csv(meta.data, sep="\t",file="input/metadata.csv", 
          row.names=TRUE,col.names=NA,quote=FALSE)
meta.data <-  read.csv(file="input/metadata.csv", row.names = 1)

# Exclude Gcg CPM> 10000, F_WT_30A2	F_WT_34B2	M_KO_374 M_WT_32A4. 
meta.data1 <- meta.data %>% filter(!sampleID %in% c("F_WT_30A2","F_WT_34B2","M_KO_374","M_WT_32A4"))
raw.counts1 <- raw.counts %>% dplyr::select(meta.data1$sampleID)
view(meta.data1)
view(raw.counts1)
# Exclude Gcg CPM> 10000, F_WT_30A2	F_WT_34B2	M_KO_374 M_WT_32A4. Also exclude F_KO_38A3 because it didn't pass QC 
meta.data2 <- meta.data %>% filter(!sampleID %in% c("F_KO_38A3","F_WT_30A2","F_WT_34B2","M_KO_374","M_WT_32A4"))
raw.counts2 <- raw.counts %>% dplyr::select(meta.data2$sampleID)
view(meta.data2)
view(raw.counts2)

# Also exclude Gcg CPM> 100 in male and female KO, F_KO_23B1, F_KO_32A1, F_KO_32C1
meta.data3 <- meta.data %>% filter(!sampleID %in% c("F_KO_38A3","F_WT_30A2","F_WT_34B2","M_KO_374","M_WT_32A4",
                                                    "F_KO_23B1", "F_KO_32A1", "F_KO_32C1")) %>%
  filter(sex.genotype %in% c("F.KO","M.KO"))
raw.counts3 <- raw.counts %>% dplyr::select(meta.data3$sampleID)
view(meta.data3)
view(raw.counts3)

#
meta.data.s <- read.csv(file="input/metadata.s.csv",row.names = 1)
view(meta.data.s)
meta.data.g <- read.csv(file="input/metadata.g.csv",row.names = 1)
view(meta.data.g)
meta.data.s <- meta.data.s %>% filter(!sampleID %in% c("F_KO_38A3","F_WT_30A2","F_WT_34B2","M_KO_374","M_WT_32A4"))
meta.data.g <- meta.data.g %>% filter(!sampleID %in% c("F_KO_38A3","F_WT_30A2","F_WT_34B2","M_KO_374","M_WT_32A4"))

### create DESeq matrix
count.data.set = DESeqDataSetFromMatrix(countData=raw.counts1, 
                                        colData=meta.data1, design= ~ genotype) 
count.data.set = DESeqDataSetFromMatrix(countData=raw.counts1, 
                                        colData=meta.data1, design= ~ sex.genotype) 

count.data.set = DESeqDataSetFromMatrix(countData=raw.counts2, 
                                        colData=meta.data2, design= ~ genotype) 
count.data.set = DESeqDataSetFromMatrix(countData=raw.counts2, 
                                        colData=meta.data2, design= ~ sex.genotype) 

count.data.set = DESeqDataSetFromMatrix(countData=raw.counts3, 
                                        colData=meta.data3, design= ~ sex.genotype) 
# Filter low count
nrow(count.data.set)
keep <- rowSums(counts(count.data.set)>5) >=10 # genes counts more than 5 in at least 3 samples
keep <- rowSums(counts(count.data.set)>5) >=4 # genes counts more than 5 in at least 3 samples
count.filter <- count.data.set[keep,]
nrow(count.filter)

# create DESeq object
#count.data.set.object <- DESeq(count.data.set)
count.data.set.object <- DESeq(count.filter)
#count.data.set.object
# 'vst' normalization (varianceStabilizingTransformation)
vsd <- vst(count.data.set.object)

### extract normalized counts
norm.data = assay(vsd)
head(norm.data)
dim(norm.data)
write.table(norm.data, sep="\t",file="data/Norm_data_noAlpha-1.txt", 
            row.names=TRUE,col.names=NA,quote=FALSE)
### hierarchical clustering analyses and to plot a dendrogram. Evaluate dissimilarities (calculate Euclidean distance) between all eight replicates based on their normalized gene counts.
sampleDists <- dist(t(norm.data),  method = "euclidean")

### Having the distance (dissimilarity) we can finally perform hierarchical cluster analysis using hclust function
clusters=hclust(sampleDists)
plot(clusters)

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
#rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pDist<-pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(pDist, "pDist.png")
dev.off()
### plots PCA for the first two principal components
plotPCA(vsd, intgroup=c("sex.genotype")) +
  theme_bw()
plotPCA(vsd, intgroup=c("genotype")) +
  theme_bw()

#=================================================
# extra
getMethod("plotPCA","DESeqTransform")

plotPCA.pc23 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[2:3]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC2", y = "PC3", 
                                color = "group")) + geom_point(size = 3) + 
      xlab(paste0("PC2: ", round(percentVar[2] * 
                                   100), "% variance")) + ylab(paste0("PC3: ", 
                                                                      round(percentVar[3] * 100), "% variance")) + 
      coord_fixed()
  }
  .local(object, ...)
}
plotPCA.pc23(vsd,intgroup=c("group"))+ # 4x3.5in
  theme_bw()
##
plotPCA.pc34 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[3:4]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", 
                                color = "group")) + geom_point(size = 3) + 
      xlab(paste0("PC3: ", round(percentVar[3] * 
                                   100), "% variance")) + ylab(paste0("PC4: ", 
                                                                      round(percentVar[4] * 100), "% variance")) + 
      coord_fixed()
  }
  .local(object, ...)
}
plotPCA.pc34(vsd,intgroup=c("group"))+ # 4x3.5in
  theme_bw()

#=========================================================
### DEG
##What cooksCutoff filtering is... This is simply a cutoff for excluding genes that have higher variability of counts between replicates than between the compared condition. 
##Another test that is performed is excluding genes that have too low counts relative to the mean counts of the experiment. This test is called independentFiltering.
## Order the results based on increasing padj value
res <- results(count.data.set.object, 
               contrast=c("genotype","KO","WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("genotype","HET","WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("genotype","KO","HET"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("sex.genotype","M.KO","M.WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("sex.genotype","M.HET","M.WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("sex.genotype","M.KO","M.HET"),
               alpha = 0.05)

res <- results(count.data.set.object, 
               contrast=c("sex.genotype","F.KO","F.WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("sex.genotype","F.HET","F.WT"),
               alpha = 0.05)
res <- results(count.data.set.object, 
               contrast=c("sex.genotype","F.KO","F.HET"),
               alpha = 0.05)


res <- results(count.data.set.object, 
               contrast=c("sex.genotype","F.KO","M.KO"),
               alpha = 0.05)

summary(res)
out <- capture.output(summary(res))
##
##
#cat("result summary", out, file="data/res_summary_noAlpha_KO_WT.txt", sep="\n", append=TRUE)
#cat("result summary", out, file="data/res_summary_noAlpha_KO_WT_male.txt", sep="\n", append=TRUE)
#cat("result summary", out, file="data/res_summary_noAlpha_KO_WT_female.txt", sep="\n", append=TRUE)

cat("result summary", out, file="data/res_summary_noAlpha-1_KO_WT.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_HET_WT.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_KO_HET.txt", sep="\n", append=TRUE)

cat("result summary", out, file="data/res_summary_noAlpha-1_KO_WT_male.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_HET_WT_male.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_KO_HET_male.txt", sep="\n", append=TRUE)

cat("result summary", out, file="data/res_summary_noAlpha-1_KO_WT_female.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_HET_WT_female.txt", sep="\n", append=TRUE)
cat("result summary", out, file="data/res_summary_noAlpha-1_KO_HET_female.txt", sep="\n", append=TRUE)


cat("result summary", out, file="data/res_summary_noAlpha100_KO_female_male.txt", sep="\n", append=TRUE)
#
res = na.omit(res)
res.ord = res[order(res$padj),]
res.ord <- as.data.frame(res.ord)
##
##
#write.table(res.ord, sep="\t",file="data/Results_noAlpha_KO_WT.txt", row.names=TRUE,col.names=NA,quote=FALSE)
#write.table(res.ord, sep="\t",file="data/Results_noAlpha_KO_WT_male.txt", row.names=TRUE,col.names=NA,quote=FALSE)
#write.table(res.ord, sep="\t",file="data/Results_noAlpha_KO_WT_female.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_WT.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_HET_WT.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_HET.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_WT_male.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_HET_WT_male.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_HET_male.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_WT_female.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_HET_WT_female.txt", row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.ord, sep="\t",file="data/Results_noAlpha-1_KO_HET_female.txt", row.names=TRUE,col.names=NA,quote=FALSE)


write.table(res.ord, sep="\t",file="data/Results_noAlpha100_KO_female_male.txt", row.names=TRUE,col.names=NA,quote=FALSE)

#res.ord=read.table(file="data/Results_noAlpha_KO_WT.txt",row.names = 1)
#res.ord=read.table(file="data/Results_noAlpha_KO_WT_male.txt",row.names = 1)
#res.ord=read.table(file="data/Results_noAlpha_KO_WT_female.txt",row.names = 1)

res.sig = res.ord[res.ord$padj <= 0.05,]

res.sig$ensembl <- gsub("\\..*","",rownames(res.sig))
res.sig$entrez = mapIds(org.Mm.eg.db, keys=res.sig$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res.sig$symbol =mapIds(org.Mm.eg.db, keys=res.sig$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res.sig$description =mapIds(org.Mm.eg.db, keys=res.sig$ensembl, column="GENENAME", keytype="ENSEMBL", multiVals="first")
res.sig$description <- gsub(x = res.sig$description, pattern = "\\ ",replacement = "_") 
res.sig$description <- gsub(x = res.sig$description, pattern = "\\,",replacement = ".")
# KO vs WT, both sex: ENSMUSG00000095403, Gm21092, cadherin 11 pseudogen
#res.sig[res.sig$ensembl=="ENSMUSG00000095403",colnames(res.sig)=="symbol"] <- "Gm21092"
#res.sig[res.sig$ensembl=="ENSMUSG00000095403",colnames(res.sig)=="description"] <- "cadherin 11 pseudogen"

# add gene name for KO female vs male
res.sig[res.sig$ensembl=="ENSMUSG00000098078",colnames(res.sig)=="symbol"] <- "Gm26992"

View(res.sig)

dup=res.sig %>% group_by(symbol) %>% dplyr::filter(n() > 1)
View(dup)

write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_WT.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_HET_WT.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_HET.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)

write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_WT_male.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_HET_WT_male.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_HET_male.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)

write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_WT_female.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_HET_WT_female.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(res.sig,file="data/Results_sig_noAlpha-1_KO_HET_female.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)


write.table(res.sig,file="data/Results_sig_noAlpha100_KO_female_male.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
#
sig.KO.WT <- read.table(file="data/Results_sig_noAlpha-1_KO_WT.txt", row.names=1)
sig.HET.WT <- read.table(file="data/Results_sig_noAlpha-1_HET_WT.txt",  row.names=1)
sig.KO.HET <- read.table(file="data/Results_sig_noAlpha-1_KO_HET.txt",  row.names=1)

sig.KO.WT.m <- read.table(file="data/Results_sig_noAlpha-1_KO_WT_male.txt",  row.names=1)
sig.HET.WT.m <- read.table(file="data/Results_sig_noAlpha-1_HET_WT_male.txt",  row.names=1)
sig.KO.HET.m <- read.table(file="data/Results_sig_noAlpha-1_KO_HET_male.txt",  row.names=1)

sig.KO.WT.f <- read.table(file="data/Results_sig_noAlpha-1_KO_WT_female.txt",  row.names=1)
sig.HET.WT.f <- read.table(file="data/Results_sig_noAlpha-1_HET_WT_female.txt",  row.names=1)
sig.KO.HET.f <- read.table(file="data/Results_sig_noAlpha-1_KO_HET_female.txt",  row.names=1)

sig <- Reduce(rbind, list(sig.KO.WT,sig.HET.WT,sig.KO.HET,sig.KO.WT.m,sig.HET.WT.m,sig.KO.HET.m,sig.KO.WT.f,sig.HET.WT.f,sig.KO.HET.f))
sig1 <- sig[!duplicated(sig$ensembl),]
view(sig1)

norm.data <- read.table(file="data/Norm_data_noAlpha-1.txt", row.names=1)
norm.data$ensembl <- gsub("\\..*","",rownames(norm.data))
#
#
head(norm.data)
norm.sig <- norm.data %>% filter(norm.data$ensembl %in% sig1$ensembl)

colnames(sig.KO.WT) <- paste(colnames(sig.KO.WT),"KO.WT", sep = "_")
colnames(sig.HET.WT) <- paste(colnames(sig.HET.WT),"HET.WT", sep = "_")
colnames(sig.KO.HET) <- paste(colnames(sig.KO.HET),"KO.HET", sep = "_")
colnames(sig.KO.WT.m) <- paste(colnames(sig.KO.WT.m),"KO.WT.Male", sep = "_")
colnames(sig.HET.WT.m) <- paste(colnames(sig.HET.WT.m),"HET.WT.Male", sep = "_")
colnames(sig.KO.HET.m) <- paste(colnames(sig.KO.HET.m),"KO.HET.Male", sep = "_")
colnames(sig.KO.WT.f) <- paste(colnames(sig.KO.WT.f),"KO.WT.Female", sep = "_")
colnames(sig.HET.WT.f) <- paste(colnames(sig.HET.WT.f),"HET.WT.Female", sep = "_")
colnames(sig.KO.HET.f) <- paste(colnames(sig.KO.HET.f),"KO.HET.Female", sep = "_")
view(sig.KO.WT)

#====
IDorder <- match(meta.data.s$sampleID,colnames(norm.sig)) # https://hbctraining.github.io/Intro-to-R/lessons/06_matching_reordering.html
IDorder
norm.sig <- norm.sig[,IDorder]
all(colnames(norm.sig) == meta.data.s$sampleID)
#===
IDorder <- match(meta.data.g$sampleID,colnames(norm.sig)) 
IDorder
norm.sig <- norm.sig[,IDorder]
all(colnames(norm.sig) == meta.data.g$sampleID)
#====

norm.sig$ensembl <- gsub("\\..*","",rownames(norm.sig))
head(norm.sig)
norm.sig.p<- norm.sig %>% 
  left_join(sig.KO.WT[,c(7,6,2)], by =c("ensembl"="ensembl_KO.WT")) %>% 
  left_join(sig.HET.WT[,c(7,6,2)], by =c("ensembl"="ensembl_HET.WT"))%>% 
  left_join(sig.KO.HET[,c(7,6,2)], by =c("ensembl"="ensembl_KO.HET"))%>% 
  left_join(sig.KO.WT.m[,c(7,6,2)], by =c("ensembl"="ensembl_KO.WT.Male")) %>% 
  left_join(sig.HET.WT.m[,c(7,6,2)], by =c("ensembl"="ensembl_HET.WT.Male"))%>% 
  left_join(sig.KO.HET.m[,c(7,6,2)], by =c("ensembl"="ensembl_KO.HET.Male"))%>% 
  left_join(sig.KO.WT.f[,c(7,6,2)], by =c("ensembl"="ensembl_KO.WT.Female")) %>% 
  left_join(sig.HET.WT.f[,c(7,6,2)], by =c("ensembl"="ensembl_HET.WT.Female"))%>% 
  left_join(sig.KO.HET.f[,c(7,6,2)], by =c("ensembl"="ensembl_KO.HET.Female"))

norm.sig.p$entrez =mapIds(org.Mm.eg.db, keys=norm.sig.p$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
norm.sig.p$symbol =mapIds(org.Mm.eg.db, keys=norm.sig.p$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

norm.sig.p$description =mapIds(org.Mm.eg.db, keys=norm.sig.p$ensembl, column="GENENAME", keytype="ENSEMBL", multiVals="first")
norm.sig.p$description <- gsub(x = norm.sig.p$description, pattern = "\\ ",replacement = "_") 
norm.sig.p$description <- gsub(x = norm.sig.p$description, pattern = "\\,",replacement = ".")
dup=norm.sig.p %>% group_by(ensembl) %>% dplyr::filter(n() > 1) 
head(dup)
write.table(norm.sig.p,file="data/norm.sig.s.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(norm.sig.p,file="data/norm.sig.g.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)


#$$$$$$$$$$$$$$$$$$$$$
norm.sig.list <- norm.sig.p %>%  mutate(logFC = coalesce(log2FoldChange_KO.WT,log2FoldChange_HET.WT,log2FoldChange_KO.HET,
                                                         log2FoldChange_KO.WT.Male,log2FoldChange_HET.WT.Male,log2FoldChange_KO.HET.Male,
                                                         log2FoldChange_KO.WT.Female,log2FoldChange_HET.WT.Female,log2FoldChange_KO.HET.Female))
view(norm.sig.list)
write.table(norm.sig.list,file="data/norm.sig.list.txt", sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
#$$$$$$$$$$$$$$$$$$$$$

norm.sig.p <- read.table(file="data/norm.sig.s.txt",row.names = 1)
norm.sig.p <- read.table(file="data/norm.sig.g.txt",row.names = 1)

norm.sig.p <- norm.sig.p %>% dplyr::filter(!symbol %in% c(NA)) 
view(norm.sig.p)
#norm.sig.p <- norm.sig.p[,c(25:33,13:24,1:12,34:45)] # genotype together format, reorder WT-HET-KO
# Heatmap
norm.sig.p.noH <- norm.sig.p[!is.na(norm.sig.p$padj_KO.WT)| !is.na(norm.sig.p$padj_KO.WT.Male) | !is.na(norm.sig.p$padj_KO.WT.Female),
                             c(1:8,16:23,29:33,35,38,41,44,45) #sex together format
                             #c(1:12,25:33,35,38,41,44,45) # genotype together format
                             ]
view(norm.sig.p.noH)

m <- as.matrix(as.data.frame(lapply(norm.sig.p[,c(1:33)], as.numeric),check.names=F))
m <- as.matrix(as.data.frame(lapply(norm.sig.p.noH[,c(1:21)], as.numeric),check.names=F)) #sex together format
m
library(gplots)
par(oma=c(2,3,2,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m,
                 labRow = norm.sig.p$symbol,
                 #labRow = norm.sig.p.noH$symbol,
                 #labCol = NA,
                 scale = "row", 
                 col=my_palette, 
                 trace = "none", 
                 density.info = "none",
                 cexRow = 0.9,
                 cexCol = 0.8,
                 offsetRow = -0.2,
                 offsetCol = 0,
                 distfun = dist,
                 #hclustfun = hclust,
                 hclustfun=hclust.m,
                 Colv=FALSE,
                 dendrogram='row',
                 key=TRUE, keysize=0.75, key.title = NA,key.xlab=NA,cex.lab=5.0, cex.axis=1.0,
                 #key.par=list(mar=c(1,1,1,1)),
                 lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
                 key.par = list(mar=c(5,1,0,8),cex=0.6),
                 key.xtickfun=function() {
                   cex <- par("cex")*par("cex.axis")
                   side <- 1
                   line <- 0
                   col <- par("col.axis")
                   font <- par("font.axis")
                   mtext("low", side=side, at=0, adj=0,
                         line=line, cex=cex, col=col, font=font)
                   mtext("high", side=side, at=1, adj=1,
                         line=line, cex=cex, col=col, font=font)
                   return(list(labels=FALSE, tick=FALSE))
                 },
                 #colsep=c(8,15,19,23,28)
                 colsep=c(4,9,16,21,29) #genotype together, wt left
                 #colsep=c(8,12,16)
) # all 5x11 in
dev.off()
out$rowInd
norm.sig.hm <- norm.sig.p[out$rowInd,] %>% map_df(rev)
norm.sig.hm <- norm.sig.p.noH[out$rowInd,] %>% map_df(rev)
view(norm.sig.hm)

write.table(norm.sig.hm,file="data/heatmap-all-s.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-all-g.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-noHet-g.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-noHet-s.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 

# male
norm.sig.m <- norm.sig.p[,c(1:33,38:40,44,45)] %>% dplyr::filter(!symbol %in% c(NA)) #from genotype together format?
norm.sig.m <- norm.sig.m[!is.na(norm.sig.m$padj_KO.WT.Male) | !is.na(norm.sig.m$padj_HET.WT.Male) | !is.na(norm.sig.m$padj_KO.HET.Male), c(9:12,20:24,29:38)]
view(norm.sig.m)
m <- as.matrix(as.data.frame(lapply(norm.sig.m[,c(1:14)], as.numeric),check.names=F))

norm.sig.m <- norm.sig.p.noH[!is.na(norm.sig.p.noH$padj_KO.WT.Male),c(13:21,23,25,26)]  #from sex together format
view(norm.sig.m)
m <- as.matrix(as.data.frame(lapply(norm.sig.m[,c(1:9)], as.numeric),check.names=F))

library(gplots)
par(oma=c(2,2,3,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m,
                 labRow = norm.sig.m$symbol,
                 #labCol = NA,
                 scale = "row", 
                 col=my_palette, 
                 trace = "none", 
                 density.info = "none",
                 cexRow = 0.9,
                 cexCol = 0.8,
                 offsetRow = -0.2,
                 offsetCol = 0,
                 distfun = dist,
                 #hclustfun = hclust,
                 hclustfun=hclust.m,
                 Colv=FALSE,
                 dendrogram='row',
                 key=TRUE, keysize=0.75, key.title = NA,key.xlab=NA,cex.lab=5.0, cex.axis=1.0,
                 #key.par=list(mar=c(1,1,1,1)),
                 lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
                 key.par = list(mar=c(5,1,0,8),cex=0.6),
                 key.xtickfun=function() {
                   cex <- par("cex")*par("cex.axis")
                   side <- 1
                   line <- 0
                   col <- par("col.axis")
                   font <- par("font.axis")
                   mtext("low", side=side, at=0, adj=0,
                         line=line, cex=cex, col=col, font=font)
                   mtext("high", side=side, at=1, adj=1,
                         line=line, cex=cex, col=col, font=font)
                   return(list(labels=FALSE, tick=FALSE))
                 },
                 #colsep=c(8,15,19,23,28)
                 colsep=c(4,9)
) # male 5x11 in
dev.off()
out$rowInd
norm.sig.hm <- norm.sig.m[out$rowInd,] %>% map_df(rev)
view(norm.sig.hm)

write.table(norm.sig.hm,file="data/heatmap-male.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-noHET-male.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 

# female

norm.sig.f <- norm.sig.p[,c(1:33,41:45)] %>% dplyr::filter(!symbol %in% c(NA))
norm.sig.f <- norm.sig.f[!is.na(norm.sig.f$padj_KO.WT.Female) | !is.na(norm.sig.f$padj_HET.WT.Female) | !is.na(norm.sig.f$padj_KO.HET.Female), c(1:8,13:19,25:28,34:38)]
view(norm.sig.f)
m <- as.matrix(as.data.frame(lapply(norm.sig.f[,c(1:19)], as.numeric),check.names=F))
m

norm.sig.f <- norm.sig.p.noH[!is.na(norm.sig.p.noH$padj_KO.WT.Female),c(1:12,24,25,26)]  #from sex together format
view(norm.sig.f)
m <- as.matrix(as.data.frame(lapply(norm.sig.f[,c(1:12)], as.numeric),check.names=F))

library(gplots)
par(oma=c(2,2,3,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m,
                 labRow = norm.sig.f$symbol,
                 #labCol = NA,
                 scale = "row", 
                 col=my_palette, 
                 trace = "none", 
                 density.info = "none",
                 cexRow = 0.9,
                 cexCol = 0.8,
                 offsetRow = -0.2,
                 offsetCol = 0,
                 distfun = dist,
                 #hclustfun = hclust,
                 hclustfun=hclust.m,
                 Colv=FALSE,
                 dendrogram='row',
                 key=TRUE, keysize=0.75, key.title = NA,key.xlab=NA,cex.lab=5.0, cex.axis=1.0,
                 #key.par=list(mar=c(1,1,1,1)),
                 lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
                 key.par = list(mar=c(5,1,0,8),cex=0.6),
                 key.xtickfun=function() {
                   cex <- par("cex")*par("cex.axis")
                   side <- 1
                   line <- 0
                   col <- par("col.axis")
                   font <- par("font.axis")
                   mtext("low", side=side, at=0, adj=0,
                         line=line, cex=cex, col=col, font=font)
                   mtext("high", side=side, at=1, adj=1,
                         line=line, cex=cex, col=col, font=font)
                   return(list(labels=FALSE, tick=FALSE))
                 },
                 #colsep=c(8,15,19,23,28)
                 colsep=c(8,15)
) # male 5x11 in
dev.off()
out$rowInd
norm.sig.hm <- norm.sig.f[out$rowInd,] %>% map_df(rev)
view(norm.sig.hm)

write.table(norm.sig.hm,file="data/heatmap-female.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-noHet-female.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 

# both sex, mixed
view(norm.sig.p)
norm.sig.b <- norm.sig.p[,c(1:37,44,45)] %>% dplyr::filter(!symbol %in% c(NA))
norm.sig.b <- norm.sig.b[!is.na(norm.sig.b$padj_KO.WT) | !is.na(norm.sig.b$padj_HET.WT) | !is.na(norm.sig.b$padj_KO.HET), ]
view(norm.sig.b)
m <- as.matrix(as.data.frame(lapply(norm.sig.b[,c(1:33)], as.numeric),check.names=F))

view(norm.sig.p.noH)
norm.sig.b <- norm.sig.p.noH[,c(1:22,25,26)] %>% dplyr::filter(!padj_KO.WT %in% c(NA))
view(norm.sig.b)

m <- as.matrix(as.data.frame(lapply(norm.sig.b[,c(1:21)], as.numeric),check.names=F))

library(gplots)
par(oma=c(2,2,3,2))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
hclust.m <- function(x) hclust(x, method="ward.D2")
out <- heatmap.2(m,
                 labRow = norm.sig.b$symbol,
                 #labCol = NA,
                 scale = "row", 
                 col=my_palette, 
                 trace = "none", 
                 density.info = "none",
                 cexRow = 0.9,
                 cexCol = 0.8,
                 offsetRow = -0.2,
                 offsetCol = 0,
                 distfun = dist,
                 #hclustfun = hclust,
                 hclustfun=hclust.m,
                 Colv=FALSE,
                 dendrogram='row',
                 key=TRUE, keysize=0.75, key.title = NA,key.xlab=NA,cex.lab=5.0, cex.axis=1.0,
                 #key.par=list(mar=c(1,1,1,1)),
                 lhei=c(0.1,12,1), lwid=c(5,30),lmat = rbind(c(0,3),c(2,1),c(0,4)),  # https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
                 key.par = list(mar=c(5,1,0,8),cex=0.6),
                 key.xtickfun=function() {
                   cex <- par("cex")*par("cex.axis")
                   side <- 1
                   line <- 0
                   col <- par("col.axis")
                   font <- par("font.axis")
                   mtext("low", side=side, at=0, adj=0,
                         line=line, cex=cex, col=col, font=font)
                   mtext("high", side=side, at=1, adj=1,
                         line=line, cex=cex, col=col, font=font)
                   return(list(labels=FALSE, tick=FALSE))
                 },
                 #colsep=c(8,15,19,23,28)
                 colsep=c(12,24)
) # male 5x11 in
dev.off()
out$rowInd
norm.sig.hm <- norm.sig.b[out$rowInd,] %>% map_df(rev)
view(norm.sig.hm)

write.table(norm.sig.hm,file="data/heatmap-mixedSex.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 
write.table(norm.sig.hm,file="data/heatmap-mixedSex-noHET.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE) 




### Vocano plots using "EnhancedVolcano" package

#read_xlsx(path="data/Results_sig_noAlpha_KO_WT.xlsx")
res.ord=read.table(file="data/Results_noAlpha-1_KO_WT.txt",row.names = 1)
res.ord.m=read.table(file="data/Results_noAlpha-1_KO_WT_male.txt",row.names = 1)
res.ord.f=read.table(file="data/Results_noAlpha-1_KO_WT_female.txt",row.names = 1)

head(res.ord)
is.data.frame(res.ord)
res.ord <- as.data.frame(res.ord)
res.ord.m <- as.data.frame(res.ord.m)
res.ord.f <- as.data.frame(res.ord.f)

res.ord.f$ensembl <- gsub("\\..*","",rownames(res.ord.f))
res.ord.f$entrez = mapIds(org.Mm.eg.db, keys=res.ord.f$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res.ord.f$symbol =mapIds(org.Mm.eg.db, keys=res.ord.f$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
view(res.ord.f)
res.ord.f <- na.omit(res.ord.f)

res.ord.m$ensembl <- gsub("\\..*","",rownames(res.ord.m))
res.ord.m$entrez = mapIds(org.Mm.eg.db, keys=res.ord.m$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res.ord.m$symbol =mapIds(org.Mm.eg.db, keys=res.ord.m$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
view(res.ord.m)
res.ord.m <- na.omit(res.ord.m)

res.ord$ensembl <- gsub("\\..*","",rownames(res.ord))
res.ord$entrez = mapIds(org.Mm.eg.db, keys=res.ord$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res.ord$symbol =mapIds(org.Mm.eg.db, keys=res.ord$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
view(res.ord)
res.ord <- na.omit(res.ord)

# for KO female vs male:
res.ord[res.ord$ensembl=="ENSMUSG00000098078",colnames(res.ord)=="symbol"] <- "Gm26992"

#
EnhancedVolcano(res.ord,
                lab = res.ord$symbol,
                x = "log2FoldChange",
                y = "padj",
                ylab = bquote(~-Log[10]~ 'Adj. P'),
                FCcutoff = 0.01,
                pCutoff = 0.05,
                xlim = c(-2.5, 3),
                ylim = c(-0.2, 3),
                #selectLab = c("Atp5md","Crybb3"),
                pointSize = 3,
                labSize = 4.0,
                labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                col=c('grey20', 'grey20', 'grey20', 'red3'),
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = "KO vs WT",
                legendPosition = 'none')
##
EnhancedVolcano(res.ord.f,
                lab = res.ord.f$symbol,
                x = "log2FoldChange",
                y = "padj",
                ylab = bquote(~-Log[10]~ 'Adj. P'),
                FCcutoff = 0.01,
                pCutoff = 0.05,
                #xlim = c(-2.5, 2.5),
                ylim = c(-0.2, 15),
                #selectLab = c("Atp5md","Crybb3"),
                pointSize = 3,
                labSize = 4.0,
                labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                col=c('grey20', 'grey20', 'grey20', 'red3'),
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = "Female KO vs WT",
                legendPosition = 'none')

#
EnhancedVolcano(res.ord.m,
                lab = res.ord.m$symbol,
                x = "log2FoldChange",
                y = "padj",
                ylab = bquote(~-Log[10]~ 'Adj. P'),
                FCcutoff = 0.01,
                pCutoff = 0.05,
                xlim = c(-10, 8),
                ylim = c(-0.2, 5),
                #selectLab = c("Atp5md","Crybb3"),
                pointSize = 3,
                labSize = 4.0,
                labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                col=c('grey20', 'grey20', 'grey20', 'red3'),
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = "Male KO vs WT",
                legendPosition = 'none')

#
EnhancedVolcano(res.ord,
                lab = res.ord$symbol,
                x = "log2FoldChange",
                y = "padj",
                ylab = bquote(~-Log[10]~ 'Adj. P'),
                FCcutoff = 0.01,
                pCutoff = 0.05,
                xlim = c(-10, 25),
                #ylim = c(-0.2, 5),
                ylim = c(-1, 20),
                #selectLab = c("Atp5md","Crybb3"),
                pointSize = 3,
                labSize = 5.0,
                labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                col=c('grey20', 'grey20', 'grey20', 'red3'),
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = "KO Female vs Male",
                legendPosition = 'none')
library("EnhancedVolcano")
#==================================================================
