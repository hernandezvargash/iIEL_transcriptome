
# title -------------------------------------------------------------------

# Transcriptome analysis of mouse gamma-delta T cells of different origin

rm(list=ls())


# libraries ---------------------------------------------------------------


setwd()
dir()

suppressPackageStartupMessages({
  library(oligo)
  library(stringr)
  library(limma)
  library(qqman)
  library(pd.mogene.2.0.st)
  library(mogene20sttranscriptcluster.db)
  library(annotate)
  library(xlsx)
  library(viridis)                  
  library(NMF)
  library(RColorBrewer)
  library(ggpubr)
  library(ggfortify)
  library(org.Mm.eg.db)
  library(stringr)
  library(eulerr)
  library(pals)
  library(ggpubr)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(pathview)
  library(ReactomePA)
  library(msigdbr)
  library(fgsea)
  library(pathfindR)
  library(EnhancedVolcano)
  
})



# loading -----------------------------------------------------------------


filenames <- list.files("~/raw", pattern="CEL", recursive = T, full.names = T)
cel<-read.celfiles(filenames)

ID <- str_sub(filenames, 46, -5)
group<-c(rep("IEL",3), rep("CD44lo_CD27+",3), rep("T1", 2), rep("T17",3), "T1")
pData<-pData(cel)
pData$group <- group
pData$ID <- ID
head(pData)

group.colors <- c("#BABABA", "darkorange", "#1F78B4", "#66A61E")

#hist(cel,which="all",main="raw data")
par(mar=c(12,5,5,5), mfrow=c(2,3))
boxplot(cel, which="all", main="raw data", las=2)
mdsPlot(exprs(cel),numPositions=1000,sampGroups=pData$group,
        legendPos="top",legendNCol=1, #sampNames=pData$group,
        pal = group.colors, pch = 16,
        main="MDS after RMA
        1000 most variable positions")
plotSampleRelation(exprs(cel), subset = 10000, method = "cluster",
                   main= "Raw data
                   Sample relations based on 10000 selected genes")

eset<-rma(cel)
exprs<-exprs(eset)
head(exprs)
dim(exprs) # 41345 genes and 12 samples

# save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="data.txt")
pData(eset)<-pData
save(eset,file="eset.RData")



# inspection --------------------------------------------------------------


#MAplot(eset[,1:4],cex=1)
#boxplot(eset,main="after RMA")
#hist(eset,main="after RMA")
pData <- pData(eset)
boxplot(exprs(eset), which="all", main="RMA data", las=2)

par(mfrow=c(1,1), mar=c(4,4,4,4))
mdsPlot(exprs(eset),numPositions=1000,sampGroups=pData$group,
        legendPos="bottomright",legendNCol=1, #sampNames=pData$group,
        pal = group.colors, pch = 16,
        main="MDS after RMA
        1000 most variable positions")
plotSampleRelation(exprs(eset), subset = 10000, method = "cluster",
                   main= "After RMA
                   Sample relations based on 10000 selected genes")


# PCA

t.counts <- t(exprs)
t.counts[1:5, 1:10]

pca_res <- prcomp(t.counts, scale. = F)
autoplot(pca_res, data = pData, colour = 'group', size = 3) +
  theme_minimal() +
  scale_color_manual(values = group.colors) +
  ggtitle("PCA after normalization")



# differential expression -------------------------------------------------

rm(list=ls())

load("eset.RData")

par(mar=c(5,5,5,5), mfrow=c(1,3))
exprs <- exprs(eset)
pData <- pData(eset)


# Model 1: IEL pairwise comparisons

mod<-model.matrix(~0+as.factor(group),data=pData)
fit = lmFit(exprs,mod)
colnames(fit)

contrast.matrix<-c(1,-1,0,0) # Naive vs IEL
contrast.matrix<-c(0,-1,0,1) # T17 vs IEL
contrast.matrix<-c(0,-1,1,0) # Type1 vs IEL

fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
#top<-topTable(eb, adjust="BH",number=100000, p.value=0.05, lfc=1)
#head(top)
#lmPvals = sample(eb$p.value,10000); qq(lmPvals)
#chisq <- qchisq(1-eb$p.value,1) # Calculating ChiSq values
#lmbd=median(chisq)/qchisq(0.5,1) # Calculating Lambda
#mtext(paste("model= ~as.factor(group)"), side=3, line=-2, at=1, cex=1)
#mtext(expression(paste(lambda,"=")), side=3, line=-3, at=0.5, cex=1)
#mtext(paste(round(lmbd, digits=2)), side=3, line=-3, at=1, cex=1)

# volcano
top<-topTable(eb, adjust="BH",number=100000, p.value=1, lfc=0)
volcanoData <- cbind(top$logFC, -log10(top$adj.P.Val))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch= 20, col = ifelse((volcanoData[, 2] > 1.3) & (volcanoData[, 1] > 1) ,'red',
                                        ifelse((volcanoData[, 2] > 1.3) & (volcanoData[, 1] < -1) ,'blue','gray')),
#     main = "CD44lo_CD27+ vs IEL")
#     main = "T17 vs IEL")
main = "T1 vs IEL")


write.csv(top, "Naive.like_vs_IEL_all.csv")
write.csv(top, "T17_vs_IEL_all.csv")
write.csv(top, "T1_vs_IEL_all.csv")

top<-topTable(eb, adjust="BH",number=100000, p.value=0.05, lfc=0)

# Load annotation library
annodb <- "mogene20sttranscriptcluster.db"
ID<-rownames(top) # for the significant probes
#ID     <- featureNames(eset) # for all the dataset
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))
#my_frame <- data.frame(exprs(eset))
top$Symbol<-Symbol
top$Name<-Name
top$Entrez<-Entrez
head(top)

# Write out to a file:
write.xlsx(top,file="Naive.like_vs_IEL.xlsx")
write.csv(top,file="Naive.like_vs_IEL.csv")

write.xlsx(top,file="T17_vs_IEL.xlsx")
write.csv(top,file="T17_vs_IEL.csv")

write.xlsx(top,file="Type1_vs_IEL.xlsx")
write.csv(top,file="Type1_vs_IEL.csv")


# Model 2: each vs all comparisons

group <- factor(pData$group)

group2 <- ifelse(group=="IEL", paste(group), "All_Others")
group2 <- ifelse(group=="CD44lo_CD27+", paste(group), "All_Others")
group2 <- ifelse(group=="T17", paste(group), "All_Others")
group2 <- ifelse(group=="T1", paste(group), "All_Others")

#group <- relevel(group, "CD44lo_CD27+")
#mod<-model.matrix(~group)
#fit = lmFit(exprs,mod)
#contrast.matrix<-c(0,1,1,1) # One vs All according to the old analysis

mod<-model.matrix(~group2)
fit = lmFit(exprs,mod)
colnames(fit)
contrast.matrix<-c(0,1) # One vs All
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)

top<-topTable(eb, adjust="BH",number=100000, p.value=1, lfc=0)

write.csv(top,file="IEL_vs_All_all.probes.csv")
write.csv(top,file="Naive.like_vs_All_all.probes.csv")
write.csv(top,file="T17_vs_All_all.probes.csv")
write.csv(top,file="T1_vs_All_all.probes.csv")

top<-topTable(eb, adjust="BH",number=100000, p.value=0.05, lfc=0)
head(top)

# Load annotation library
annodb <- "mogene20sttranscriptcluster.db"
ID<-rownames(top) # for the significant probes
#ID     <- featureNames(eset) # for all the dataset
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))
#my_frame <- data.frame(exprs(eset))
top$Symbol<-Symbol
top$Name<-Name
top$Entrez<-Entrez
head(top)

# Write out to a file:
write.csv(top,file="IEL_vs_All.csv")
write.csv(top,file="Naive.like_vs_All.csv")
write.csv(top,file="T17_vs_All.csv")
write.csv(top,file="Type1_vs_All.csv")

write.xlsx(top,file="IEL_vs_All.xlsx")
write.xlsx(top,file="Naive.like_vs_All.xlsx")
write.xlsx(top,file="T17_vs_All.xlsx")
write.xlsx(top,file="Type1_vs_All.xlsx")




# volcano IEL vs All ------------------------------------------------------
# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

rm(list=ls())

top <- read.csv(file="IEL_vs_All_all.probes.csv", row.names = 1)
head(top)
sel <- read.csv(file = "IEL_vs_All_volcano_selected.csv")
head(sel)
EnhancedVolcano(top,
                lab = top$Symbol,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 1e-03,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                legendPosition = "top",
                selectLab = sel$Symbol,
                drawConnectors = T,
                boxedLabels = T,
                title = "Volcano plot",
                subtitle = bquote(italic("IELs versus All"))
                )



# Overlaps ----------------------------------------------------------

rm(list=ls())

# Overlaps pairwise

list1 <- read.csv(file="Naive.like_vs_IEL.csv", row.names = 1)
list2 <- read.csv(file="T17_vs_IEL.csv", row.names = 1)
list3 <- read.csv(file="Type1_vs_IEL.csv", row.names = 1)
head(list3)

gene.lists <- list(Naive.like = unique(list1$Symbol[!list1$Symbol=="NA"]), 
                   T17 = unique(list2$Symbol[!list2$Symbol=="NA"]),
                   T1 = unique(list3$Symbol[!list3$Symbol=="NA"]))

group.colors <- c("#BABABA", "#66A61E","#1F78B4")
plot(euler(gene.lists, shape = "circle"), quantities = list(TRUE, cex=1.5, font = 3), 
     fills = group.colors, edges = group.colors,
     labels = list(labels=names(gene.lists), cex = 1.5),
     main = list(label = "\nOverlap between IEL pairwise comparisons", cex = 2.5, font = 1))
plot(euler(gene.lists, shape = "circle"), quantities = list(TRUE, cex=0, font = 3), 
     fills = group.colors, edges = group.colors,
     labels = list(labels=names(gene.lists), cex = 0),
     main = list(label = "\nOverlap between IEL pairwise comparisons", cex = 2.5, font = 1))


# Overlaps one vs all
IEL <- read.csv("IEL_vs_All.csv", row.names = 1)
head(IEL)
Naive.like <- read.csv("Naive.like_vs_All.csv", row.names = 1)
T17 <- read.csv("T17_vs_All.csv", row.names = 1)
T1 <- read.csv("Type1_vs_All.csv", row.names = 1)

gene.lists <- list(IEL = unique(IEL$Symbol[!IEL$Symbol=="NA"]), 
                   Naive.like = unique(Naive.like$Symbol[!Naive.like$Symbol=="NA"]),
                   T17 = unique(T17$Symbol[!T17$Symbol=="NA"]),
                   T1 = unique(T1$Symbol[!T1$Symbol=="NA"]))

group.colors <- c("darkorange", "#BABABA", "#66A61E","#1F78B4")

plot(euler(gene.lists, shape = "circle"), quantities = list(TRUE, cex=1.5, font = 3), 
     fills = group.colors, edges = group.colors,
     labels = list(labels=names(gene.lists), cex = 1),
     main = list(label = "\nOverlap between comparisons", cex = 1.5, font = 1))


int.all <- intersect(intersect(list1$Symbol, list2$Symbol), list3$Symbol)
#int.all <- Reduce(intersect, list(list1$Symbol, list2$Symbol, list3$Symbol))

setdiff(int.all, IEL$Symbol)
table(IEL$logFC > 0)



# Main heatmap ----------------------------------------------------------------

rm(list=ls())

IEL <- read.csv("IEL_vs_All.csv", row.names = 1)
head(IEL)

load("eset.RData")
exprs <- exprs(eset)

# for GEO submission (rma-normalized data)
# exprs[1:5,1:10]

pData <- pData(eset)

#top2 <- rbind(IEL, Naive.like, T17, T1)
top2 <- IEL

head(top2)
Sel<-eset[featureNames(eset)%in%rownames(top2), ]
pData<- pData(Sel)
exprs <- exprs(Sel)
colnames(exprs) <- paste0(pData$group, ".", letters[1:3])
head(exprs)
classvecCol <- brewer.pal(8, "Dark2")
group.colors <- c("#BABABA", "darkorange", "#1F78B4", "#66A61E")

par(mfrow=c(1,1), mar=c(5,5,5,5))
aheatmap(exprs, scale = "row", distfun = "euclidean",
        # cellwidth = 20, cellheight = 0.3, treeheight = 20,
         annCol = list(group = pData$group), annColors = list(group.colors),
         annLegend = T,
         main="DEGs\nall cell subtypes")



# Heatmaps with selected genelists ----------------------------------------------------------------

rm(list=ls())

load("eset.RData")
eset2 <- eset[, c(1:8,12,9:11)]
exprs <- exprs(eset2)
pData <- pData(eset2)

annodb <- "mogene20sttranscriptcluster.db" # Load annotation libraries
ID     <- featureNames(eset2) # for all the dataset
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))

rownames(exprs) <- Symbol
colnames(exprs) <- paste0(pData$group, ".", letters[1:3])
head(exprs)

genelist <- read.csv("genelists_heatmaps.csv", header = T)
head(genelist)

# Figure 3
sel.genes <- intersect(str_to_title(genelist$Fig3A.NK_Treg), rownames(exprs))
sel.exprs <- exprs[sel.genes, ]
head(exprs)
classvecCol <- brewer.pal(8, "Dark2")
group.colors <- c("#BABABA", "darkorange", "#1F78B4", "#66A61E")
#annRow2 <- c(rep("Transcription Factor", 41), rep("Chemokine, Chemokine Receptor", 13))
hm1 <- aheatmap(sel.exprs, scale = "row", distfun = "euclidean",
         Rowv = NA, Colv = NA,
         cellwidth = 12, cellheight = 13, treeheight = 20,
         annCol = list(group = pData$group), annColors = list(group.colors),
#         annRow = list(category = annRow2),
         annLegend = T, fontsize = 12, cexRow = 1, cexCol = 1,
         main="")

# Figure 4
head(genelist)
sel.genes <- intersect(str_to_title(genelist$Fig4A), rownames(exprs))
sel.exprs <- exprs[sel.genes, ]
head(sel.exprs)
hm1 <- aheatmap(sel.exprs, scale = "row", distfun = "euclidean",
                Rowv = NA, Colv = NA,
                cellwidth = 12, cellheight = 13, treeheight = 20,
                annCol = list(group = pData$group), annColors = list(group.colors),
                annLegend = T, fontsize = 12, cexRow = 1, cexCol = 1,
                main="")

# Figure 5
head(genelist)
sel.genes <- intersect(str_to_title(genelist$Fig5A), rownames(exprs))
sel.exprs <- exprs[sel.genes, ]
head(sel.exprs)
hm1 <- aheatmap(sel.exprs, scale = "row", distfun = "euclidean",
                Rowv = NA, Colv = NA,
                cellwidth = 12, cellheight = 13, treeheight = 20,
                annCol = list(group = pData$group), annColors = list(group.colors),
                annLegend = T, fontsize = 12, cexRow = 1, cexCol = 1,
                main="")

# Figure Supp1
head(genelist)
sel.genes <- intersect(str_to_title(genelist$Supp.1), rownames(exprs))
sel.exprs <- exprs[sel.genes, ]
head(sel.exprs)
hm1 <- aheatmap(sel.exprs, scale = "row", distfun = "euclidean",
                Rowv = NA, Colv = NA,
                cellwidth = 12, cellheight = 13, treeheight = 20,
                annCol = list(group = pData$group), annColors = list(group.colors),
                annLegend = T, fontsize = 12, cexRow = 1, cexCol = 1,
                main="")

# Figure Supp2
head(genelist)
sel.genes <- intersect(str_to_title(genelist$Supp.2), rownames(exprs))
sel.exprs <- exprs[sel.genes, ]
head(sel.exprs)
hm1 <- aheatmap(sel.exprs, scale = "row", distfun = "euclidean",
                Rowv = NA, Colv = NA,
                cellwidth = 12, cellheight = 13, treeheight = 20,
                annCol = list(group = pData$group), annColors = list(group.colors),
                annLegend = T, fontsize = 12, cexRow = 1, cexCol = 1,
                main="")




# GSEA --------------------------------------------------------------------

rm(list=ls())

IEL <- read.csv("IEL_vs_All_all.probes.csv", row.names = 1)
IEL$neg.logP <- -log10(IEL$P.Value)
IEL$gsea <- IEL$neg.logP * IEL$logFC
IEL <- IEL[order(IEL$gsea, decreasing = T), ]
head(IEL, 20)
tail(IEL, 20)

geneList = IEL[,"gsea"]
names(geneList) = as.character(IEL[,"Entrez"])
geneList = sort(geneList, decreasing = T)
head(geneList)
tail(geneList)


# GSEA (Biological Process)

gse_GO.BP <- gseGO(geneList = geneList,
                   OrgDb        = org.Mm.eg.db,
                   ont          = "BP", # using Biological Process terms
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
head(gse_GO.BP)
dotplot(gse_GO.BP, showCategory = 20) + ggtitle("GSEA dotplot for Gene Ontology : Biological Process")

## convert gene ID to Symbol
edox <- setReadable(gse_GO.BP, 'org.Mm.eg.db', 'ENTREZID')

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
p1+ggtitle("GSEA: Biological Process\nIEL vs All")

cnetplot(edox, categorySize="pvalue", foldChange=geneList, layout = "kk") +
  ggtitle("GSEA Gene-Concept Network for IEL vs All\nGene Ontology : Biological Process")

test <- as.data.frame(edox)
head(test)
write.xlsx(test, file = "GSEA_GO.BP_IEL.vs.All.xlsx")

for(i in 1:10){
  jpeg(filename = paste0("GSEA_GO_BP_IEL.vs.All_", gse_GO.BP$Description[i], ".jpeg"), width = 480, height = 480, quality = 100)
  p <- gseaplot2(gse_GO.BP, gse_GO.BP$ID[i], title = paste0("GSEA GO:BP in IEL vs. All\n", gse_GO.BP$Description[i]))
  plot(p)
  dev.off()
  rm(p)
}



# GSEA (KEGG)

gsea_KEGG <- gseKEGG(geneList     = geneList,
                     organism     = 'mmu',
                     minGSSize    = 100,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)
head(gsea_KEGG)
dotplot(gsea_KEGG, showCategory = 20) + ggtitle("GSEA dotplot for KEGG")

edox <- setReadable(gsea_KEGG, 'org.Mm.eg.db', 'ENTREZID')

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
p1+ggtitle("GSEA: KEGG\nIEL vs All")

cnetplot(edox, categorySize="pvalue", foldChange=geneList, layout = "kk") +
  ggtitle("GSEA Gene-Concept Network for IEL vs All\nKEGG")

test <- as.data.frame(edox)
head(test)
write.xlsx(test, file = "GSEA_KEGG_IEL.vs.All.xlsx")

for(i in 1:10){
  jpeg(filename = paste0("GSEA_KEGG_IEL.vs.All_", gsea_KEGG$Description[i], ".jpeg"), width = 480, height = 480, quality = 100)
  p <- gseaplot2(gsea_KEGG, gsea_KEGG$ID[i], title = paste0("GSEA KEGG in IEL vs. All\n", gsea_KEGG$Description[i]))
  plot(p)
  dev.off()
  rm(p)
}



# GSEA (Reactome)

gsea_Rx <- gsePathway(geneList, 
                      organism = "mouse",
                      minGSSize = 100,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE)
head(gsea_Rx)
dotplot(gsea_Rx, showCategory = 20) + ggtitle("GSEA dotplot for Reactome")

edox <- setReadable(gsea_Rx, 'org.Mm.eg.db', 'ENTREZID')

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
p1+ggtitle("GSEA: Reactome\nIEL vs All")

cnetplot(edox, categorySize="pvalue", foldChange=geneList, layout = "kk") +
  ggtitle("GSEA Gene-Concept Network for IEL vs All\nReactome")

test <- as.data.frame(edox)
head(test)
write.xlsx(test, file = "GSEA_Reactome_IEL.vs.All.xlsx")

for(i in 1:10){
  jpeg(filename = paste0("GSEA_Reactome_IEL.vs.All_", gsea_Rx$Description[i], ".jpeg"), width = 480, height = 480, quality = 100)
  p <- gseaplot2(gsea_Rx, gsea_Rx$ID[i], title = paste0("GSEA Reactome in IEL vs. All\n", gsea_Rx$Description[i]))
  plot(p)
  dev.off()
  rm(p)
}



# MSIGDB

immunesigdb <- msigdbr(species = "mouse", category = "C7", subcategory = "IMMUNESIGDB") %>% 
  dplyr::select(gs_name, entrez_gene)
head(immunesigdb)

em2 <- GSEA(geneList, TERM2GENE = immunesigdb)
head(em2)

gseaplot2(em2, "GSE14415_NATURAL_TREG_VS_TCONV_DN", title = paste0("GSEA Reactome in IEL vs. All\n", "GSE14415_NATURAL_TREG_VS_TCONV_DN"))
gseaplot2(em2, "GSE14415_NATURAL_TREG_VS_TCONV_UP", title = paste0("GSEA Reactome in IEL vs. All\n", "GSE14415_NATURAL_TREG_VS_TCONV_UP"))




# if using fgsea:

immunesigdb = msigdbr(species = "mouse", category = "C7", subcategory = "IMMUNESIGDB")
head(immunesigdb)

names(geneList) = as.character(IEL[,"Symbol"])

msigdbr_list = split(x = immunesigdb$gene_symbol, f = immunesigdb$gs_name)
fgsea.res <- fgsea(pathways = msigdbr_list, geneList)

head(fgsea.res)

fgsea.sel <- fgsea.res[fgsea.res$padj < 0.05, ]
head(fgsea.sel)

plotEnrichment(msigdbr_list[["GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN"]],
               geneList) + labs(title="GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN")

df <- apply(fgsea.sel,2,as.character)

write.csv(df, "MSIGDB_IEL.vs.All.csv")


plotEnrichment(msigdbr_list[["GSE14415_NATURAL_TREG_VS_TCONV_DN"]],
               geneList) + labs(title="GSE14415_NATURAL_TREG_VS_TCONV_DN")




# TF motif/target enrichment ----------------------------------------------------------------------

# prepared symbols genelist for IEL.vs.All, up and down separately, and saved as tab delimited for HOMER

rm(list=ls())

top <- read.csv(file="IEL_vs_All_all.probes.csv", row.names = 1)
head(top)

input_df <- top[, c(7,1,5)]
head(input_df)

# get gene lists
gsets_list <- get_gene_sets_list(source = "MSigDB",
                                 species = "Mus musculus",
                                 collection = "C3",
                                 subcollection = "TFT:GTRD") # alternatively use "TFT:TFT_Legacy"
length(gsets_list)
lapply(gsets_list, length)
head(gsets_list)[[1]][100]
head(gsets_list)[[2]][280]

output_df <- run_pathfindR(input_df, gene_sets = "Custom", custom_genes = gsets_list[[1]], custom_descriptions = gsets_list[[2]])
head(output_df)
dim(output_df)

write.csv(output_df, "pathfindR.output.csv")




# end ---------------------------------------------------------------------
sessionInfo()





