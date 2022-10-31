---
title: "Single-Cell RNA Sequencing Analysis"
author: "Huijian Feng"
date: "2020/8/5"
output:
  html_document:
    theme: cosmo
    highlight: pygments
    toc: yes
    toc_float: yes
---
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "D:/Tutorial/scRNA_seq_tutorial")

# rcran
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# install bioconductor library.
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
# bioconductor
options(useHTTPS=FALSE, BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

# A function: intalling the packages to R ->IPR
IPR <- function(packages, 
                sources = c("RCran", "Bioconductor", "Github")){
  installedpackages <- installed.packages()[ ,"Package"]
  existpackages <- intersect(packages, installedpackages)
  if(length(existpackages) != 0){
    print(paste0("Packages ", existpackages, " have been installed"))
    }
  newpackages <- setdiff(packages, installedpackages)
  if(length(newpackages) != 0){
    if(sources == "RCran"){
      install.packages(newpackages)
      }
    if(sources == "Bioconductor"){
      BiocManager::install(newpackages)
      }
    if(sources == "Github"){
      devtools::install_github(newpackages)
      }
    }
}

# R CRAN packages
# rpkg <- c("dplyr","ggplot2","matrixStats","pheatmap","DT","data.table", "Seurat")
# IPR(packages = rpkg, sources = "RCran")

# Bioconductor packages
# biopkg <- c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "genefilter", "ReactomePA")
# IPR(packages = biopkg, sources = "Bioconductor")
   

library(Seurat)
library(dplyr) 
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)
library(genefilter)
library(ReactomePA)
library(DT)
library(plotly)

## Example for set working directory
## setwd("D:/singlecell_pratice")

## loading the data 
adata <- Read10X("./counts/filtered_feature_bc_matrix")

## create the seurat object
adata <- CreateSeuratObject(counts = adata, 
                            project = "sc_tutorial",
                            min.cells = 3,
                            min.features = 200)
# parameter counts : input the "raws gene expression marix"
# parameter project : input the "name of our project"
# parameter min.cells : exclude some "genes" which are expressed in less than 3 cells
# parameter min.features : exclude some "cells" which are expressing less then 200 genes

## view the data
adata

## An object of class Seurat 
## 19037 features across 5100 samples within 1 assay 
## Active assay: RNA (19037 features)

## caculate the percentage of mitochondrial genes
adata[["percent.mt"]] <- PercentageFeatureSet(adata, 
                                             pattern = "^MT-")
# The [[ operator can add columns to object metadata. You can use "adata[[]]" to check object metadata.
# parameter pattern : match regular expression. "^MT-" means that the gene with prefix "MT-"

## function head() allow us to check the first 6 rows of data frame
head(adata[[]])
##                   orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACCCAAGACAGCTG sc_tutorial       8384         2550   9.148378
## AAACCCAAGTTAACGA sc_tutorial      11768         3017   6.296737
## AAACCCACAGTCGCAC sc_tutorial       6947         2033   5.570750
## AAACCCAGTAGCTTGT sc_tutorial       9354         2193   7.750695
## AAACCCATCTGCCCTA sc_tutorial       8970         2630   7.781494
## AAACGAAAGCAATTAG sc_tutorial       9291         2286  11.634916

## Check the distribution of the three QC metric
VlnPlot(adata, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

## parameter
nfeature_min = 200
nfeature_max = 5000
percent_mt = 20

## Check the distribution of the three QC metric
plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  geom_hline(yintercept = percent_mt, linetype="dotted", size = 1, colour = "#ae0404")

plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = nfeature_min, linetype="dotted", size = 1, colour = "#04ae59")+ 
  geom_hline(yintercept = nfeature_max, linetype="dotted", size = 1, colour = "#04ae59")
## plot 
plot_grid(plot1, plot2, nrow = 2)

## QC 
# 200< nFeature_RNA < 5000 & percent.mt < 20%
adata <- subset(adata, 
                subset = nFeature_RNA < nfeature_max & nFeature_RNA > nfeature_min & percent.mt < percent_mt)
## check data
adata
## An object of class Seurat 
## 19037 features across 4591 samples within 1 assay 
## Active assay: RNA (19037 features)

# normalizing data 
adata <- NormalizeData(adata, 
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)

##  caculates hvgs in single-cell data
hvgs = 3000
adata <- FindVariableFeatures(adata, 
                              selection.method = "vst", 
                              nfeatures = hvgs)
# "vst" is a method that model the mean-variance relationship in single-cell data.
# parameter nfeatures is to select top genes

## plot top 3000 HVGs
plot1 <- VariableFeaturePlot(adata)
plot1

# select top 3000 HVGs and perform the scale function
tgenes <- adata[["RNA"]]@var.features
# tgenes <- rownames(adata[["RNA"]])
adata <- ScaleData(adata, 
                   features = tgenes)

## perform the PCA by input top 3000 HVGs
# tgenes <- adata[["RNA"]]@var.features
npcs = 50
adata <- RunPCA(adata, 
                features = tgenes,
                npcs = npcs, 
                verbose = FALSE)
# parameter features to input user-defined genes
# parameter npcs, the number of PCA 

# adata <- JackStraw(adata, num.replicate = 100, dims = 30)
ElbowPlot(adata, 
          ndims = npcs)

# selcet the first 10 PCs
pcs = 20
# create nn graph 
adata <- FindNeighbors(adata, 
                       dims = 1:pcs, 
                       k.param = 20)
# louvain cluster
adata <- FindClusters(adata, 
                      resolution = 0.4)

## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 4591
## Number of edges: 162469
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9151
## Number of communities: 13
## Elapsed time: 0 seconds

# check the metadata
head(adata[[]])

##                   orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACCCAAGACAGCTG sc_tutorial       8384         2550   9.148378
## AAACCCAAGTTAACGA sc_tutorial      11768         3017   6.296737
## AAACCCACAGTCGCAC sc_tutorial       6947         2033   5.570750
## AAACCCAGTAGCTTGT sc_tutorial       9354         2193   7.750695
## AAACCCATCTGCCCTA sc_tutorial       8970         2630   7.781494
## AAACGAAAGCAATTAG sc_tutorial       9291         2286  11.634916
##                  RNA_snn_res.0.4 seurat_clusters
## AAACCCAAGACAGCTG               2               2
## AAACCCAAGTTAACGA               1               1
## AAACCCACAGTCGCAC               0               0
## AAACCCAGTAGCTTGT               8               8
## AAACCCATCTGCCCTA               6               6
## AAACGAAAGCAATTAG               6               6

# UMAP
adata <- RunUMAP(adata,
                 dims = 1:pcs)
# individual clusters 
p1 <- DimPlot(adata, 
              reduction = "umap", 
              label = TRUE)
ggplotly(p1)

# TSNE
adata <- RunTSNE(adata,
                 dims = 1:pcs)
# individual clusters 
p1 <- DimPlot(adata, 
              reduction = "tsne", 
              label = TRUE)
 
ggplotly(p1)

# find markers for every cluster compared to all remaining cells, report only the positive ones
adata_marker <- FindAllMarkers(adata,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# select top 10 markers for each cluster.
top10 <- adata_marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# heatmap the top 10 markers for each cluster.
DoHeatmap(adata, 
          features = top10$gene, 
          raster = T,
          label = F,
          size = 2) + 
  theme(axis.text.y = element_blank())

# B, NK 
FeaturePlot(adata, 
            features = c("MS4A1", # B
                         "GNLY", "NKG7"), # NK
            cols = c("grey", "red", "red"), 
            ncol = 2)

# CD14+Mono, Naive CD4+T
FeaturePlot(adata, 
            features = c("CD14", "LYZ", # CD14+Mono
                         "IL7R", "CCR7"), # Naive CD4+T
            cols = c("grey", "red", "red"), 
            ncol = 2)

# DC, FCGR3A+ Mono
FeaturePlot(adata, 
            features = c("FCER1A", "CST3", # DC 
                         "MS4A7", "FCGR3A"), # FCGR3A+ Mono
            cols = c("grey", "red", "red"), 
            ncol = 2)

# Memory CD4+, CD8+ T, Platelet
FeaturePlot(adata, 
            features = c("IL7R", "S100A4", # Memory CD4+
                         "CD8A", # CD8+T
                         "PPBP"), # Platelet
            cols = c("grey", "red", "red"), 
            ncol = 2)

adata[["cell_types"]] <- "undetermined"
adata@meta.data$cell_labels <-  as.character(adata@meta.data$seurat_clusters)
celltype_list = list("Naive CD4+ T"=c("0","8","9"),
                     "Memory CD4+"=c("3"),
                     "CD14+ Mono"=c("1"),
                     "B"=c("5","6"),
                     "CD8+ T"=c("2"),
                     "FCGR3A+ Mono"=c("7"),
                     "NK"=c("4"),
                     "DC"=c("10","12"),
                     "Platelet"=c("11"))
for(i in names(celltype_list)){
  adata@meta.data$cell_types[adata@meta.data$cell_labels %in% celltype_list[[i]]] <- i
}


p1 <- DimPlot(adata, 
              reduction = "umap", 
              label = T, 
              group.by = "cell_types")
ggplotly(p1)

# At first, we select the some cluster to do GO analysis
# 11:Platelet, 2:CD8+ T
marker <- adata_marker[as.character(adata_marker$cluster) %in% c("2","11"), ]
marker$cluster <- as.character(marker$cluster)
marker <-  split(marker$gene, marker$cluster)

# we can select the top 100 genes.
# t100 <- adata_marker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
# t100 <- split(top20$gene, top20$cluster)
entrezid <-  lapply(marker, function(gr) as.numeric(bitr(gr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID))

# PBMC GO
pvalueCutoff = 0.01
qvalueCutoff = 0.01
pbmc_go <- compareCluster(entrezid,
                          OrgDb='org.Hs.eg.db',
                          fun='enrichGO',
                          pvalueCutoff = pvalueCutoff,
                          qvalueCutoff = qvalueCutoff,
                          ont = "BP",
                          readable=T)
# plot the go 
dotplot(pbmc_go,
        title = paste0("PBMC Gene Ontology (qval < ", qvalueCutoff, ")"))

# view the go table data
go_table = pbmc_go@compareClusterResult
datatable(go_table, filter="top", options=list(pageLength = 10))
# we can save the GO table as the .html file
# datatable(go_table, filter="top", options=list(pageLength = 10)) %>% saveWidget("pbmc_go.html")

# And then we can do the pathway analysis
pbmc_pa<- compareCluster(entrezid,
                         organism = "human",
                         fun='enrichPathway',
                         pvalueCutoff = pvalueCutoff,
                         qvalueCutoff = qvalueCutoff,
                         readable=T)
# plot the PA 
dotplot(pbmc_pa,
        title = paste0("PBMC Pathway (qval < ", qvalueCutoff, ")"))

# pa_table = pbmc_pa@compareClusterResult
