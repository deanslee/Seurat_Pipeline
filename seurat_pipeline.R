# This is vastly simplified code to illustrate the work flow I execute for single-cell RNA-seq data generated from my lab.

# There are four primary goals.

# The first goal is to combine data from both experimental and control samples so that at the end of this analysis pipeline, cells from both samples are subjected to the same cluster assignment. This makes it possible to test for cluster-specific effects of the experimental condition. This combination step is done through canonical correlation analysis.

# The second goal is to define biologically meaningful clusters of cells from the scRNA-seq data. These clusters should correspond to actual cell types in the tissue of interest. Ex. neurons, glia, oligodendrocytes, subtypes of neurons, subtypes of glia, subtypes of oligodendrocytes.

# The third goal is to use these cluster assignments to perform differential gene expression analysis. Each cluster/cell type now has a list of marker genes that characterize it.

# The fourth goal is to test for the effect of the experimental condition, compared to the control condition, on gene expression in each cluster. This reveals to us what the experimental condition is doing in each cell type that we have identified.

# install libraries

library(Seurat)
library(dplyr)
library(Matrix)
library(MAST)
library(plotly)
library(reticulate)
library(scales)

# read 10x Genomics data: experimental sample and control sample

EXP.data = Read10X(data.dir = "/Users/deanlee/Dropbox/10X_EXP/cellranger/counts/mm10-1.2.0_premrna_EXP/")
CNTRL.data = Read10X(data.dir = "/Users/deanlee/Dropbox/10X_EXP/cellranger/counts/mm10-1.2.0_premrna_CNTRL/")

# create Seurat object

EXP = CreateSeuratObject(raw.data = EXP.data)
CNTRL = CreateSeuratObject(raw.data = CNTRL.data)

# set up EXP and CNTRL object, filter out cells contaminated by mitochondrial, ribosomal, and apoptotic transcripts, filter out cells with abnormally high UMI

EXP = NormalizeData(EXP)
EXP = SetAllIdent(EXP, id = 'orig.ident')

mito.genes = grep("^mt-", rownames(EXP@data), value = T)
EXP.percent.mito = colSums(as.array(expm1(EXP@data[mito.genes,])))/colSums(as.array(expm1(EXP@data)))
EXP = AddMetaData(EXP, EXP.percent.mito, "percent.mito")
ribo.genes = grep("^Rps|Rpl", rownames(EXP@data), value = T)
EXP.percent.ribo = colSums(as.array(expm1(EXP@data[ribo.genes,])))/colSums(as.array(expm1(EXP@data)))
EXP = AddMetaData(EXP, EXP.percent.ribo, "percent.ribo")
apop.genes = c("Apaf1","Cycs","Bax","Tnf","Hspa2","Junb")
EXP.percent.apop = colSums(as.array(expm1(EXP@data[apop.genes,])))/colSums(as.array(expm1(EXP@data)))
EXP = AddMetaData(EXP, EXP.percent.apop, "percent.apop")

EXP = FilterCells(EXP, subset.names = c("nUMI", "percent.mito", "percent.ribo", "percent.apop"), high.thresholds = c(3*median(EXP@meta.data$nUMI), 0.05, 0.01, 0.0025))

CNTRL = NormalizeData(CNTRL)
CNTRL = SetAllIdent(CNTRL, id = 'orig.ident')

CNTRL.percent.mito = colSums(as.array(expm1(CNTRL@data[mito.genes,])))/colSums(as.array(expm1(CNTRL@data)))
CNTRL = AddMetaData(CNTRL, CNTRL.percent.mito, "percent.mito")
CNTRL.percent.ribo = colSums(as.array(expm1(CNTRL@data[ribo.genes,])))/colSums(as.array(expm1(CNTRL@data)))
CNTRL = AddMetaData(CNTRL, CNTRL.percent.ribo, "percent.ribo")
CNTRL.percent.apop = colSums(as.array(expm1(CNTRL@data[apop.genes,])))/colSums(as.array(expm1(CNTRL@data)))
CNTRL = AddMetaData(CNTRL, CNTRL.percent.apop, "percent.apop")

CNTRL = FilterCells(CNTRL, subset.names = c("nUMI", "percent.mito", "percent.ribo", "percent.apop"), high.thresholds = c(3*median(CNTRL@meta.data$nUMI), 0.04, 0.009, 0.0025))

# find most highly variable genes for canonical correlation analysis
EXP = FindVariableGenes(EXP, y.cutoff = 1, x.low.cutoff = 0.0125, x.high.cutoff = 3)
CNTRL = FindVariableGenes(CNTRL, y.cutoff = 1, x.low.cutoff = 0.0125, x.high.cutoff = 3)

EXP = ScaleData(EXP, genes.use = EXP@var.genes, vars.to.regress = c("nUMI","percent.mito"), model.use = 'negbinom', block.size = 5000, display.progress = F)
CNTRL = ScaleData(CNTRL, genes.use = CNTRL@var.genes, vars.to.regress = c("nUMI","percent.mito"), model.use = 'negbinom', block.size = 5000, display.progress = F)

EXP@meta.data[,"condition"] = "EXP"
CNTRL@meta.data[,"condition"] = "CNTRL"

EXPvg = head(rownames(EXP@hvg.info), 1500)
CNTRLvg = head(rownames(CNTRL@hvg.info), 1500)
genes.use = unique(c(EXPvg, CNTRLvg))
genes.use = intersect(genes.use, rownames(EXP@scale.data))
genes.use = intersect(genes.use, rownames(CNTRL@scale.data))

CCA = RunCCA(EXP, CNTRL, genes.use = genes.use, num.cc = 50, add.cell.id1 = "EXP", add.cell.id2 = "CNTRL")
CCA = AlignSubspace(CCA, reduction.type = "cca", grouping.var = "condition", dims.align = 1:50)

# cluster analysis

dims = 1:30 
CCA = FindClusters(CCA, reduction.type = "cca.aligned", dims.use = dims, k.param = 15, save.SNN = T, n.start = 20, algorithm = 3, resolution = 2, print.output = F, force.recalc = T)
CCA = ReorderIdent(CCA, feature = "ACC1", reorder.numeric = T)
CCA <- RunUMAP(CCA, reduction.use = "cca.aligned", dims.use = dims, max.dim = 3)

# visualize clusters in 2D space

set.seed(1)
cols = sample(scales::hue_pal()(length(levels(CCA@ident))), length(levels(CCA@ident)))
DimPlot(CCA, reduction.use = "umap", do.label = T, cols.use = cols, pt.size = 0.1)

# visualize clusters in 3D space

umap1 = CCA@dr$umap@cell.embeddings[,1]
umap2 = CCA@dr$umap@cell.embeddings[,2]
umap3 = CCA@dr$umap@cell.embeddings[,3]

set.seed(1)
cols = sample(scales::hue_pal()(length(levels(CCA@ident))), length(levels(CCA@ident)))
plot_ly(x = ~umap1, y = ~umap2, z = ~umap3, color = ~CCA@ident, colors = cols, marker = list(size = 2)) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'UMAP1'), yaxis = list(title = 'UMAP2'), zaxis = list(title = 'UMAP3')))

# differential gene expression analysis to identify marker genes for each cluster/cell type

CCA.markers = FindAllMarkers(CCA, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T, min.diff.pct = .1, min.pct = .1, latent.vars = c('nUMI','percent.mito','percent.ribo'), do.print = T)
write.table(CCA.markers, "/Users/deanlee/Dropbox/10X_EXP/CCA/CCA_markers.tsv", quote = F, sep = "\t", row.names = F)

# visualize Cell Type 1 in 3D plot, where two marker genes overlap

gene1 = "Siglech"
gene2 = "Runx1"
seuObj = CCA

neither = seuObj@scale.data[gene1,]<0 & seuObj@scale.data[gene2,]<0
neither = as.numeric(neither)
neither[neither==1] = 0
gene1only = seuObj@scale.data[gene1,]>0 & seuObj@scale.data[gene2,]<0
gene1only = as.numeric(gene1only)
gene1only[gene1only==1] = 1
gene2only = seuObj@scale.data[gene1,]<0 & seuObj@scale.data[gene2,]>0
gene2only = as.numeric(gene2only)
gene2only[gene2only==1] = 2
both = seuObj@scale.data[gene1,]>0 & seuObj@scale.data[gene2,]>0
both = as.numeric(both)
both[both==1] = 3
two_genes = factor(neither + gene1only + gene2only + both)

plot_ly(x = ~umap1, y = ~umap2, z = ~umap3, color = ~two_genes, colors = c("grey","red","blue","purple"), marker = list(size = 2)) %>% add_markers() %>% layout(title = paste(gene1,",",gene2), scene = list(xaxis = list(title = 'UMAP1'), yaxis = list(title = 'UMAP2'), zaxis = list(title = 'UMAP3')))