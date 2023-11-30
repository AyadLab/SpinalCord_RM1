#packages to load
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(SPOTlight)
#library(SingleCellExperiment)
#library(SpatialExperiment)
#library(scater)
library(scran)
#devtools::install_github("meichendong/SCDC") #to install SCDC
suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

#cleaning of spatial data

#obj <- readRDS("/data/ccbb/ayad_lab/spatial_RM1/sp.rds")
#This object contains only 4 images that need's to be deconvoluted and downstream analysis is done..
#obj <- readRDS("/data/ccbb/ayad_lab/spatial_RM1/subset.RDS")
#images_to_remove <- c("Young_3dpi_rep2","Aged_3dpi_rep1", "Aged_3dpi_rep2", "Young_3dpi_rep1.1","Aged_3dpi_rep1.1","Aged_3dpi_rep2.1","Young_3dpi_rep1.2","Young_3dpi_rep2.2","Aged_3dpi_rep2.2","Young_3dpi_rep1.3","Young_3dpi_rep2.3","Aged_3dpi_rep1.3")
#image_metadata <- obj@images
#image_metadata <- image_metadata[!names(image_metadata) %in% images_to_remove]
#obj@images <- image_metadata
#obj <- ScaleData(obj, verbose = FALSE)
#obj <- FindVariableFeatures(obj)
#obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
#obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
sc_rna <- FindNeighbors(sc_rna, reduction = "pca", dims = 1:30)
sc_rna <- FindClusters(sc_rna , resolution = 0.5)
#saveRDS(obj,"subset.RDS")
#Young_3dpi_rep1
#Young_3dpi_rep2.1
#Aged_3dpi_rep1.2
#Aged_3dpi_rep2.3

obj <- readRDS("/data/ayad_lab/spatial_RM1/subset.RDS")
#sc_rna <- readRDS("/data/ccbb/ayad_lab/spatial_RM1/ScRNA_Mouse_T8_UI-3dpi_processed_CF.rds")
sc_rna <- readRDS("/data/ayad_lab/tabulaParalytica/tabula_rnaseqdata.RDS")
image_keys <- unique(obj@meta.data$SampleName)

sub_obj <- obj[, obj$SampleName %in% image_keys[[2]]]
image_remove <- c("Aged_3dpi_rep1.2","Young_3dpi_rep1","Aged_3dpi_rep2.3")
image_metadata <- obj@images
image_metadata <- image_metadata[!names(image_metadata) %in% image_remove]
sub_obj@images <- image_metadata
sc_rna@meta.data
sc_rna@meta.data <- sc_rna@meta.data[, colSums(is.na(sc_rna@meta.data)) == 0]

ScaleData(sc_rna, verbose = FALSE)
RunUMAP(sc_rna, reduction = "pca", dims = 1:30)
#FindNeighbors(sc_rna, reduction = "pca", dims = 1:30)
rownames(sc_rna@meta.data) <- colnames(sc_rna)
p1 <- DimPlot(obj, reduction = "umap")
p2 <- DimPlot(sc_rna, reduction = "umap",group.by = "seurat_clusters",raster = FALSE)

pdf("tabula_results/umap_tabula_cluster.pdf",width = 50,height = 20)
p2
dev.off()

# reference <- SCTransform(sc_rna, ncells = 3000, verbose = FALSE, method = "poisson") %>%
#   RunPCA(verbose = FALSE) %>%
#   RunUMAP(dims = 1:30)
# reference <- SCTransform(sc_rna)
unique(sc_rna@meta.data$layer3)
markers_sc <- FindAllMarkers(sc_rna, only.pos = TRUE, logfc.threshold = 0.1,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05,assay = "RNA")

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(sub_obj), ]

# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(150, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(30, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))
#eset_SC :create a reference expression set
eset_SC <- ExpressionSet(assayData = as.matrix(sc_rna@assays$RNA@counts[m_feats,
]), phenoData = AnnotatedDataFrame(sc_rna@meta.data))
#eset_ST : create a spatial expression set
eset_ST <- ExpressionSet(assayData = as.matrix(sub_obj@assays$Spatial@counts[m_feats,
]), phenoData = AnnotatedDataFrame(sub_obj@meta.data))
#Calculate proportions of cell type using SCDC_prop function
#deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "SingleR2",
#                                     ct.sub = as.character(unique(eset_SC$SingleR2)))

deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "layer3",
                                     ct.sub = as.character(unique(eset_SC$layer3)))
head(deconvolution_crc$prop.est.mvw)
#annot <- sub_obj$layer
#coords <- GetTissueCoordinates(sub_obj, which = "centroids")
#colnames(coords) <- c("x", "y")
#coords[is.na(colnames(pos))] <- NULL
#proportions table
deconprop<- deconvolution_crc$prop.est.mvw

image_4 <- plotInteractions(deconprop, which = "heatmap", metric = "prop")

ct <- colnames(deconprop)
deconprop[deconprop < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

scatterpie_image_1<- plotSpatialScatterpie(
  x = subset_brain_7,
  y = deconprop,
  cell_types = colnames(deconprop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.3) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))
pdf("RM1_spinalcord_outs/spinaldecov_image1_test.pdf")
pdf("spinaldecov_image1_test_with_tabuladata.pdf")

scatterpie_image_1
dev.off()

pdf("tabula_results/plotinteractions_Young_3dpi_rep1.pdf")
image_4
dev.off()
