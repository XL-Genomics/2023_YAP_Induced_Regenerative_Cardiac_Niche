####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- '91_NCR_Revision_1'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2022_yap5sa_st_rli/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/stRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject(Machine = 'Rivendell', Ver = Ver, Part = Step, Catagory = 'mouse',
                Project_dir = '2022_yap5sa_st_rli', Data_drive = 'ithil')
library('ggbeeswarm')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load objects and update metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drop.srt <- readRDS('external/fransico_data/reference_scrnaseq/Reference_combined_compendium_&_yap5sa_CMs_v2.rds')
drop.srt$Cell_type <- revalue(drop.srt$celltype_bulk, c('MAC' = 'MP', 'EC_Flt1+' = 'EC2',
                                                        'EC_Pecam1+' = 'EC1', 'SMC' = 'Mural'))
drop.srt$Cell_type <- factor(drop.srt$Cell_type, levels = c('CM', 'CF', 'MP', 'EC1', 'EC2', 'Mural'))
drop.srt$Cell_type2 <- revalue(drop.srt$celltype, c('MAC' = 'MP', 'EC_Flt1+' = 'EC2',
                                                    'EC_Pecam1+' = 'EC1', 'SMC' = 'Mural',
                                                    'CM1' = 'aCM1', 'CM2' = 'aCM2'))
drop.srt$Cell_type2 <- factor(drop.srt$Cell_type2, levels = c('aCM1', 'aCM2', 'CF', 'MP', 'EC1', 'EC2', 'Mural'))
drop.srt$Sample <- revalue(drop.srt$orig.ident, c('YAP5SA.1' = 'YAP5SA 3', ## yap5sa_wh_s1
                                                  'YAP5SA.2' = 'YAP5SA 2', ## yap5sa_wh_s2
                                                  'YAP5SA.cm' = 'YAP5SA 1', ## yap5sa_cm_s1
                                                  'CTRL.1.cm' = 'Control 1', ## ctrl_cm_s1
                                                  'CTRL.3' = 'Control 2', ## ctrl_wh_s1
                                                  'CTRL.2' = 'Control 3')) ## ctrl_wh_s2
drop.srt$Sample <- factor(drop.srt$Sample, levels = c(paste('Control', 1:3), paste('YAP5SA', 1:3)))

drop_cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
new_cm_dim <- read.csv(paste0('~/Documents/Bioinformatics/project/2022_yap5sa_st_rli/',
                              'meta/mouse_v0/90_Plot/09.cm_umap_coordinates.csv'))
colnames(new_cm_dim) <- paste0('NEWUMAP_', 1:2)
rownames(new_cm_dim) <- Cells(drop_cm.srt)
drop_cm.srt@reductions$newumap <- CreateDimReducObject(embeddings = as.matrix(new_cm_dim),
                                                       assay = 'SCT', key = 'NEWUMAP_')
drop_cm.srt$Sample <-  drop.srt$Sample[Cells(drop_cm.srt)]

st.srt <- readRDS('individual/part01.YAP5SA_st_ventricle.srt.rds')
st_wh.srt <- readRDS('individual/part01.YAP5SA_st.srt.rds')

sc.srt <- readRDS('/Volumes/ithil/project/2022_nsp13_fmeng/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')

Color_cell_type <- mycol_10[3:8]
Color_Zone <- mycol_20[c(4, 9, 13)]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #1. Reviewer #1 Q1 (using harmonized PCs -- Not used) ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st_wh.srt <- RunPCA(st_wh.srt, assay = 'SCT')
st_wh.srt <- RunHarmony(st_wh.srt, assay.use = 'SCT', group.by.vars = 'Genotype')
st_wh.srt <- RunUMAP(st_wh.srt, dims = 1:50, reduction = 'harmony', min.dist = 0.2)
st_wh.srt <- FindNeighbors(st_wh.srt, reduction = 'harmony', dims = 1:50)
st_wh.srt <- FindClusters(st_wh.srt, resolution = 0.3)

DimPlot2(st_wh.srt, group.by = 'SCT_snn_res.0.3', label = T, repel = T, cols = mycol_10)

st_wh.srt$Niche <- factor(as.numeric(st_wh.srt$SCT_snn_res.0.3), levels = 1:8)

Idents(st_wh.srt) <- 'Niche'
SpatialDimPlot(st_wh.srt, image.alpha = 0) & scale_fill_manual(values = mycol_10)

mtx <- st_wh.srt@meta.data[, grepl('^Decon|^Niche', colnames(st_wh.srt@meta.data))]
dim(mtx)
mtx <- mtx[! is.na(mtx$Decon_CF),]

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
for(i in 1:7){mtx2[, i] <- scale(mtx2[, i])}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p1 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Decon %\nRow Z-score', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p2 <- DimPlot2(st_wh.srt, group.by = 'Niche', label = T, repel = F, cols = mycol_10)/
        DimPlot2(st_wh.srt, group.by = 'Genotype')
p3 <- SpatialDimPlot(st_wh.srt, group.by = 'Niche', image.alpha = 0, ncol = 1, stroke = 0) &
        scale_fill_manual(values = mycol_10)
p3[[1]] <- p3[[1]] +NoLegend()
p3 | p2 | p1
PlotPDF('01.1.st_niches', 16, 10)
p3 | p2 | p1
dev.off()

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
for(i in 1:5){mtx2[i , ] <- scale(mtx2[i, ])}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p4 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Decon %\nColumn Z-score', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p4

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p5 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Decon %', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p5

mtx_c <- mtx[grepl('^Control_', rownames(mtx)), ]
mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx_c[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx_c[mtx_c$Niche == c(1:5)[i], j], )
        }
}
mtx3 <- matrix(NA, 5, 7)
mtx_y <- mtx[grepl('^YAP5SA_', rownames(mtx)), ]
colnames(mtx3) <- str_split(colnames(mtx_y[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx3) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx3[i, j] <- mean(mtx_y[mtx_y$Niche == c(1:5)[i], j], )
        }
}
data <- reshape2::melt(mtx2)
data$G <- 'Ctrl'
data2 <- reshape2::melt(mtx3)
data2$G <- 'Y5sa'
data <- rbind(data, data2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p6 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Decon %', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis() +
        facet_wrap(~G)
p6

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
mtx2 <- log10(mtx2)
mtx2 <- cbind(mtx2, mtx2[,'CM2']/mtx2[,'CM1'])
colnames(mtx2)[8] <- 'CM2_CM1'
data <- reshape2::melt(mtx2[, 1:7])
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p7 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Log10 Decon %', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p7

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
mtx2 <- log10(mtx2)
mtx2 <- cbind(mtx2, mtx2[,'CM2']/mtx2[,'CM1'])
colnames(mtx2)[8] <- 'CM2_CM1'
data <- reshape2::melt(mtx2[, c(1, 4:8)])
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CF', 'EndoC', 'CM2_CM1'))
p8 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        theme_classic() +
        theme(aspect.ratio = 7/5) +
        labs(fill = 'Log10 Decon %', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p8

mtx2 <- matrix(NA, 5, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:5)
for(i in 1:5){
        for(j in 1:7){
                mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:5)[i], j], )
        }
}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p9.1 <- ggplot(data[data$Var2 %in% c('CM1', 'CM2'), ]) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        theme_classic() &
        theme(aspect.ratio = 2/5)
p9.2 <- ggplot(data[!data$Var2 %in% c('CM1', 'CM2'), ]) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        theme_classic() &
        theme(aspect.ratio = 5/5)
p9 <- p9.1/p9.2 &
        scale_fill_distiller(palette = 'RdBu') &
        labs(fill = 'Decon %', y = 'Cell Type', x = 'Niche') &
        scale_y_discrete(limits = rev) &
        RotatedAxis()
p9

meta <- st_wh.srt@meta.data[, grepl('^Decon|Niche|Genotype|Zone', colnames(st_wh.srt@meta.data))]
meta$Cell_ID <- colnames(st_wh.srt)

Idents(st_wh.srt) <- 'Niche'
tmp.srt <- PrepSCTFindMarkers(st_wh.srt)
mk <- FindAllMarkers(tmp.srt, only.pos = T)
mk <- mk[mk$p_val_adj < 0.05, ]
colnames(mk)[6] <- 'niche'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WriteCSV(meta, 'part91.niche_decon_meta_table')
WriteCSV(mk, 'part91.niche_signature_genes')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #1. Reviewer #1 Q1 (using deconvolved compostion -- Didn't work)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('compositions')

mat <- st_wh.srt@meta.data[, grep('^Decon', colnames(st_wh.srt@meta.data), value = T)]
mat <- mat[!is.na(mat[, 1]), ]

baseILR <- ilrBase(x = mat, method = "basic")
cell_ilr <- as.matrix(ilr(mat, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

st_v.srt <- st_wh.srt[, !is.na(st_wh.srt$Decon_CF)]
st_v.srt@assays$ILR <- CreateAssayObject(data = t(cell_ilr), min.cells = 0, min.features = 0, check.matrix = T)
st_v.srt@assays$ILR@key <- "ilr_"
DefaultAssay(st_v.srt) <- 'ILR'
st_v.srt <- ScaleData(st_v.srt, features = rownames(st_v.srt))
st_v.srt <- RunPCA(st_v.srt, npcs = 6, features = rownames(st_v.srt))
st_v.srt <- RunUMAP(st_v.srt, dims = 1:5, reduction = 'pca')
st_v.srt <- FindNeighbors(st_v.srt, reduction = 'pca', dims = 1:5) |> FindClusters(res = 0.2)

DimPlot2(st_v.srt)
SpatialDimPlot(st_v.srt, image.alpha = 0) &
        scale_fill_manual(values = mycol_14)

comp_umap <- uwot::umap(cell_ilr, n_neighbors = 300, n_epochs = 1000) %>%
        as.data.frame() %>%
        mutate(row_id = rownames(cell_ilr))

meta <- st_wh.srt@meta.data[rownames(cell_ilr), ]
meta$row_id <- rownames(meta)

comp_umap %>%
        left_join(meta, by = c("row_id")) %>%
        ggplot(aes(x = V1, y = V2, color = Niche)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = mycol_20) +
        theme_classic() +
        xlab("UMAP1") +
        ylab("UMAP2")


## For Hassan:
st_v.srt@reductions$spatial@assay.used <- 'ILR'
st_v2.srt <- DietSeurat(st_v.srt, assays = 'ILR', dimreducs = 'spatial')
saveRDS(st_v2.srt, 'tmp/vent_st_with_ilr_matrix.seurat.rds')
st_v.srt <- st_wh.srt[, !is.na(st_wh.srt$Decon_CF)]
st_v.srt@assays$Decon <- CreateAssayObject(data = t(mat), min.cells = 0, min.features = 0, check.matrix = T)
st_v.srt@assays$Decon@key <- "decon_"
DefaultAssay(st_v.srt) <- 'Decon'
st_v.srt@reductions$spatial@assay.used <- 'Decon'
st_v2.srt <- DietSeurat(st_v.srt, assays = 'Decon', dimreducs = 'spatial')
saveRDS(st_v2.srt, 'tmp/vent_st_with_decon_pct_matrix.seurat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #1. Reviewer #1 Q1 (using PCs but split Ctrl and Yap5sa sections)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st_wh.srt <- AddModuleScore2(st_wh.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT',
                             names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)), return_z = T)

st_wh_ctrl.srt <- RunPCA(st_wh.srt[, st_wh.srt$Sample == 'Control'], assay = 'SCT')
st_wh_ctrl.srt <- RunUMAP(st_wh_ctrl.srt, dims = 1:50, reduction = 'pca', min.dist = 0.2)
st_wh_ctrl.srt <- FindNeighbors(st_wh_ctrl.srt, reduction = 'pca', dims = 1:50)
st_wh_ctrl.srt <- FindClusters(st_wh_ctrl.srt, resolution = seq(0.1, 1, 0.1))
st_wh_ctrl.srt@images$YAP5SA <- NULL
p1 <- DimPlot2(st_wh_ctrl.srt, group.by = paste0('SCT_snn_res.', seq(0.1, 1, 0.1)),
               label = T, repel = T, cols = mycol_20, ncol = 4)
p2 <- SpatialDimPlot(st_wh_ctrl.srt, group.by = paste0('SCT_snn_res.', seq(0.1, 1, 0.1)),
                     image.alpha = 0, stroke = 0, ncol = 4) &
        scale_fill_manual(values = mycol_20)
PlotPDF('01.3.st_wh_ctrl_niche', 20, 15)
p1
p2
dev.off()

st_wh_y5sa.srt <- RunPCA(st_wh.srt[, st_wh.srt$Sample != 'Control'], assay = 'SCT')
st_wh_y5sa.srt <- RunUMAP(st_wh_y5sa.srt, dims = 1:50, reduction = 'pca', min.dist = 0.2)
st_wh_y5sa.srt <- FindNeighbors(st_wh_y5sa.srt, reduction = 'pca', dims = 1:50)
st_wh_y5sa.srt <- FindClusters(st_wh_y5sa.srt, resolution = seq(0.1, 1, 0.1))
st_wh_y5sa.srt@images$Control <- NULL
p1 <- DimPlot2(st_wh_y5sa.srt, group.by = paste0('SCT_snn_res.', seq(0.1, 1, 0.1)),
               label = T, repel = T, cols = mycol_20, ncol = 4)
p2 <- SpatialDimPlot(st_wh_y5sa.srt, group.by = paste0('SCT_snn_res.', seq(0.1, 1, 0.1)),
                     image.alpha = 0, stroke = 0, ncol = 4) &
        scale_fill_manual(values = mycol_20)
PlotPDF('01.4.st_wh_y5sa_niche', 20, 15)
p1
p2
dev.off()

st_wh_ctrl.srt$Niche_ctrl <- factor(as.numeric(st_wh_ctrl.srt$SCT_snn_res.0.5), levels = 1:7)
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('5' = '0'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('7' = '5'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('4' = '7'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('6' = '4'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('0' = '6'))
st_wh_ctrl.srt$Niche_ctrl <- factor(st_wh_ctrl.srt$Niche_ctrl, levels = 1:7)
DimPlot2(st_wh_ctrl.srt, group.by = 'Niche_ctrl', label = T, repel = T, cols = mycol_10)
p1 <- SpatialDimPlot(st_wh_ctrl.srt, group.by = 'Niche_ctrl', image.alpha = 0, stroke = 0) +
        scale_fill_manual(values = mycol_10)

st_wh_y5sa.srt$Niche_y5sa <- factor(as.numeric(st_wh_y5sa.srt$SCT_snn_res.0.4), levels = 1:7)
DimPlot2(st_wh_y5sa.srt, group.by = 'Niche_y5sa', label = T, repel = T, cols = mycol_10)
p2 <- SpatialDimPlot(st_wh_y5sa.srt, group.by = 'Niche_y5sa', image.alpha = 0, stroke = 0) +
        scale_fill_manual(values = mycol_10)
p1 | p2
PlotPDF('01.5.st_wh_split_niche', 10, 5)
p1 | p2
dev.off()

mtx <- st_wh_ctrl.srt@meta.data[, grepl('^Decon|^Niche_ctrl', colnames(st_wh_ctrl.srt@meta.data))]
dim(mtx)
mtx <- mtx[! is.na(mtx$Decon_CF),]
Table(mtx$Niche_ctrl)
mtx2 <- matrix(NA, 7, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:7)
for(i in 1:7){ for(j in 1:7){ mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:7)[i], j], ) } }
#for(i in 1:7){mtx2[, i] <- scale(mtx2[, i])}
data <- reshape2::melt(mtx2[c(1:4, 6), ])
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p1 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = log2(value*100+1)))
p1

mtx <- st_wh_y5sa.srt@meta.data[, grepl('^Decon|^Niche_y5sa', colnames(st_wh_y5sa.srt@meta.data))]
dim(mtx)
mtx <- mtx[! is.na(mtx$Decon_CF),]
Table(mtx$Niche_y5sa)
mtx2 <- matrix(NA, 7, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:7)
for(i in 1:7){ for(j in 1:7){ mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:7)[i], j], ) } }
#for(i in 1:7){mtx2[, i] <- scale(mtx2[, i])}
data <- reshape2::melt(mtx2[c(1:4, 6), ])
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
p2 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = log2(value*100+1)))
p <- p1 + p2 &
        scale_fill_viridis_c(limits = c(0, 6.2)) &
        theme_classic() &
        theme(aspect.ratio = 7/5) &
        labs(fill = 'Log Abundance', y = 'Cell Type', x = 'Niche') &
        scale_y_discrete(limits = rev) &
        RotatedAxis()
p
PlotPDF('01.6.st_wh_split_niche_cell_type_abundance', 10, 5)
p
dev.off()

## Check Triad
mtx <- st_wh_ctrl.srt@meta.data[, grepl('^Zscore|^Niche_ctrl', colnames(st_wh_ctrl.srt@meta.data))]
dim(mtx)
Table(mtx$Niche_ctrl)
mtx2 <- matrix(NA, 7, 4)
colnames(mtx2) <- str_split(colnames(mtx[1:4]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:7)
for(i in 1:7){ for(j in 1:4){ mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:7)[i], j], ) } }
data <- reshape2::melt(mtx2[c(1:4, 6), ])
data$Var2 <- factor(data$Var2, levels = c('CM1', 'CM2', 'C3ar1', 'C3'))
data$Var2 <- revalue(data$Var2, replace = c('CM1'='CM1_marker', 'CM2'='CM2_marker',
                                            'C3ar1'='C3ar1+MP_marker', 'C3'='C3+CF_marker'))
p3 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value))
p3
mtx <- st_wh_y5sa.srt@meta.data[, grepl('^Zscore|^Niche_y5sa', colnames(st_wh_y5sa.srt@meta.data))]
dim(mtx)
Table(mtx$Niche_y5sa)
mtx2 <- matrix(NA, 7, 4)
colnames(mtx2) <- str_split(colnames(mtx[1:4]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- paste('Niche', 1:7)
for(i in 1:7){ for(j in 1:4){ mtx2[i, j] <- mean(mtx[mtx$Niche == c(1:7)[i], j], ) } }
data <- reshape2::melt(mtx2[c(1:4, 6), ])
data$Var2 <- factor(data$Var2, levels = c('CM1', 'CM2', 'C3ar1', 'C3'))
data$Var2 <- revalue(data$Var2, replace = c('CM1'='CM1_marker', 'CM2'='CM2_marker',
                                            'C3ar1'='C3ar1+MP_marker', 'C3'='C3+CF_marker'))
p4 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value))
p4
p <- p3 + p4 &
        scale_fill_viridis_c(limits = c(-1, 1.6)) &
        theme_classic() &
        theme(aspect.ratio = 4/5) &
        labs(fill = 'Z-score', y = 'Marker Exp.', x = 'Niche') &
        scale_y_discrete(limits = rev) &
        RotatedAxis()
p
PlotPDF('01.7.st_wh_split_niche_marker_expr', 10, 5)
p
dev.off()

p <- DimPlot2(st_wh_ctrl.srt, group.by = 'Niche_ctrl', label = T, repel = T, cols = mycol_20, ncol = 4)[[1]] +
        DimPlot2(st_wh_y5sa.srt, group.by = 'Niche_y5sa', label = T, repel = T, cols = mycol_20, ncol = 4)[[1]]
PlotPDF('01.8.st_wh_split_niche_umap', 9, 4)
p
dev.off()

## Combined version:
identical(Cells(st_wh.srt), c(Cells(st_wh_ctrl.srt), Cells(st_wh_y5sa.srt)))
st_wh.srt$Niche_split <- c(paste0('Niche C', as.vector(st_wh_ctrl.srt$Niche_ctrl)),
                           paste0('Niche Y', as.vector(st_wh_y5sa.srt$Niche_y5sa)))
st_wh.srt$Niche_split <- revalue(st_wh.srt$Niche_split, replace = c('Niche C5' = 'x', 'Niche Y5' = 'y'))
st_wh.srt$Niche_split <- revalue(st_wh.srt$Niche_split, replace = c('Niche C6' = 'Niche C5', 'Niche Y6' = 'Niche Y5'))
st_wh.srt$Niche_split <- revalue(st_wh.srt$Niche_split, replace = c('x' = 'Niche C6', 'y' = 'Niche Y6'))
st_wh.srt$Niche_split <- factor(st_wh.srt$Niche_split)
p1 <- SpatialDimPlot(st_wh.srt, group.by = 'Niche_split', image.alpha = 0, stroke = 0) & labs(fill = '')
p1[[1]] <- p1[[1]] + scale_fill_manual(values = c(mycol_20[1:6], 'grey30'))
p1[[2]] <- p1[[2]] + scale_fill_manual(values = c(mycol_20[c(8:9, 12:13, 7, 15)], 'grey70'))
p1
PlotPDF('01.9.st_wh_split_niche', 9, 4)
p1
dev.off()

st_wh_ctrl.srt$Niche_split <- c(paste0('Niche C', as.vector(st_wh_ctrl.srt$Niche_ctrl)))
st_wh_y5sa.srt$Niche_split <- c(paste0('Niche Y', as.vector(st_wh_y5sa.srt$Niche_y5sa)))
st_wh_ctrl.srt$Niche_split <- revalue(st_wh_ctrl.srt$Niche_split, replace = c('Niche C5' = 'x'))
st_wh_ctrl.srt$Niche_split <- revalue(st_wh_ctrl.srt$Niche_split, replace = c('Niche C6' = 'Niche C5'))
st_wh_ctrl.srt$Niche_split <- revalue(st_wh_ctrl.srt$Niche_split, replace = c('x' = 'Niche C6'))
st_wh_y5sa.srt$Niche_split <- revalue(st_wh_y5sa.srt$Niche_split, replace = c('Niche Y5' = 'x'))
st_wh_y5sa.srt$Niche_split <- revalue(st_wh_y5sa.srt$Niche_split, replace = c('Niche Y6' = 'Niche Y5'))
st_wh_y5sa.srt$Niche_split <- revalue(st_wh_y5sa.srt$Niche_split, replace = c('x' = 'Niche Y6'))

p2.1 <- DimPlot2(st_wh_ctrl.srt, group.by = 'Niche_split', label = F, repel = T,
               cols = c(mycol_20[1:6], 'grey30'))[[1]] +
        labs(title = 'Control', x = 'UMAP1', y = 'UMAP2')

p2.2 <- DimPlot2(st_wh_y5sa.srt, group.by = 'Niche_split', label = F, repel = T,
                 cols = c(mycol_20[c(8:9, 12:13, 7, 15)], 'grey70'))[[1]] +
        labs(title = 'YAP5SA', x = 'UMAP1', y = 'UMAP2')
p2 <- p2.1 + p2.2
PlotPDF('01.10.st_wh_split_niche_umap', 9, 4)
p2
dev.off()

mtx <- st_wh.srt@meta.data[, grepl('^Decon|^Niche_split', colnames(st_wh.srt@meta.data))]
dim(mtx)
mtx <- mtx[! is.na(mtx$Decon_CF),]
mtx <- mtx[! mtx$Niche_split %in% c('Niche C6', 'Niche C7', 'Niche Y6', 'Niche Y7'), ]
Table(mtx$Niche_split)
mtx2 <- matrix(NA, 10, 7)
colnames(mtx2) <- str_split(colnames(mtx[1:7]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))
for(i in 1:10){ for(j in 1:7){ mtx2[i, j] <-
        mean(mtx[mtx$Niche_split == c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))[i], j], ) }
        }
for(i in 1:7){mtx2[, i] <- scale(mtx2[, i])}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('Mac', 'SMC', 'EC', 'CM1', 'CM2', 'CF', 'EndoC'))
data$value <- Range01(data$value)
p1 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c() +
        theme_classic() +
        theme(aspect.ratio = 7/10) +
        labs(fill = 'Relative Abundance', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p1
PlotPDF('01.11.st_wh_split_niche_cell_type_relative_abundance', 10, 5)
p1
dev.off()

## Check Triad
mtx <- st_wh.srt@meta.data[, grepl('^Zscore|^Niche_split', colnames(st_wh.srt@meta.data))]
dim(mtx)
mtx <- mtx[! mtx$Niche_split %in% c('Niche C6', 'Niche C7', 'Niche Y6', 'Niche Y7'), ]
Table(mtx$Niche_split)
mtx2 <- matrix(NA, 10, 4)
colnames(mtx2) <- str_split(colnames(mtx[2:5]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))
for(i in 1:10){ for(j in 1:4){
        mtx2[i, j] <- mean(mtx[mtx$Niche_split == c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))[i], j+1], ) }
        }
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('CM1', 'CM2', 'C3ar1', 'C3'))
data$Var2 <- revalue(data$Var2, replace = c('CM1'='CM1_marker', 'CM2'='CM2_marker',
                                            'C3ar1'='C3ar1+MP_marker', 'C3'='C3+CF_marker'))
p3 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c(limits = c(-1, 1.6)) +
        theme_classic() +
        theme(aspect.ratio = 4/10) +
        labs(fill = 'Z-score', y = 'Marker Exp.', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p3
PlotPDF('01.12.st_wh_split_niche_marker_expr', 10, 5)
p3
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(st_wh_ctrl.srt, 'analysis/part91.ST_YAP5SA_wholeheart_ctrl_niche_annotated.seurat.rds')
saveRDS(st_wh_y5sa.srt, 'analysis/part91.ST_YAP5SA_wholeheart_y5sa_niche_annotated.seurat.rds')
saveRDS(st_wh.srt, 'analysis/part91.ST_YAP5SA_wholeheart_niche_annotated.seurat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #2. Reviewer #1 Q2  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Idents(drop.srt) <- 'Cell_type2'
# ct_marker <- FindAllMarkers(drop.srt[, !is.na(drop.srt$Sample)],
#                             logfc.threshold = 0.5, only.pos = T, return.thresh = 0.001)
# ct_marker <- ct_marker[ct_marker$pct.1-ct_marker$pct.2 >= 0.3, ]
# top_ct_marker <- ct_marker |> group_by(cluster) |> top_n(wt = avg_log2FC, n = 30)
# top_ct_marker <- top_ct_marker[!duplicated(top_ct_marker$gene),]
#
# shown_genes <- c(
#         'Sorbs2', 'Flnc', 'Ahnak',
#         'Ndufv2', 'Etfa', 'Atp5j',
#         'Rgs5', 'Prkg1', 'Fstl1',
#         'Smoc2', 'Lyz2', 'Plcb4',
#         'Nrp1', 'Tm4sf1',
#         'Igfbp5', 'Mgp')
# shown_genes[!shown_genes %in% top_ct_marker$gene]

p1.1 <- DimPlot2(drop.srt, group.by = 'Cell_type', reduction = 'umap',
                 cols = Color_cell_type, label = T, raster = T) +
        labs(title = '', x = 'UMAP 1', y = 'UMAP 2')
p1.2 <- CountCellBarPlot(drop.srt[, !is.na(drop.srt$Sample)],
                         group.var = 'Sample', stack.var = 'Cell_type',
                         percentage = T, stack.color = Color_cell_type, width = 0.8)$plot +
        labs(x = 'Samples', y = 'Fraction of All Cells', fill = 'Cell Type') +
        theme(aspect.ratio = 1)

p1.1 | p1.2
PlotPDF('02.1.global_umap_sample_dist', 8, 4)
p1.1 | p1.2
dev.off()

meta <- U(drop_cm.srt@meta.data[! is.na(drop_cm.srt$Sample), c('Sample', 'CM_State')])
meta <- meta[order(meta$Sample),]
mtx <- Table(drop_cm.srt$Sample, drop_cm.srt$CM_State)
mtx <- as.matrix(mtx[, colSums(mtx)>0])
mtx <- mtx/as.vector(Table(drop_cm.srt$Sample))
data1 <- data.frame(Frac = mtx[1:6, 'CM1'],
                 Cell_state = 'aCM1',
                 Genotype = rep(c('Control','YAP5SA'), each = 3))
data2 <- data.frame(Frac = mtx[1:6, 'CM2'],
                    Cell_state = 'aCM2',
                    Genotype = rep(c('Control','YAP5SA'), each = 3))
p1.3 <- ggplot(data1, aes(y = Frac, x = Genotype, color = Genotype)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'aCM1', y = 'Fraction of CMs', x = '')
p1.4 <- ggplot(data2, aes(y = Frac, x = Genotype, color = Genotype)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'aCM2', y = 'Fraction of CMs', x = '')
p1.5 <- p1.3 + p1.4 &
        scale_y_continuous(limits = c(0.2, 0.8)) &
        theme_classic() &
        theme(aspect.ratio = 2) &
        RotatedAxis() &
        NoLegend()
p1.5
PlotPDF('02.2.cm1_cm2_dist_across_sample', 4, 4)
p1.5
dev.off()

df <- CountCellBarPlot(drop.srt, group.var = 'Cell_type', stack.var = 'Sample', percentage = T)$data
df$Genotype <- str_split(df$StackVar, pattern = ' ', simplify = T)[, 1]
df <- df[! is.na(df$Genotype), ]
df$Fraction <- df$Count/Table(drop.srt$Sample)[df$StackVar]

p1.6 <- ggplot(df, aes(x = Genotype, y = Fraction, color = Genotype)) +
        geom_boxplot() +
        geom_point(color = 'black') +
        theme_classic() +
        #scale_y_continuous(minor_breaks = seq(-2, 3.5, 0.5), breaks = seq(-2, 3.5, 1), limits = c(-2, 3.5)) +
        theme(aspect.ratio = 3, panel.grid.major.y = element_line()) +
        RotatedAxis() +
        labs(x = '', y = 'Fraction', color = 'Genotype') +
        NoLegend() +
        facet_wrap(~GroupVar)
p1.6
PlotPDF('02.3.celltype_composition_across_samples', 3, 5)
p1.6
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #3. Reviewer #1 Q3  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For Hassan to answer Q3:
x <- readRDS('~/Downloads/mouse_v1.part01.ST_YAP5SA_wholeheart.srt.rds')
meta <- read.csv(paste0('~/Documents/Bioinformatics/project/2022_yap5sa_st_rli/meta/mouse_v0/91_NCR_Revision_1/',
                        'part91.niche_decon_meta_table.csv'))
meta$Cell_ID <- str_replace(meta$Cell_ID, pattern = '^Control', replacement = 'D005S001')
meta$Cell_ID <- str_replace(meta$Cell_ID, pattern = '^YAP5SA', replacement = 'D005S002')
meta$Coord_x <- x$Coord_x_slide
meta$Coord_y <- x$Coord_y_slide

jong_anno_ctrl.csv <- read.csv('external/control.csv')
jong_anno_ctrl.csv$Cell_ID <- paste0('D005S001_', jong_anno_ctrl.csv$Barcode)
jong_anno_ctrl.csv$Genotype <- 'WT'
jong_anno_y5sa.csv <- read.csv('external/y5sa.csv')
jong_anno_y5sa.csv$Cell_ID <- paste0('D005S002_', jong_anno_y5sa.csv$Barcode)
jong_anno_y5sa.csv$Genotype <- 'YAP5SA'
jong_anno.csv <- rbind(jong_anno_ctrl.csv, jong_anno_y5sa.csv)
jong_anno.csv$Manual[jong_anno.csv$Manual == ''] <- 'Unanno'

O(meta$Cell_ID, jong_anno.csv$Cell_ID)
rownames(jong_anno.csv) <- jong_anno.csv$Cell_ID
meta$Region <- jong_anno.csv[meta$Cell_ID, 'Manual']

x$Region <- jong_anno.csv[Cells(x), 'Manual']
SpatialDimPlot(x, group.by = 'Region', image.alpha = 0, stroke = 0) &
        scale_fill_manual(values = mycol_20)
p <- ggplot(meta, aes(x = Coord_x, y = Coord_y, color = Region)) +
        scale_color_manual(values = mycol_20) +
        geom_point() +
        facet_wrap(~Genotype) +
        theme_classic() +
        theme(aspect.ratio = 1)
PlotPDF('03.1.st_region_annotation', 12, 5)
p
dev.off()

O(Cells(x), Cells(st_wh.srt))
tmp <- str_replace(Cells(x), pattern = 'D005S001_', replacement = 'Control_')
tmp <- str_replace(tmp, pattern = 'D005S002_', replacement = 'YAP5SA_')
identical(Cells(st_wh.srt), tmp)
x <- RenameCells(x, new.names = Cells(st_wh.srt))
x <- AddMetaData(x, metadata = st_wh.srt@meta.data[, c('Niche_split', 'Zscore_CM1', 'Zscore_CM2',
                                                       'Zscore_C3ar1_Mac', 'Zscore_C3_FB')])

#### Try colocalization ####
ct <- str_split(grep('^Decon', colnames(st.srt@meta.data), value = T), pattern = '_', simplify = T)[, 2]
mtx1 <- matrix(NA, 7, 7)
rownames(mtx1) <- ct
colnames(mtx1) <- ct
mtx2 <- mtx1

for(i in 1:L(ct)){
        for(j in 1:L(ct)){
                if(i !=j ){
                        x <- GetColocalProb(st.srt, meta_features = paste0('Decon_', c(ct[i], ct[j])))
                        y <- x>-log10(0.05)
                        z <- split(y, st.srt$Sample)
                        mtx1[i, j] <- sum(z$Control)/L(z$Control)
                        mtx2[i, j] <- sum(z$YAP5SA)/L(z$YAP5SA)
                }
        }
}
data1 <- reshape2::melt(mtx1)
data1$Sample <- 'Control'
data2 <- reshape2::melt(mtx2)
data2$Sample <- 'YAP5SA'
data <- rbind(data1, data2)
data$value[data$value > 0.2] <- 0.2
p <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c(na.value = 0) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        facet_wrap(~Sample)
PlotPDF('03.2.heatamp.st_colocal_pairwise', 10, 5)
p
dev.off()

st.srt$Colocal_CM1_EC <- GetColocalProb(st.srt, meta_features = c('Decon_EC', 'Decon_CM1'))
st.srt$Colocal_EndoC_MP <- GetColocalProb(st.srt, meta_features = c('Decon_EndoC', 'Decon_Mac'))
st.srt$Colocal_CM2_MP <- GetColocalProb(st.srt, meta_features = c('Decon_CM2', 'Decon_Mac'))

p <- FeaturePlotST(srt = st.srt,
                   features = c('Colocal_CM1_EC',
                                'Colocal_EndoC_MP',
                                'Colocal_CM2_MP'),
                   title = c('EC-CM1',
                             'EndoC_MP',
                             'CM2_MP'),
                   minvals = rep(0, 6), maxvals = rep(3, 6), pt.sizes = st.srt@misc$spot_scale*0.35,
                   ncol = 2) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 3), direction = 0)
p[[2]] <- p[[2]] + RestoreLegend() + labs(color='-Log10(p)')
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[6]] <- p[[6]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('03.3.st_colocal_pairwise_example', 8, 12)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(x, 'analysis/part91.ST_YAP5SA_wholeheart_niche_annotated.FOR_HASSAN.seurat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #4. Reviewer #1 Q6  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Calculate Y5SA dataset CM1/2-C3-C3ar1 scores ####
st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = T,
                          names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)))

st.srt <- AddModuleScore2(st.srt, features = Yap_target[1], assay = 'SCT', return_z = T,
                          names = c('Zscore_Y5SA_Target'))

st.srt <- AddModuleScore2(st.srt, features = list('C3', 'C3ar1'), assay = 'SCT', return_z = T,
                names = paste0('Zscore_', c('C3', 'C3ar1')))

RidgePlot(st.srt, features = c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3', 'Zscore_C3ar1',
                               'Zscore_Y5SA_Target'))

#### Calculate Y5SA dataset CM1/2-C3-C3ar1 colocalization ####
st.srt$Colocal_CM2_C3ar1_C3 <- GetColocalProb(st.srt, meta_features = c('Zscore_CM2', 'Zscore_C3', 'Zscore_C3ar1'))
st.srt$Colocal_CM2_C3 <- GetColocalProb(st.srt, meta_features = c('Zscore_CM2', 'Zscore_C3'))
p <- FeaturePlotST(srt = st.srt,
                   features = c('Colocal_CM2_C3',
                                'Colocal_CM2_C3ar1_C3'),
                   title = c('CM2 + C3 Hi',
                             'CM2 + C3ar1 Hi + C3 Hi'),
                   minvals = rep(0, 4), maxvals = rep(3, 4), pt.sizes = st.srt@misc$spot_scale*0.35,
                   ncol = 2) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 3), direction = 0)
p[[2]] <- p[[2]] + RestoreLegend() + labs(color='-Log10(p)')
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('04.1.feat_st.cm2_c3_c3ar1_colocal', 6, 5)
print(p)
dev.off()

#### Quantify Y5SA dataset MC2-C3-C3ar1 colocalization ####
st.srt$Psig_CM2_C3 <- factor('p<0.05', levels = c('Not sig.', 'p<0.05'))
st.srt$Psig_CM2_C3_C3ar1 <- factor('p<0.05', levels = c('Not sig.', 'p<0.05'))
st.srt$Psig_CM2_C3[st.srt$Colocal_CM2_C3 < -log10(0.05)] <- 'Not sig.'
st.srt$Psig_CM2_C3_C3ar1[st.srt$Colocal_CM2_C3ar1_C3 < -log10(0.05)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st.srt$Psig_CM2_C3, st.srt$Sample)),
            as.data.frame(Table(st.srt$Psig_CM2_C3_C3ar1, st.srt$Sample)))
df$Group <- rep(c('CM2_C3', 'CM2_C3_C3ar1'), each = 2*LU(st.srt$Sample))
df$Fraction = df$Freq/table(st.srt$Sample)[df$Var2]
p <- ggplot(df[df$Var1=='p<0.05',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~Group) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('04.2.bar.cm2_c3_c3ar1_colocolization_quantification', 4, 5)
print(p)
dev.off()

meta <- st.srt@meta.data[, c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB',
                             'Zscore_Y5SA_Target', 'Zscore_C3', 'Zscore_C3ar1',
                             'Colocal_CM2_C3ar1_C3', 'Colocal_CM2_C3',
                             'Psig_CM2_C3', 'Psig_CM2_C3_C3ar1')]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(meta, 'analysis/part91.ST_YAP5SA_venticle_colocal_meta.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### Calculate Olson dataset CM1/2-C3CF-C3ar1MP scores ####
st2.srt <- readRDS('individual/part01.2021_NatComm_EOlson.srt.rds')
st2.srt <- AddModuleScore2(st2.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'ST', return_z = T,
                           names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)))
st2.srt <- AddModuleScore2(st2.srt, features = Yap_target[1], assay = 'ST', return_z = T,
                           names = 'Zscore_Target')
st2.srt$TMP <- 'ALL'
RidgePlot(st2.srt, features = c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB', 'Zscore_Target'),
          group.by = 'TMP')

minvals <- rep(-1, 5)
maxvals <- rep(2, 5)
scale <- c(2.48, 4, 1.3, 2.25)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB', 'Zscore_Target'),
                   title = c('CM1 score', 'CM2 score', 'C3ar1+ Mac score', 'C3+ FB score', 'Yap target score'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.4,
                   ncol = 4)
p[[20]] <- p[[20]] + RestoreLegend() + labs(color='Z score')
p
PlotPDF('04.3.feat_st.olson_st_yap5sa_marker_zscore', 12, 10)
print(p)
dev.off()

#### Calculate Olson dataset CM1/2-C3CF-C3ar1MP colocalization ####
st2.srt$Colocal_CM1_C3ar1Mac_C3FB <- GetColocalProb(st2.srt,
                                                    meta_features = c('Zscore_CM1', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB'))
st2.srt$Colocal_CM2_C3ar1Mac_C3FB <- GetColocalProb(st2.srt,
                                                    meta_features = c('Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB'))

minvals <- rep(0, 4)
maxvals <- rep(5, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_C3ar1Mac_C3FB',
                                'Colocal_CM2_C3ar1Mac_C3FB'),
                   title = c('CM1, C3ar1+ Mac, C3+ FB',
                             'CM2, C3ar1+ Mac, C3+ FB'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.6,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 5), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('04.4.feat_st.olson_st_yap5sa_marker_colocolization', 12, 5)
print(p)
dev.off()

#### Quantify Olson dataset CM1/2-C3CF-C3ar1MP colocalization ####
st2.srt$Psig_CM1_C3CF_C3ar1MP <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM2_C3CF_C3ar1MP <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM1_C3CF_C3ar1MP[st2.srt$Colocal_CM1_C3ar1Mac_C3FB < -log10(0.01)] <- 'Not sig.'
st2.srt$Psig_CM2_C3CF_C3ar1MP[st2.srt$Colocal_CM2_C3ar1Mac_C3FB < -log10(0.01)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1_C3CF_C3ar1MP, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2_C3CF_C3ar1MP, st2.srt$Sample)))
df$Group <- rep(c('CM1_C3CF_C3ar1MP', 'CM2_C3CF_C3ar1MP'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.01',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~Group) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('04.4.bar.olson_st_yap5sa_marker_colocolization_on_olson_st', 4, 5)
print(p)
dev.off()

#### Calculate Olson dataset CM1/2-YapTarget colocalization ####
st2.srt$Colocal_CM1_YapTarget <- GetColocalProb(st2.srt, meta_features = c('Zscore_CM1', 'Zscore_Target'))
st2.srt$Colocal_CM2_YapTarget <- GetColocalProb(st2.srt, meta_features = c('Zscore_CM2', 'Zscore_Target'))

minvals <- rep(0, 4)
maxvals <- rep(5, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_YapTarget',
                                'Colocal_CM2_YapTarget'),
                   title = c('CM1, YapHigh',
                             'CM2, YapHigh'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.6,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 5), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('04.5.feat_st.olson_st_yap5sa_marker_yap_target_colocolization', 12, 5)
print(p)
dev.off()

#### Quantify Olson dataset CM1/2-YapTarget colocalization ####
st2.srt$Psig_CM1_YapHi <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM2_YapHi <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM1_YapHi[st2.srt$Colocal_CM1_YapTarget < -log10(0.01)] <- 'Not sig.'
st2.srt$Psig_CM2_YapHi[st2.srt$Colocal_CM2_YapTarget < -log10(0.01)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1_YapHi, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2_YapHi, st2.srt$Sample)))
df$Group <- rep(c('CM1_YapHi', 'CM2_YapHi'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.01',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~Group) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('04.6.bar.olson_st_yap5sa_marker_yap_target_colocolization', 4, 5)
print(p)
dev.off()

meta <- st2.srt@meta.data[, c(
        'Zscore_CM1', 'Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB', 'Zscore_Target',
        'Colocal_CM1_C3ar1Mac_C3FB', 'Colocal_CM2_C3ar1Mac_C3FB',
        'Psig_CM1_C3CF_C3ar1MP', 'Psig_CM2_C3CF_C3ar1MP',
        'Colocal_CM1_YapTarget', 'Colocal_CM2_YapTarget',
        'Psig_CM1_YapHi', 'Psig_CM2_YapHi')
]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(meta, 'analysis/part91.ST_Olson_colocal_meta.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## For hassan:
data <- data.frame(Spot_ID = names(sort(st2.srt$Zscore_CM1)),
                   CDF_cm1 = sort(st2.srt$Zscore_CM1),
                   CDF_cm2 = st2.srt$Zscore_CM2[names(x)],
                   CDF_yaptarget = st2.srt$Zscore_Target[names(x)])
WriteCSV(data, 'olson_st_cm1_cm2_yaptarget_cdf_values')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #5. Reviewer #1 Q8  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Try Velocyto == Failed due to too few splicing events  ####

## Add CM sub umap
O(Cells(drop.srt), Cells(drop_cm.srt))
drop.srt@reductions$sub_umap <- drop.srt@reductions$umap
drop.srt@reductions$sub_umap@cell.embeddings[,] <- NA
colnames(drop.srt@reductions$sub_umap@cell.embeddings) <- c('SUBUMAP_1', 'SUBUMAP_2')
drop.srt@reductions$sub_umap@cell.embeddings[Cells(drop_cm.srt), 1:2] <-
        drop_cm.srt@reductions$newumap@cell.embeddings[, 1:2]
drop.srt@reductions$sub_umap@key <- 'SUBUMAP_'

## Generating valid barcode files for velocyto
bc1 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'Control 1']), pattern = '_', simplify = T)[, 7]
L(bc1)
bc2 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'Control 2']), pattern = '_', simplify = T)[, 8]
L(bc2)
bc3 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'Control 3']), pattern = '_', simplify = T)[, 7]
L(bc3)
bc4 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'YAP5SA 1']), pattern = '_', simplify = T)[, 8]
L(bc4)
bc5 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'YAP5SA 2']), pattern = '_', simplify = T)[, 7]
L(bc5)
bc6 <- str_split(Cells(drop.srt[, drop.srt$Sample == 'YAP5SA 3']), pattern = '_', simplify = T)[, 7]
L(bc6)
write.table(bc1, file = '~/Desktop/ctrl_cm_s1.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(bc2, file = '~/Desktop/ctrl_wh_s1.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(bc3, file = '~/Desktop/ctrl_wh_s2.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(bc4, file = '~/Desktop/yap5sa_cm_s1.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(bc5, file = '~/Desktop/yap5sa_wh_s2.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(bc6, file = '~/Desktop/yap5sa_wh_s1.bc.txt', quote = F, sep = '\t', row.names = F, col.names = F)

## Run velocyto on lorien
## conda activate velocyto
## velocyto run-dropest \
## -b yap5sa_wh_s2.bc.txt \
## -o velocyto_yap5sa_wh_s1/ \
## -e yap5sa_wh_s2 \
## -m /lorien/genome/velocyto/hg38_rmsk.gtf \
## -@ 14 \
## --samtools-memory 5120 \
## yap5sa_wh_s2.clean.bam \
## /lorien/genome/index_cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

## rename cells to match loom files
barcode <- list(bc1, bc2, bc3, bc4, bc5, bc6)
id <- split(Cells(drop.srt),  drop.srt$Sample)
for(i in 1:L(barcode)){
        barcode[[i]] <- paste0(barcode[[i]], '_', i)
        names(barcode[[i]]) <- as.vector(id[[i]])
}
new_name <- unlist(barcode)
head(new_name)
tail(new_name)

drop2.srt <- RenameCells(drop.srt[, !is.na(drop.srt$Sample)], new.names = new_name)

##  Export Seurat data for scVelo
##  Save metadata table:
drop2.srt$barcode <- Cells(drop2.srt)
drop2.srt$UMAP_1 <- drop2.srt@reductions$umap@cell.embeddings[,1]
drop2.srt$UMAP_2 <- drop2.srt@reductions$umap@cell.embeddings[,2]
drop2.srt$SUBUMAP_1 <- drop2.srt@reductions$sub_umap@cell.embeddings[,1]
drop2.srt$SUBUMAP_2 <- drop2.srt@reductions$sub_umap@cell.embeddings[,2]
drop2.srt$group1 <- drop2.srt$Sample
drop2.srt$group2 <- drop2.srt$Experiment
drop2.srt$Cell_state <- drop2.srt$Cell_type2

meta <- drop2.srt@meta.data[, c('group1', 'group2', 'Cell_type', 'Cell_state',
                                'barcode', 'UMAP_1', 'UMAP_2', 'SUBUMAP_1', 'SUBUMAP_2')]
WriteCSV(meta, title = 'full.srt_meta')

# write dimesnionality reduction matrix
harmony <- drop2.srt@reductions$harmony@cell.embeddings
WriteCSV(as.data.frame(harmony), title = 'full_harmony.srt_dim')

# write gene names
WriteCSV(as.data.frame(rownames(drop2.srt)), title = 'full.gene_names', col.names = F)

## Write expression counts matrix ## HUGE FILE!
counts_matrix <- GetAssayData(drop2.srt, assay = 'SCT', slot = 'counts')
writeMM(counts_matrix, file = paste0(Meta_dir, 'full.counts_mtx'))

## Manually gzip all 4 txt files and move to data drive rdata/analysis


#### Try Slingshot  ####
library('slingshot')
sce <- as.SingleCellExperiment(drop_cm.srt, assay = 'SCT')
sce <- slingshot(sce, reducedDim = 'HARMONY', clusterLabels = 'CM_State', start.clus = 'CM1', use.median = T)
# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)

# Plot Slingshot pseudotime vs cell stage.
drop_cm.srt$SlingPseudotime <- sce$slingPseudotime_1
drop_cm.srt$SlingRank <- rank(sce$slingPseudotime_1)
data <- data.frame(PST = rank(sce$slingPseudotime_1), State = sce$CM_State)
p1 <- ggplot(data, aes(x = PST, y = State, colour = State)) +
        geom_quasirandom(groupOnX = F, width = 0.4, size = 0.2) +
        scale_color_manual(values = mycol_10) +
        labs(x = "Pseudotime", y = "CM State", title = "Infered CMs Transition") +
        scale_y_discrete(limits = rev) +
        theme_classic() +
        theme(aspect.ratio = 1)
p2 <- FeaturePlot2(drop_cm.srt, features = 'SlingRank', reduction = 'newumap', pt.size = 0.2) +
        scale_color_distiller(palette = 'Spectral')
p3 <- DimPlot2(drop_cm.srt, group.by = 'CM_State', reduction = 'newumap', cols = mycol_10, pt.size = 0.2)

p <- wrap_plots(p3, p2, p1, ncol = 2)
p
PlotPDF('05.1.slingshot_pseudotime', 10, 10)
print(p)
dev.off()

#### Try DiffusionMap ####
library('destiny')
library('ggbeeswarm')
umi <- logcounts(sce)
State <- sce$CM_State
hpca <- reducedDim(sce, "HARMONY")
dm <- DiffusionMap(hpca)
dpt <- DPT(dm, tips = which.max(eigenvectors(dm)[, 1]))
## tip cell: 2021_Unpublished_JMartin__Ventricule__GCTTTAGGAGGT_YAP5SA.2

drop_cm.srt$DiffDim1 <- rank(eigenvectors(dm)[, 1])
drop_cm.srt$DPT <- rank(-dpt$dpt)
drop_cm.srt@reductions$diffusion <- CreateDimReducObject(embeddings = dm@eigenvectors, assay = 'RNA', key = 'DC')

data <- data.frame(DPT = drop_cm.srt$DPT,
                   DC1 = drop_cm.srt@reductions$diffusion@cell.embeddings[, 1],
                   DC2 = drop_cm.srt@reductions$diffusion@cell.embeddings[, 2],
                   State = factor(drop_cm.srt$CM_State))

p1 <- DimPlot2(drop_cm.srt, reduction = 'diffusion', dims = 1:2, group.by = 'CM_State', cols = mycol_10, pt.size = 0.2)
p2 <- FeaturePlot2(drop_cm.srt, reduction = 'newumap', features = 'DPT', pt.size = 0.2)
p3 <-  ggplot(data, aes(x = DC1, y = DC2, colour = DPT)) +
        geom_point(size = 0.2) +
        scale_color_distiller(palette = 'Spectral') +
        labs(x = "DC 1", y = "DC 2", title = "CMs Colored by Pseudotime", color = 'Pseudotime') +
        theme_classic() +
        theme(aspect.ratio = 1)
p4 <- ggplot(data, aes(x = DPT, y = State, colour = State)) +
        geom_quasirandom(groupOnX = F, width = 0.4, size = 0.2) +
        scale_color_manual(values = mycol_10) +
        labs(x = "Pseudotime", y = "CM State", title = "Infered CMs Transition") +
        scale_y_discrete(limits = rev) +
        theme_classic() +
        theme(aspect.ratio = 1)
p <- wrap_plots(p1, p2, p3, p4, ncol = 2)
p
PlotPDF('05.2.diffusion_pseudotime', 10, 10)
print(p)
dev.off()

PlotPDF('05.3.feat_scat.comparison_sling_vs_dpt', 5, 5)
FeatureScatter(drop_cm.srt, feature1 ='DPT', feature2 = 'SlingRank',
               group.by = 'CM_State', cols = mycol_10, pt.size = 0.2) +
        theme(aspect.ratio = 1)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(sce, 'analysis/part91.drop_cm_slingshot.sce.rds')
saveRDS(dm, 'analysis/part91.drop_cm_diffusion_map.dm.rds')
saveRDS(dpt, 'analysis/part91.drop_cm_diffusion_pseudotime.dpt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Find temporal correlated genes using DPT ####
Y <- drop_cm.srt@assays$SCT@data
var3K <- drop_cm.srt@assays$SCT@var.features[1:3000]
Y <- Y[var3K, ]  # only counts for variable genes

t <- drop_cm.srt$DPT
pcc.df <- data.frame(PCC = rep(NA, 3000), Pval = rep(NA, 3000), row.names = rownames(Y))
for(i in 1:3000){
        pcc.df[i, ] <- c(cor(Y[i,], t),
                         cor.test(Y[i,], t)$p.value)
}
hist(pcc.df$PCC, breaks = 100)
pcc.df <- pcc.df[pcc.df$Pval <= 0.01 & abs(pcc.df$PCC) >= 0.2, ]
module1 <- split(rownames(pcc.df), pcc.df$PCC > 0)
str(module1)
names(module1) <- c('Neg', 'Pos')

## Find temporal correlated genes using Slignshot Pseudotime ####
Y <- drop_cm.srt@assays$SCT@data
var3K <- drop_cm.srt@assays$SCT@var.features[1:3000]
Y <- Y[var3K, ]  # only counts for variable genes

t <- drop_cm.srt$SlingRank
pcc.df <- data.frame(PCC = rep(NA, 3000), Pval = rep(NA, 3000), row.names = rownames(Y))
for(i in 1:3000){
        pcc.df[i, ] <- c(cor(Y[i,], t),
                         cor.test(Y[i,], t)$p.value)
}
pcc.df <- pcc.df[! is.na(pcc.df$PCC), ]
hist(pcc.df$PCC, breaks = 100)
pcc.df2 <- pcc.df[pcc.df$Pval <= 0.01 & abs(pcc.df$PCC) >= 0.2, ]
module2 <- split(rownames(pcc.df2), pcc.df2$PCC > 0)
str(module2)
names(module2) <- c('Neg', 'Pos')

pcc.df3 <- pcc.df2[pcc.df2$PCC > 0, ]
colnames(pcc.df3) <- c('Pearson Correlation Coefficient', 'P value')
pcc.df3$`Gene Symbol` <- rownames(pcc.df3)
pcc.df3 <- pcc.df3[, c(3, 1, 2)]
WriteCSV(pcc.df3,  title = 'part91.cm1_cm2_transition_correlated_genes')

## Plot temporal correlated genes ####
FlattenExpr <- function(expr_scale, genes){
        expr_flat <- expr_scale[genes,]
        expr_flat <- as.data.frame(matrix(expr_flat, ncol = 1))
        expr_flat <- cbind(expr_flat,
                           rep(1:ncol(expr_scale), each = length(genes)),
                           rep(genes, times = ncol(expr_scale)))
        colnames(expr_flat) <- c("Expression", "Cells", "Genes")
        expr_flat <- expr_flat[order(expr_flat$Genes), ]
        expr_flat <- cbind(expr_flat,
                           colMeans(expr_scale[genes, ]),
                           colMedians(expr_scale[genes, ]),
                           colQuantiles(expr_scale[genes, ], probs = 0.75),
                           colQuantiles(expr_scale[genes, ], probs = 0.25))
        colnames(expr_flat)[4:7] <- c("Mean Expr", "Median Expr",
                                      "Upper Quartile", "Lower Quartile")
        return(expr_flat)
}

all_module_genes <-  U(unlist(c(module1, module2)))
drop_cm.srt <- ScaleData(drop_cm.srt, features = all_module_genes, assay = 'SCT')
expr_mat <- drop_cm.srt@assays$SCT@scale.data[all_module_genes, colnames(drop_cm.srt)[order(drop_cm.srt$DPT)]]
for(i in 1:nrow(expr_mat)){expr_mat[i, ] <- smooth.spline(expr_mat[i, ], spar = 1.2)$y}

## Plot DPT modules
df_all <- mapply(function(x, y){
        expr_mat <- t(apply(expr_mat[x,], 1, Range01)) ## normalize to [0-1]
        df <- FlattenExpr(expr_mat[, seq(1, ncol(drop_cm.srt), 1)], x)
        df$DPT <- Range01(df$Cells)
        df$mod_id <- names(module1)[y]
        return(df)},
        x = module1, y = seq_along(module1), SIMPLIFY = F)
df_all <- bind_rows(df_all)
name <- paste0("Module ", 1:2, ' (', lapply(module1, function(x) length(x)), ' genes)')
p <- ggplot(data = df_all[df_all$mod_id == 'Pos', ]) +
        geom_ribbon(aes(x = DPT, ymax = `Upper Quartile`, ymin = `Lower Quartile`, fill = mod_id),
                    alpha = 0.1, show.legend = F) +
        geom_line(aes(x = DPT, y = `Median Expr`, color = mod_id),
                  show.legend = T, size = 1, alpha = 1) +
        labs(y = 'Normalized expression') +
        scale_color_manual(values = mycol_10, labels = name) +
        scale_fill_manual(values = mycol_10, labels = name)  +
        scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
        theme_classic() +
        theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.text.x = element_blank(),
              text = element_text(color = "black"),
              line = element_line(color = "black"),
              axis.line = element_line(color = "black"),
              legend.key.size = unit(10, 'pt'),
              legend.text = element_text(size = 10),
              legend.box.margin = margin(rep(1,1)), legend.title = element_blank(),
              legend.direction = 'vertical', legend.position = 'right',
              #panel.grid.major.x = element_line(color = 'black', linetype = 'dotted')
        )
p1
PlotPDF('05.4.trend.cm1_cm2_diffusion_pseudotime_module', 4, 4)
p1
dev.off()

## Plot Slingshot modules
df_all <- mapply(function(x, y){
        expr_mat <- t(apply(expr_mat[x,], 1, Range01)) ## normalize to [0-1]
        df <- FlattenExpr(expr_mat[, seq(1, ncol(drop_cm.srt), 1)], x)
        df$SlingRank <- Range01(df$Cells)
        df$mod_id <- names(module2)[y]
        return(df)},
        x = module2, y = seq_along(module2), SIMPLIFY = F)
df_all <- bind_rows(df_all)
name <- paste0("Module ", 1:2, ' (', lapply(module2, function(x) length(x)), ' genes)')
p2 <- ggplot(data = df_all[df_all$mod_id == 'Pos', ]) +
        geom_ribbon(aes(x = SlingRank, ymax = `Upper Quartile`, ymin = `Lower Quartile`),
                    alpha = 0.1, show.legend = F) +
        geom_line(aes(x = SlingRank, y = `Median Expr`),
                  show.legend = T, size = 1, alpha = 1) +
        labs(y = 'Normalized expression', x = 'aCM1 to aCM2 Pseudotime',
             title = 'Temporally-correlated Genes', caption = '62 Genes') +
        # scale_color_manual(values = mycol_10, labels = name) +
        # scale_fill_manual(values = mycol_10, labels = name)  +
        scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
        theme_classic() +
        theme(aspect.ratio = 0.7,
              axis.text.x = element_blank(),
              text = element_text(color = "black"),
              line = element_line(color = "black"),
              axis.line = element_line(color = "black"),
              axis.ticks.x.bottom = element_blank()
              #panel.grid.major.x = element_line(color = 'black', linetype = 'dotted')
        )
p2
PlotPDF('05.5.trend.cm1_cm2_slingshot_module', 4, 4)
p2
dev.off()

data <- data.frame(PST = rank(sce$slingPseudotime_1), State = sce$CM_State)
p1 <- ggplot(data, aes(x = PST, y = State, colour = State)) +
        geom_quasirandom(groupOnX = F, width = 0.4, size = 0.01) +
        scale_color_manual(values = mycol_10) +
        labs(x = "aCM1 to aCM2 Pseudotime", y = "") +
        scale_y_discrete(limits = rev) +
        theme_classic() +
        NoLegend() +
        theme(aspect.ratio = 0.2,
              axis.text.x = element_blank(),
              text = element_text(color = "black"),
              line = element_line(color = "black"),
              axis.line = element_line(color = "black"),
              )
p2/p1
PlotPDF('05.6.trend_beeswarm.cm1_cm2_slingshot_module', 4, 4)
p2/p1
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveRDS(list('DPT' = module1, 'Slingshot' = module2), 'analysis/part91.cm1_cm2_transition_modules.list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #6. Reviewer #2 Q5  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cm2 <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$CM2, from = 'mouse', 'human')
cm2 <- revalue(cm2, replace = c('CCN2' = 'CTGF'))

hs_ref_cm.srt <- readRDS('external/fransico_data/reference_scrnaseq/2020_Nature_STeichmann.cm.srt.rds')
hs_ref_cm.srt <- NormalizeData(hs_ref_cm.srt, verbose = F)
Idents(hs_ref_cm.srt) <- 'pub_id'
hs_ref_cm.srt <- subset(hs_ref_cm.srt, idents = c("Ventricular_Cardiomyocyte")) ## exclude for atrial CMs
hs_ref_cm.srt <- FindNeighbors(hs_ref_cm.srt, reduction = "pca", dims = 1:50) |>
        FindClusters(resolution = 0.5)

hs_ref_cm.srt <- AddModuleScore2(hs_ref_cm.srt, features = list(cm2), names = 'Score_aCM2', return_z = T)

p <- BoxPlot(hs_ref_cm.srt, group.by = 'seurat_clusters', feature = 'Score_aCM2', cols = rep('grey90', 14)) +
        NoLegend() +
        theme(axis.text.x = element_text())
p
PlotPDF('06.1.box.teichmann_vCM_cluster_cm2_score', 4, 4)
p
dev.off()

saveRDS(hs_ref_cm.srt@meta.data, 'analysis/part91.Teichmann_vCM_clustered.srt_meta.rds')

## Get CM8 Markers
Idents(hs_ref_cm.srt) <- 'seurat_clusters'
deg <- FindMarkers(hs_ref_cm.srt,
                   ident.1 = '8',
                   min.pct = 0.25,
                   random.seed = 123,
                   logfc.threshold = 0.1)
deg <- deg[order(deg$avg_log2FC, decreasing = T), ]
deg$gene <- rownames(deg)
saveRDS(deg, 'analysis/part91.Teichmann_CM8_marker.srt_marker.rds')

deg$DEG_FOR_GO <- NA
deg$DEG_FOR_GO[deg$avg_log2FC >  0.6 & deg$p_val < 0.05] <- 'CM8_Up'
deg$DEG_FOR_GO[deg$avg_log2FC < -0.6 & deg$p_val < 0.05] <- 'CM8_Down'
WriteCSV(deg, 'part91.teichmann_cm8_vs_other_vcm_deg')

Up <- split(deg$gene, deg$DEG_FOR_GO)[['CM8_Up']]
Dn <- split(deg$gene, deg$DEG_FOR_GO)[['CM8_Down']]

terms <- ModuleEnrichment(module_list = list(Up = Up, Dn = Dn), human_or_mouse = 'human')
go_up <- terms$GO$GO_Up
go_dn <- terms$GO$GO_Dn

WriteCSV(go_up, 'part91.teichmann_cm8_up_enriched_go_bp')
WriteCSV(go_dn, 'part91.teichmann_cm8_dn_enriched_go_bp')

drop_cm_deg <- FindMarkers(drop_cm.srt, ident.1 = 'YAP5SA', ident.2 = 'Control')
drop_cm_deg <- drop_cm_deg[drop_cm_deg$p_val_adj < 0.05, ]
Table(drop_cm_deg$avg_log2FC>0)
CM1 <- rownames(drop_cm_deg)[drop_cm_deg$avg_log2FC < 0]
CM2 <- rownames(drop_cm_deg)[drop_cm_deg$avg_log2FC > 0]
CM1 <- ConvertGeneSpecies(CM1, from = 'mouse', 'human')
CM2 <- ConvertGeneSpecies(CM2, from = 'mouse', 'human')

S(intersect(Dn, CM1))
S(intersect(Up, CM2))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #7. Reviewer #3 Q0 (Re-run inference)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('CellChat')
## Build LR db
CellChatDB.my <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
CellChatDB.my$interaction[2022:2025, ] <- NA
CellChatDB.my$interaction$interaction_name[2022:2025] <- c('ADAM15-ITGB1','ADAM15-ITGB3','ADAM15-ITGA5','ADAM15-ITGA9')
rownames(CellChatDB.my$interaction)[2022:2025] <- c('ADAM15-ITGB1', 'ADAM15-ITGB3', 'ADAM15-ITGA5', 'ADAM15-ITGA9')
CellChatDB.my$interaction$pathway_name[2022:2025] <- 'ADAM15'
CellChatDB.my$interaction$ligand[2022:2025] <- 'Adam15'
CellChatDB.my$interaction$receptor[2022:2025] <- c('Itgb1', 'Itgb3', 'Itga5', 'Itga9')
CellChatDB.my$interaction$agonist[2022:2025] <- ''
CellChatDB.my$interaction$antagonist[2022:2025] <- ''
CellChatDB.my$interaction$co_A_receptor[2022:2025] <- ''
CellChatDB.my$interaction$co_I_receptor[2022:2025] <- ''
CellChatDB.my$interaction$evidence[2022:2025] <- ''
CellChatDB.my$interaction$annotation[2022:2025] <- ''
CellChatDB.my$interaction$interaction_name_2[2022:2025] <-c('Adam15-Itgb1','Adam15-Itgb3','Adam15-Itga5','Adam15-Itga9')
c('Adam15', 'Itgb1', 'Itgb3', 'Itga5', 'Itga9') %in% CellChatDB.my$geneInfo$Symbol

DoCellChat <- function(srt, group.by, CellChatDB = CellChatDB.mouse, LR.type = 'all', species = 'mouse', trim = 0.1){
        # use CellChatDB.human if running on human data
        if(LR.type == 'all') {CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling",
                                                                                "Cell-Cell Contact",
                                                                                "ECM-Receptor"))
        } else {CellChatDB.use <- subsetDB(CellChatDB, search = LR.type)}
        ####  Build CellChat Object
        cch <- createCellChat(object = srt, group.by = group.by)
        cch <- addMeta(cch, meta = srt@meta.data)
        groupSize <- as.numeric(table(cch@idents))
        cch@DB <- CellChatDB.use # set the used database in the object
        ####  Preprocessing the expression data for cell-cell communication analysis
        cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
        cch <- identifyOverExpressedGenes(cch)
        cch <- identifyOverExpressedInteractions(cch)
        if(species == 'mouse') {
                cch <- projectData(cch, PPI.mouse)
        } else if(pecies == 'human') {
                cch <- projectData(cch, PPI.human)
        } else {stop('Species not supported')}
        ####  Inference of cell-cell communication network
        cch <- computeCommunProb(cch, type = "truncatedMean", trim = trim)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}

## Run Cell Chat on RNA data ####
drop.ctrl.srt <- drop.srt[, drop.srt$Sample %in% paste('Control', 1:3)]
drop.y5sa.srt <- drop.srt[, drop.srt$Sample %in% paste('YAP5SA', 1:3)]
DefaultAssay(drop.ctrl.srt) <- 'RNA'
drop.ctrl.srt <- DietSeurat(drop.ctrl.srt, assays = 'RNA')
drop.ctrl.srt$Cell_type2 <- factor(drop.ctrl.srt$Cell_type2,
                                   levels = c('aCM2', 'aCM1', 'CF', 'MP', 'EC1', 'EC2', 'Mural'))
drop.ctrl.cch <- DoCellChat(
        drop.ctrl.srt,
        group.by = 'Cell_type2',
        CellChatDB = CellChatDB.mouse,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)
DefaultAssay(drop.y5sa.srt) <- 'RNA'
drop.y5sa.srt <- DietSeurat(drop.y5sa.srt, assays = 'RNA')
drop.y5sa.srt$Cell_type2 <- factor(drop.y5sa.srt$Cell_type2,
                                   levels = c('aCM2', 'aCM1', 'CF', 'MP', 'EC1', 'EC2', 'Mural'))
drop.y5sa.cch <- DoCellChat(
        drop.y5sa.srt,
        group.by = 'Cell_type2',
        CellChatDB = CellChatDB.mouse,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)
## Reproducing Francisco's result  ####
PlotPDF('07.1.cellchat.mp_to_cm1_cm2', 6, 5)
netVisual_chord_gene(drop.ctrl.cch, sources.use = c('aCM1', 'aCM2'), targets.use = c('MP'),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Control',
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
netVisual_chord_gene(drop.y5sa.cch, sources.use = c('aCM1', 'aCM2'), targets.use = c('MP'),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'YAP5SA',
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
PlotPDF('07.2.cellchat.mp_to_cf', 6, 5)
netVisual_chord_gene(drop.ctrl.cch, sources.use = c('CF'), targets.use = c('MP'),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Control',
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
netVisual_chord_gene(drop.y5sa.cch, sources.use = c('CF'), targets.use = c('MP'),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'YAP5SA',
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()

## Split C3+/- CF ####
tmp.srt <- drop.srt[, drop.srt$Cell_type2 %in% c('aCM1', 'aCM2', 'CF', 'MP')]
DefaultAssay(tmp.srt) <- 'RNA'
tmp.srt <- DietSeurat(tmp.srt, assays = 'RNA')
tmp.srt$Cell_type2 <- droplevels(tmp.srt$Cell_type2)
tmp.srt$Cell_type3 <- factor(tmp.srt$Cell_type2, levels = c('aCM2', 'aCM1', 'C3+CF', 'C3-CF', 'C3ar1+MP', 'C3ar1-MP'))
tmp.srt$Cell_type3[tmp.srt@assays$RNA@data['C3', ] >= 0.2 & tmp.srt$Cell_type2 == 'CF'] <- 'C3+CF'
tmp.srt$Cell_type3[tmp.srt@assays$RNA@data['C3', ] < 0.2 & tmp.srt$Cell_type2 == 'CF'] <- 'C3-CF'
tmp.srt$Cell_type3[tmp.srt@assays$RNA@data['C3ar1', ] >= 0.2 & tmp.srt$Cell_type2 == 'MP'] <- 'C3ar1+MP'
tmp.srt$Cell_type3[tmp.srt@assays$RNA@data['C3ar1', ] < 0.2 & tmp.srt$Cell_type2 == 'MP'] <- 'C3ar1-MP'
Table(tmp.srt$Cell_type3, tmp.srt$Experiment)

drop.ctrl.srt <- tmp.srt[, tmp.srt$Sample %in% paste('Control', 1:3)]
drop.y5sa.srt <- tmp.srt[, tmp.srt$Sample %in% paste('YAP5SA', 1:3)]
drop.ctrl2.cch <- DoCellChat(
        drop.ctrl.srt,
        group.by = 'Cell_type3',
        CellChatDB = CellChatDB.my,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)
drop.y5sa2.cch <- DoCellChat(
        drop.y5sa.srt,
        group.by = 'Cell_type3',
        CellChatDB = CellChatDB.my,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)
PlotChord <- function(source, target){
        netVisual_chord_gene(drop.ctrl2.cch, sources.use = source, targets.use = target,
                             slot.name = 'netP',
                             scale = T,
                             title.name = 'Control',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.y5sa2.cch, sources.use = source, targets.use = target,
                             slot.name = 'netP',
                             scale = T,
                             title.name = 'YAP5SA',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.ctrl2.cch, sources.use = source, targets.use = target,
                             slot.name = 'net',
                             scale = T,
                             title.name = 'Control',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.y5sa2.cch, sources.use = source, targets.use = target,
                             slot.name = 'net',
                             scale = T,
                             title.name = 'YAP5SA',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
}

## C3+/-CF to CM communication ####
PlotPDF('07.3.cellchat.c3_cf_to_cm1_cm2', 8, 8)
PlotChord(c('C3+CF', 'C3-CF'), c('aCM1', 'aCM2'))
dev.off()

## CM to C3+/-CF communication ####
PlotPDF('07.4.cellchat.cm1_cm2_to_c3_cf', 8, 8)
PlotChord(c('aCM1', 'aCM2'), c('C3+CF', 'C3-CF'))
dev.off()

## C3ar1+/-MP to C3+/-CF communication ####
PlotPDF('07.5.cellchat.c3ar1_mp_to_c3_cf', 8, 8)
PlotChord(c('C3ar1+MP', 'C3ar1-MP'), c('C3+CF', 'C3-CF'))
dev.off()

## C3+/-CF to C3ar1+/-MP communication ####
PlotPDF('07.6.cellchat.c3_cf_to_c3ar1_mp', 8, 8)
PlotChord(c('C3+CF', 'C3-CF'), c('C3ar1+MP', 'C3ar1-MP'))
dev.off()

## CM to C3ar1+/-MP communication ####
PlotPDF('07.7.cellchat.cm1_cm2_to_c3ar1_mp', 8, 8)
PlotChord(c('aCM1', 'aCM2'), c('C3ar1+MP', 'C3ar1-MP'))
dev.off()

netVisual_chord_gene(drop.y5sa2.cch, sources.use = 'C3+CF', targets.use = 'aCM2',
                     slot.name = 'net', signaling = c('VCAM', 'SEMA3', 'VEGF'),
                     scale = T,
                     title.name = 'YAP5SA',
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
plotGeneExpression(drop.y5sa2.cch, features = c('Sema3c', 'Itgb1' ,'Itga9', 'Pgf'), type = 'dot', color.use = 'RdYlBu')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(drop.ctrl2.cch, 'analysis/part91.dropseq_triad_control_cellchat.cch.rds')
saveRDS(drop.y5sa2.cch, 'analysis/part91.dropseq_triad_yap5sa_cellchat.cch.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #8. Reviewer #1 Q4  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For display
hs_ctrl_st_p1.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_control_P1.rds')
hs_ctrl_st_p17.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_control_P17.rds')
hs_mi_st_p2.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_IZ_BZ_P2.rds')
hs_mi_st_p3.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_RZ_BZ_P3.rds')
hs_mi_st_p6.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_RZ_P6.rds')
hs_mi_st_p10.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_IZ_P10.rds')

srt.list <- list(hs_ctrl_st_p1.srt, hs_ctrl_st_p17.srt, hs_mi_st_p2.srt,
                 hs_mi_st_p3.srt, hs_mi_st_p6.srt, hs_mi_st_p10.srt)

hs_merge.srt <- merge(srt.list[[1]], srt.list[2:6], merge.dr = 'spatial')
hs_merge.srt$Sample <- factor(hs_merge.srt$donor_id, levels = c('P1', 'P17', 'P6', 'P3', 'P2', 'P10'))
hs_merge.srt$Name <- factor('Control', levels = c('Control', 'MI RZ', 'MI RZ/BZ', 'MI BZ/IZ'))
hs_merge.srt$Name[hs_merge.srt$Sample == 'P6'] <- 'MI RZ'
hs_merge.srt$Name[hs_merge.srt$Sample == 'P3'] <- 'MI RZ/BZ'
hs_merge.srt$Name[hs_merge.srt$Sample == 'P2'] <- 'MI BZ/IZ'
hs_merge.srt$Name[hs_merge.srt$Sample == 'P10'] <- 'MI BZ/IZ'
hs_merge.srt$Sample <- revalue(hs_merge.srt$Sample, replace = c(
        'P1' = 'Control P1',
        'P17' = 'Control P17',
        'P6' = 'MI Remote P6',
        'P3' = 'MI Remote/Border P3',
        'P2' = 'MI Ischemic/Border P2',
        'P10' = 'MI Ischemic P10'
))
Table(hs_merge.srt$Name, hs_merge.srt$Sample)

zone_mk <- readRDS('~/Documents/Bioinformatics/r/db/adult_mouse_post_mi_zone_marker_lisftover.list.rds')
zone_mk2 <- list()
for( i in 1:4){zone_mk2[[i]] <- ConvertGeneID(zone_mk[[i]], species = 'human', from = 'symbol', to = 'id')$GENEID}
hs_merge.srt <- AddModuleScore2(hs_merge.srt, features = zone_mk2, return_z = T, names = names(zone_mk))
n <- 6
p <- FeaturePlotST(hs_merge.srt, features = c('RZ', 'BZ2', 'IZ'),
                   minvals = rep(0, n), maxvals = rep(2, n),
                   pt.sizes = c(0.42, 0.75, 0.6, 0.5, 0.5, 0.5), ncol = n) &
        labs(fill = 'Z-score') &
        theme(aspect.ratio = 1)
p[[6]] <- p[[6]] + RestoreLegend()
p[[12]] <- p[[12]] + RestoreLegend()
p[[18]] <- p[[18]] + RestoreLegend()
PlotPDF('08.1.feat_st.human_mi_zone_score', 17, 10)
p
dev.off()

## Annotate Zones
hs_merge.srt$Zone <- factor('Contol/RZ', levels = c('Contol/RZ', 'BZ', 'IZ'))
hs_merge.srt$Zone[hs_merge.srt$BZ2 > 0.5] <- 'BZ'
hs_merge.srt$Zone[hs_merge.srt$IZ > 1] <- 'IZ'
p <- DimPlotST(hs_merge.srt, group.by = 'Zone', pt.sizes = c(0.42, 0.75, 0.6, 0.5, 0.5, 0.5),
          ncol = n, legend = n, cols = Color_Zone) &
        theme(aspect.ratio = 1)
p
PlotPDF('08.2.dim_st.human_mi_zone_annotation', 17, 5)
p
dev.off()

## Score aCM1 aCM2 signature:
cm1 <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$CM1, from = 'mouse', to = 'human')
cm2 <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$CM2, from = 'mouse', to = 'human')
cm1 <- ConvertGeneID(cm1, species = 'human', from = 'symbol', to = 'id')$GENEID
cm2 <- ConvertGeneID(cm2, species = 'human', from = 'symbol', to = 'id')$GENEID
cm_mk <- list('CM1_mk' = cm1, 'CM2_mk' = cm2)
saveRDS(cm_mk, 'analysis/part91.human_cm1_cm2_markers.list.rds')

hs_merge.srt <- AddModuleScore2(hs_merge.srt, features = cm_mk, names = c('Score_aCM1', 'Score_aCM2'), return_z = T)
n <- 6
p <- FeaturePlotST(hs_merge.srt, features = c('Score_aCM1', 'Score_aCM2'),
                   minvals = rep(-1.5, n), maxvals = rep(2, n),
                   pt.sizes = c(0.42, 0.75, 0.6, 0.5, 0.5, 0.5), ncol = n) &
        theme(aspect.ratio = 1)
p[[6]] <- p[[6]] + RestoreLegend()
p[[12]] <- p[[12]] + RestoreLegend()
PlotPDF('08.3.feat_st.human_mi_cm1_cm2_score', 17, 6)
p
dev.off()
saveRDS(hs_merge.srt@meta.data, 'analysis/part91.Kramann_cm1_cm2_score_for_display.srt_meta.rds')


## For scoring
dir <- '/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/'
srt.list <- list(
        readRDS(paste0(dir, 'Visium_control_P1.rds')),
        readRDS(paste0(dir, 'Visium_control_P7.rds')),
        readRDS(paste0(dir, 'Visium_control_P8.rds')),
        readRDS(paste0(dir, 'Visium_control_P17.rds')),
        readRDS(paste0(dir, 'Visium_IZ_BZ_P2.rds')),
        readRDS(paste0(dir, 'Visium_IZ_P15.rds')),
        readRDS(paste0(dir, 'Visium_IZ_P16.rds')),
        readRDS(paste0(dir, 'Visium_IZ_P10.rds')),
        readRDS(paste0(dir, 'Visium_GT_IZ_P13.rds')),
        readRDS(paste0(dir, 'Visium_GT_IZ_P15.rds')),
        readRDS(paste0(dir, 'Visium_GT_IZ_P9_rep2.rds')),
        readRDS(paste0(dir, 'Visium_IZ_P3.rds')),
        readRDS(paste0(dir, 'Visium_RZ_P6.rds')),
        readRDS(paste0(dir, 'Visium_RZ_BZ_P12.rds')),
        readRDS(paste0(dir, 'Visium_RZ_BZ_P2.rds')),
        readRDS(paste0(dir, 'Visium_RZ_P3.rds')),
        readRDS(paste0(dir, 'Visium_RZ_P11.rds')),
        readRDS(paste0(dir, 'Visium_RZ_P9.rds'))
)
label <- c('ctrl_p1', 'ctrl_p7', 'ctrl_p8', 'ctrl_p17',
           'ibz_p2', 'ibz_p15', 'ibz_p16', 'ibz_p10', 'ibz_p13', 'ibz_p15', 'ibz_p9', 'ibz_p3',
           'rbz_p6', 'rbz_p12', 'rbz_p2', 'rbz_p3', 'rbz_p11', 'rbz_p9')
for(i in 1:L(srt.list)){
        srt.list[[i]]$Sample <- label[i]
        srt.list[[i]]$Group <- str_split(label[i], pattern = '_', simplify = T)[, 1]
}
hs_merge.srt2 <- merge(srt.list[[1]], srt.list[2:L(srt.list)], merge.dr = 'spatial')
hs_merge.srt2$Group <- factor(hs_merge.srt2$Group, levels = c('ctrl', 'rbz', 'ibz'))
hs_merge.srt2$Group <- revalue(hs_merge.srt2$Group,
                               replace = c('ctrl' = 'Control', 'rbz' = 'Remote/Border', 'ibz' = 'Ischemic/Border'))
Table(hs_merge.srt2$Sample, hs_merge.srt2$Group)
hs_merge.srt2 <- AddModuleScore2(hs_merge.srt2, features = cm_mk, names = c('Score_aCM1', 'Score_aCM2'), return_z = T)
data <- hs_merge.srt2@meta.data[, c('Group', 'Sample', 'Score_aCM1', 'Score_aCM2')] |>
        group_by(Sample) |>
        mutate(Score_aCM1 = mean(Score_aCM1), Score_aCM2 = mean(Score_aCM2))
colnames(data)[3:4] <- c('aCM1', 'aCM2')
data <- data[! duplicated(data$Sample),]
data <- reshape2::melt(data)
p <- ggplot(data, aes(x = Group, y = value, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(size = 1.5, color = 'black', ) +
        labs(title = 'aCM1 and aCM2 Signature Scores', y = 'Z-score', x = '', fill = 'Human MI Sample',
             caption = 'aCM1 ANOVA p < 0.001\naCM2 ANOVA p < 0.024') +
        scale_fill_manual(values = Color_Zone) +
        scale_y_continuous(limits = c(-1.8, 1.8)) +
        theme_classic() +
        theme(aspect.ratio = 3, axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~variable)
pval1 <- aov(data = data[data$variable == 'aCM1', ], formula = value~Group)
pval2 <- aov(data = data[data$variable == 'aCM2', ], formula = value~Group)
p
PlotPDF('08.4.box.human_mi_cm1_cm2_score_group_by_individual', 5, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #9. C3+ CF Percentage  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drop_cf.srt <- drop.srt[, drop.srt$Cell_type2 == 'CF' & !is.na(drop.srt$Sample)]
Table(drop_cf.srt$Sample)

DimPlot2(drop_cf.srt)
FeaturePlot2(drop_cf.srt, features = 'C3')
VlnPlot(drop_cf.srt, features = 'C3', assay = 'RNA')

x <- as.matrix(table(drop_cf.srt@assays$RNA@data['C3', ] > 0.2, drop_cf.srt$Sample))
total <- colSums(x)
x[1,] <- x[1, ]/total
x[2,] <- x[2, ]/total
rownames(x) <- c('C3-CF', 'C3+CF')
data <- as.data.frame(x)
data$Group <- str_split(data$Var2, pattern = ' ', simplify = T)[, 1]
p <- ggplot(data[data$Var1 == 'C3+CF', ], aes(y = Freq, x = Group, color = Group)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'C3+ CF', y = 'Fraction of CFs', x = '') +
        scale_y_continuous(limits = c(0, 0.5)) +
        theme_classic() +
        theme(aspect.ratio = 2) +
        RotatedAxis() +
        NoLegend()
PlotPDF('09.1.box.c3_cf_composition_across_samples', 2, 4)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


