####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- '92_NCR_Revision_2'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2022_yap5sa_st_rli/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/stRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject(Machine = 'Gondor', Ver = Ver, Part = Step, Catagory = 'mouse',
                Project_dir = '2022_yap5sa_st_rli', Data_drive = 'bree')
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

sc.srt <- readRDS('external/fmeng_aav9_yap5sa_cd45/11.cd45_annotated.simple.srt.rds')

Color_cell_type <- mycol_10[3:8]
Color_Zone <- mycol_20[c(4, 9, 13)]

olson_cm.srt <- readRDS('individual/part01.2020_DevCell_EOlson.srt.rds')
olson_cm.srt$Sample_pub <- revalue(olson_cm.srt$Sample, replace = c(
        'P1_1Sham' = 'P1 Sham - 1 day',
        'P1_1MI' = 'P1 MI - 1 day',
        'P1_3Sham' = 'P1 Sham - 3 day',
        'P1_3MI' = 'P1 MI - 3 day',
        'P8_1Sham' = 'P8 Sham - 1 day',
        'P8_1MI' = 'P8 MI - 1 day',
        'P8_3Sham' = 'P8 Sham - 3 day',
        'P8_3MI'  = 'P8 MI - 3 day'
))
olson_cm.srt$Sample_pub <- factor(olson_cm.srt$Sample_pub, levels = c(
        'P1 Sham - 1 day',
        'P1 MI - 1 day',
        'P1 Sham - 3 day',
        'P1 MI - 3 day',
        'P8 Sham - 1 day',
        'P8 MI - 1 day',
        'P8 Sham - 3 day',
        'P8 MI - 3 day'
))
olson_cm.srt$CM2_like <- factor('Other CMs', levels = c('Yap5sa CM2-like', 'Other CMs'))
olson_cm.srt$CM2_like[olson_cm.srt$Zscore_y5sa_CM2 >= 1.5] <- 'Yap5sa CM2-like'
olson_cm.srt$CM2_like2 <- factor(olson_cm.srt$CM2_like, levels = c('Other CMs', 'Yap5sa CM2-like'))
DefaultAssay(olson_cm.srt) <- 'RNA'
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

st_wh_y5sa.srt <- RunPCA(st_wh.srt[, st_wh.srt$Sample != 'Control'], assay = 'SCT')
st_wh_y5sa.srt <- RunUMAP(st_wh_y5sa.srt, dims = 1:50, reduction = 'pca', min.dist = 0.2)
st_wh_y5sa.srt <- FindNeighbors(st_wh_y5sa.srt, reduction = 'pca', dims = 1:50)
st_wh_y5sa.srt <- FindClusters(st_wh_y5sa.srt, resolution = seq(0.1, 1, 0.1))
st_wh_y5sa.srt@images$Control <- NULL

st_wh_ctrl.srt$Niche_ctrl <- factor(as.numeric(st_wh_ctrl.srt$SCT_snn_res.0.5), levels = 1:7)
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('5' = '0'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('7' = '5'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('4' = '7'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('6' = '4'))
st_wh_ctrl.srt$Niche_ctrl <- revalue(st_wh_ctrl.srt$Niche_ctrl, replace = c('0' = '6'))
st_wh_ctrl.srt$Niche_ctrl <- factor(st_wh_ctrl.srt$Niche_ctrl, levels = 1:7)

##~~~~ Combined both sections ####
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
PlotPDF('01.1.st_final_14_niches', 9, 4)
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
p2
PlotPDF('01.2.st_final_14_niches_umap', 9, 4)
p2
dev.off()

##~~~~ Cell Type Abundance ####
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
        labs(fill = 'Relative\nAbundance', y = 'Cell Type', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p1
PlotPDF('01.3.st_cell_type_relative_abundance_each_niche', 5, 3)
p1
dev.off()

##~~~~ Triad Abundance ####
mtx <- st_wh.srt@meta.data[, grepl('^Zscore|^Niche_split', colnames(st_wh.srt@meta.data))]
dim(mtx)
mtx <- mtx[! mtx$Niche_split %in% c('Niche C6', 'Niche C7', 'Niche Y6', 'Niche Y7'), ] ## exclude atria
Table(mtx$Niche_split)
mtx2 <- matrix(NA, 10, 4)
colnames(mtx2) <- str_split(colnames(mtx[1:4]), pattern = '_', simplify = T)[, 2]
rownames(mtx2) <- c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))
for(i in 1:10){ for(j in 1:4){
        mtx2[i, j] <- mean(mtx[mtx$Niche_split == c(paste0('Niche C', 1:5), paste0('Niche Y', 1:5))[i], j], ) }
}
data <- reshape2::melt(mtx2)
data$Var2 <- factor(data$Var2, levels = c('CM1', 'CM2', 'C3ar1', 'C3'))
data$Var2 <- revalue(data$Var2, replace = c('CM1'='CM1', 'CM2'='CM2',
                                            'C3ar1'='C3ar1+MP', 'C3'='C3+CF'))
p3 <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c(limits = c(-1, 1.6)) +
        theme_classic() +
        theme(aspect.ratio = 4/10) +
        labs(fill = 'Marker\nExpression\nZ-score', y = 'Cell State Marker', x = 'Niche') +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p3
PlotPDF('01.4.st_traid_marker_expr_each_niche', 5, 3)
p3
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(st_wh_ctrl.srt, 'analysis/part92.ST_YAP5SA_wholeheart_ctrl_niche_annotated.seurat.rds')
saveRDS(st_wh_y5sa.srt, 'analysis/part92.ST_YAP5SA_wholeheart_y5sa_niche_annotated.seurat.rds')
saveRDS(st_wh.srt, 'analysis/part92.ST_YAP5SA_wholeheart_niche_annotated.seurat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #2. Reviewer #1 Q2  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1.1 <- DimPlot2(drop.srt, group.by = 'Cell_type', reduction = 'umap',
                 cols = Color_cell_type, label = T, raster = T) +
        labs(title = '', x = 'UMAP 1', y = 'UMAP 2')
p1.2 <- CountCellBarPlot(drop.srt[, !is.na(drop.srt$Sample)],
                         group.var = 'Sample', stack.var = 'Cell_type',
                         percentage = T, stack.color = Color_cell_type, width = 0.8)$plot +
        labs(x = 'Samples', y = 'Fraction of All Cells', fill = 'Cell Type') +
        theme(aspect.ratio = 1)

p1.1 | p1.2
PlotPDF('02.1.global_umap_per_sample_dist', 8, 4)
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
p1.3 <- ggplot(data1, aes(y = Frac, x = Genotype, fill = Genotype)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'aCM1 Composition', y = 'Fraction of CMs', x = '')
p1.4 <- ggplot(data2, aes(y = Frac, x = Genotype, fill = Genotype)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'aCM2 Composition', y = '', x = '')
p1.5 <- p1.3 + p1.4 &
        scale_y_continuous(limits = c(0.2, 0.8)) &
        theme_classic() &
        theme(aspect.ratio = 3) &
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

p1.6 <- ggplot(df, aes(x = Genotype, y = Fraction, fill = Genotype)) +
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

p_val <- data.frame(
        Compare = c('aCM1', 'aCM2', 'CM', 'CF', 'MP', 'EC1', 'EC2', 'Mural'),
        PVal = c(
                wilcox.test(Frac~Genotype, data = data1)$p.value,
                wilcox.test(Frac~Genotype, data = data2)$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'CM', ])$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'CF', ])$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'MP', ])$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'EC1', ])$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'EC2', ])$p.value,
                wilcox.test(Fraction~Genotype,  data = df[df$GroupVar == 'Mural', ])$p.value
        ))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #3. Reviewer #1 Q3  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
data$value[data$value > 0.15] <- 0.15
data$Var1 <- revalue(data$Var1, replace = c('CM1' = 'aCM1', 'CM2' = 'aCM2', 'SMC' = 'Mural', 'Mac' = 'MP'))
data$Var2 <- revalue(data$Var2, replace = c('CM1' = 'aCM1', 'CM2' = 'aCM2', 'SMC' = 'Mural', 'Mac' = 'MP'))
data$Var1 <- factor(data$Var1, levels = c('aCM1', 'aCM2', 'CF', 'EC', 'EndoC', 'Mural', 'MP'))
data$Var2 <- factor(data$Var2, levels = c('aCM1', 'aCM2', 'CF', 'EC', 'EndoC', 'Mural', 'MP'))
p <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c(na.value = 0) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        facet_wrap(~Sample)
p
PlotPDF('03.1.heatamp.st_colocal_pairwise', 10, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(data, 'analysis/part92.cell_type_colocalization_matrix.dataframe.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #4. Reviewer #1 Q6  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####~~~~ Calculate Y5SA dataset CM1/2-C3-C3ar1 scores ####
st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = T,
                          names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)))

st.srt <- AddModuleScore2(st.srt, features = Yap_target[1], assay = 'SCT', return_z = T,
                          names = c('Zscore_Y5SA_Target'))

st.srt <- AddModuleScore2(st.srt, features = list('C3', 'C3ar1'), assay = 'SCT', return_z = T,
                          names = paste0('Zscore_', c('C3', 'C3ar1')))

RidgePlot(st.srt, features = c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3', 'Zscore_C3ar1',
                               'Zscore_Y5SA_Target', 'Decon_CM2'))

st.srt$Decon_CM2_for_colocal <- st.srt$Decon_CM2
st.srt$Decon_CM2_for_colocal[st.srt$Decon_CM2 == 0] <- NA

####~~~~ Calculate Y5SA dataset CM1/2-C3-C3ar1 colocalization ####
st.srt$Colocal_CM2_C3ar1_C3 <- GetColocalProb(st.srt, meta_features = c('Zscore_CM2',
                                                                        'Zscore_C3',
                                                                        'Zscore_C3ar1'))
st.srt$Colocal_CM2_C3ar1_C3[is.na(st.srt$Colocal_CM2_C3ar1_C3)] <- 0
st.srt$Colocal_CM2_C3 <- GetColocalProb(st.srt, meta_features = c('Zscore_CM2',
                                                                  'Zscore_C3'))
st.srt$Colocal_CM2_C3[is.na(st.srt$Colocal_CM2_C3)] <- 0

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

####~~~~ Quantify Y5SA dataset MC2-C3-C3ar1 colocalization (Bar) ####
st.srt$Psig_CM2_C3 <- factor('p<0.05', levels = c('Not sig.', 'p<0.05'))
st.srt$Psig_CM2_C3_C3ar1 <- factor('p<0.05', levels = c('Not sig.', 'p<0.05'))
st.srt$Psig_CM2_C3[st.srt$Colocal_CM2_C3 < -log10(0.05)] <- 'Not sig.'
st.srt$Psig_CM2_C3_C3ar1[st.srt$Colocal_CM2_C3ar1_C3 < -log10(0.05)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st.srt$Psig_CM2_C3, st.srt$Sample)),
            as.data.frame(Table(st.srt$Psig_CM2_C3_C3ar1, st.srt$Sample)))
df$Group <- rep(c('CM2 + C3 Hi', 'CM2 + C3ar1 Hi + C3 Hi'), each = 2*LU(st.srt$Sample))
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

####~~~~ Quantify Y5SA dataset MC2-C3-C3ar1 colocalization (Violin) ####
st.srt$Prob_CM2_C3 <- 1-(10^-(st.srt$Colocal_CM2_C3))
st.srt$Prob_CM2_C3_C3ar1 <- 1-(10^-(st.srt$Colocal_CM2_C3ar1_C3))
p <- VlnPlot2(st.srt, features = c('Prob_CM2_C3', 'Prob_CM2_C3_C3ar1'), group.by = 'Sample', ncol = 1, pt.size = 0.1) &
        theme_classic() &
        theme(aspect.ratio = 1) &
        RotatedAxis() &
        labs(y = 'Probability', x = '', fill = '') &
        NoLegend()
p
PlotPDF('04.3.vln.cm2_c3_c3ar1_colocolization_quantification', 3, 6)
print(p)
dev.off()

meta <- st.srt@meta.data[, c('Zscore_CM1', 'Zscore_CM2', 'Zscore_C3ar1_Mac', 'Zscore_C3_FB',
                             'Zscore_Y5SA_Target', 'Zscore_C3', 'Zscore_C3ar1',
                             'Colocal_CM2_C3ar1_C3', 'Colocal_CM2_C3',
                             'Psig_CM2_C3', 'Psig_CM2_C3_C3ar1', 'Prob_CM2_C3', 'Prob_CM2_C3_C3ar1')]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(meta, 'analysis/part92.ST_YAP5SA_venticle_colocal_meta.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~ Calculate Olson dataset CM1/2-C3CF-C3ar1MP scores ####
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

####~~~~ Calculate Olson dataset CM1/2-C3CF-C3ar1MP colocalization ####
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

####~~~~ Quantify Olson dataset CM1/2-C3CF-C3ar1MP colocalization ####
st2.srt$Psig_CM1_C3CF_C3ar1MP <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM2_C3CF_C3ar1MP <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM1_C3CF_C3ar1MP[st2.srt$Colocal_CM1_C3ar1Mac_C3FB < -log10(0.01)] <- 'Not sig.'
st2.srt$Psig_CM2_C3CF_C3ar1MP[st2.srt$Colocal_CM2_C3ar1Mac_C3FB < -log10(0.01)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1_C3CF_C3ar1MP, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2_C3CF_C3ar1MP, st2.srt$Sample)))
df$Group <- rep(c('CM1, C3ar1+ Mac, C3+ FB', 'CM2, C3ar1+ Mac, C3+ FB'), each = 2*LU(st2.srt$Sample))
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

####~~~~ Calculate Olson dataset CM1/2-YapTarget colocalization ####
st2.srt$Colocal_CM1_YapTarget <- GetColocalProb(st2.srt, meta_features = c('Zscore_CM1', 'Zscore_Target'))
st2.srt$Colocal_CM2_YapTarget <- GetColocalProb(st2.srt, meta_features = c('Zscore_CM2', 'Zscore_Target'))

minvals <- rep(0, 4)
maxvals <- rep(5, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_YapTarget',
                                'Colocal_CM2_YapTarget'),
                   title = c('CM1, Yap Hi',
                             'CM2, Yap Hi'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.6,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 5), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('04.5.feat_st.olson_st_yap5sa_marker_yap_target_colocolization', 12, 5)
print(p)
dev.off()

####~~~~ Quantify Olson dataset CM1/2-YapTarget colocalization ####
st2.srt$Psig_CM1_YapHi <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM2_YapHi <- factor('p<0.01', levels = c('Not sig.', 'p<0.01'))
st2.srt$Psig_CM1_YapHi[st2.srt$Colocal_CM1_YapTarget < -log10(0.01)] <- 'Not sig.'
st2.srt$Psig_CM2_YapHi[st2.srt$Colocal_CM2_YapTarget < -log10(0.01)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1_YapHi, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2_YapHi, st2.srt$Sample)))
df$Group <- rep(c('CM1 + Yap Hi', 'CM2 + Yap Hi'), each = 2*LU(st2.srt$Sample))
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
saveRDS(meta, 'analysis/part92.ST_Olson_colocal_meta.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #5. Reviewer #1 Q8  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~ Slingshot infered pseudotime  ####
# library('slingshot')
# sce <- as.SingleCellExperiment(drop_cm.srt, assay = 'SCT')
# sce <- slingshot(sce, reducedDim = 'HARMONY', clusterLabels = 'CM_State', start.clus = 'CM1', use.median = T)
# # Plot PC1 vs PC2 colored by Slingshot pseudotime.
# colors <- rainbow(50, alpha = 1)

##~~~~ Find temporal correlated genes using Slignshot Pseudotime ####
# Y <- drop_cm.srt@assays$SCT@data
# var3K <- drop_cm.srt@assays$SCT@var.features[1:3000]
# Y <- Y[var3K, ]  # only counts for variable genes
#
# t <- drop_cm.srt$SlingRank
# pcc.df <- data.frame(PCC = rep(NA, 3000), Pval = rep(NA, 3000), row.names = rownames(Y))
# for(i in 1:3000){
#         pcc.df[i, ] <- c(cor(Y[i,], t),
#                          cor.test(Y[i,], t)$p.value)
# }
# pcc.df <- pcc.df[! is.na(pcc.df$PCC), ]
# hist(pcc.df$PCC, breaks = 100)
# pcc.df2 <- pcc.df[pcc.df$Pval <= 0.01 & abs(pcc.df$PCC) >= 0.2, ]
# module2 <- split(rownames(pcc.df2), pcc.df2$PCC > 0)
# str(module2)
# names(module2) <- c('Neg', 'Pos')
# pcc.df3 <- pcc.df2[pcc.df2$PCC > 0, ]
# colnames(pcc.df3) <- c('Pearson Correlation Coefficient', 'P value')
# pcc.df3$`Gene Symbol` <- rownames(pcc.df3)
# pcc.df3 <- pcc.df3[, c(3, 1, 2)]
# WriteCSV(pcc.df3,  title = 'part91.cm1_cm2_transition_correlated_genes')

sce <- readRDS('analysis/part91.drop_cm_slingshot.sce.rds')
module <- readRDS('analysis/part91.cm1_cm2_transition_modules.list.rds')[[2]]

##~~~~ Plot Slingshot pseudotime vs cell stage.
drop_cm.srt$SlingPseudotime <- sce$slingPseudotime_1
drop_cm.srt$SlingRank <- Range01(rank(sce$slingPseudotime_1))
data <- data.frame(PST = rank(sce$slingPseudotime_1), State = sce$CM_State)
p2 <- FeaturePlot2(drop_cm.srt, features = 'SlingRank', reduction = 'newumap', pt.size = 0.5) +
        labs(color = 'Pseudotime', x = 'UMAP1', y = 'UMAP2', title = 'Infered Pseudotime') +
        scale_color_distiller(palette = 'Spectral')
p3 <- DimPlot2(drop_cm.srt, group.by = 'CM_State', reduction = 'newumap', cols = mycol_10, pt.size = 0.5) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'CM Cell State')
p <- wrap_plots(p3, p2, ncol = 2)
p
PlotPDF('05.1.cm1_cm2_slingshot_pseudotime', 9, 4)
print(p)
dev.off()

##~~~~ Plot temporal correlated genes ####
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

## Plot Slingshot modules
drop_cm.srt <- ScaleData(drop_cm.srt, features = unlist(module), assay = 'SCT')
expr_mat <- drop_cm.srt@assays$SCT@scale.data[all_module_genes, colnames(drop_cm.srt)[order(drop_cm.srt$SlingRank)]]
for(i in 1:nrow(expr_mat)){expr_mat[i, ] <- smooth.spline(expr_mat[i, ], spar = 1.1)$y}

df_all <- mapply(function(x, y){
        expr_mat <- t(apply(expr_mat[x,], 1, Range01)) ## normalize to [0-1]
        df <- FlattenExpr(expr_mat[, seq(1, ncol(drop_cm.srt), 1)], x)
        df$SlingRank <- Range01(df$Cells)
        df$mod_id <- names(module)[y]
        return(df)},
        x = module, y = seq_along(module), SIMPLIFY = F)
df_all <- bind_rows(df_all)
p2 <- ggplot(data = df_all[df_all$mod_id == 'Pos', ]) +
        geom_ribbon(aes(x = SlingRank, ymax = `Upper Quartile`, ymin = `Lower Quartile`),
                    alpha = 0.1, show.legend = F) +
        geom_line(aes(x = SlingRank, y = `Median Expr`),
                  show.legend = T, size = 1, alpha = 1) +
        labs(y = 'Normalized expression', x = 'aCM1 to aCM2 Pseudotime',
             title = 'Temporally-correlated Genes', caption = '62 Genes') +
        scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
        theme_classic() +
        theme(aspect.ratio = 0.7,
              axis.text.x = element_blank(),
              text = element_text(color = "black"),
              line = element_line(color = "black"),
              axis.line = element_line(color = "black"),
              axis.ticks.x.bottom = element_blank()
        )
p2

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
PlotPDF('05.2.trend_beeswarm.cm1_cm2_slingshot_module', 4, 4)
p2/p1
dev.off()
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

# hs_ref_cm.srt <- FindNeighbors(hs_ref_cm.srt, reduction = "pca", dims = 1:50) |>
#         FindClusters(resolution = 0.5)
# hs_ref_cm.srt <- AddModuleScore2(hs_ref_cm.srt, features = list(cm2), names = 'Score_aCM2', return_z = T)
# saveRDS(hs_ref_cm.srt@meta.data, 'analysis/part91.Teichmann_vCM_clustered.srt_meta.rds')
hs_ref_cm.srt <- AddMetaData(hs_ref_cm.srt, readRDS('analysis/part91.Teichmann_vCM_clustered.srt_meta.rds'))

## Get CM8 Markers
# Idents(hs_ref_cm.srt) <- 'seurat_clusters'
# deg <- FindMarkers(hs_ref_cm.srt,
#                    ident.1 = '8',
#                    min.pct = 0.25,
#                    random.seed = 123,
#                    logfc.threshold = 0.1)
# deg <- deg[order(deg$avg_log2FC, decreasing = T), ]
# deg$gene <- rownames(deg)
# saveRDS(deg, 'analysis/part91.Teichmann_CM8_marker.srt_marker.rds')
deg <- readRDS('analysis/part91.Teichmann_CM8_marker.srt_marker.rds')

deg$DEG_FOR_GO <- NA
deg$DEG_FOR_GO[deg$avg_log2FC >  0.6 & deg$p_val < 0.05] <- 'CM8_Up'
deg$DEG_FOR_GO[deg$avg_log2FC < -0.6 & deg$p_val < 0.05] <- 'CM8_Down'
WriteCSV(deg, 'part92.teichmann_cm8_vs_other_vcm_deg')

Up <- split(deg$gene, deg$DEG_FOR_GO)[['CM8_Up']]
Dn <- split(deg$gene, deg$DEG_FOR_GO)[['CM8_Down']]

terms <- ModuleEnrichment(module_list = list(Up = Up, Dn = Dn), human_or_mouse = 'human')
go_up <- terms$GO$GO_Up
go_dn <- terms$GO$GO_Dn

WriteCSV(go_up, 'part92.teichmann_cm8_up_enriched_go_bp')
WriteCSV(go_dn, 'part92.teichmann_cm8_dn_enriched_go_bp')
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

drop.ctrl.cch <- DoCellChat(
        drop.ctrl.srt,
        group.by = 'Cell_type3',
        CellChatDB = CellChatDB.my,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)
drop.y5sa.cch <- DoCellChat(
        drop.y5sa.srt,
        group.by = 'Cell_type3',
        CellChatDB = CellChatDB.my,
        LR.type = 'all',
        species = 'mouse',
        trim = 0.1
)

PlotChord <- function(source, target){
        netVisual_chord_gene(drop.ctrl.cch, sources.use = source, targets.use = target,
                             slot.name = 'netP',
                             scale = T,
                             title.name = 'Control',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.y5sa.cch, sources.use = source, targets.use = target,
                             slot.name = 'netP',
                             scale = T,
                             title.name = 'YAP5SA',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.ctrl.cch, sources.use = source, targets.use = target,
                             slot.name = 'net',
                             scale = T,
                             title.name = 'Control',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
        netVisual_chord_gene(drop.y5sa.cch, sources.use = source, targets.use = target,
                             slot.name = 'net',
                             scale = T,
                             title.name = 'YAP5SA',
                             link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
}

## C3+/-CF to CM communication ####
PlotPDF('07.1.cellchat.c3_cf_to_cm1_cm2', 8, 8)
PlotChord(c('C3+CF', 'C3-CF'), c('aCM1', 'aCM2'))
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(drop.ctrl.cch, 'analysis/part92.dropseq_triad_control_cellchat.cch.rds')
saveRDS(drop.y5sa.cch, 'analysis/part92.dropseq_triad_yap5sa_cellchat.cch.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #8. Reviewer #1 Q4  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For display
# hs_ctrl_st_p1.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_control_P1.rds')
# hs_ctrl_st_p17.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_control_P17.rds')
# hs_mi_st_p2.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_IZ_BZ_P2.rds')
# hs_mi_st_p3.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_RZ_BZ_P3.rds')
# hs_mi_st_p6.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_RZ_P6.rds')
# hs_mi_st_p10.srt <- readRDS('/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/Visium_IZ_P10.rds')

hs_ctrl_st_p1.srt <- readRDS('~/Downloads/Visium_control_P1.rds')
hs_ctrl_st_p17.srt <- readRDS('~/Downloads/Visium_control_P17.rds')
hs_mi_st_p2.srt <- readRDS('~/Downloads/Visium_IZ_BZ_P2.rds')
hs_mi_st_p3.srt <- readRDS('~/Downloads/Visium_RZ_BZ_P3.rds')
hs_mi_st_p6.srt <- readRDS('~/Downloads/Visium_RZ_P6.rds')
hs_mi_st_p10.srt <- readRDS('~/Downloads/Visium_IZ_P10.rds')

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

## Score aCM1 aCM2 C3CF and C3ar1MP signature:
cm1 <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$CM1, from = 'mouse', to = 'human')
cm2 <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$CM2, from = 'mouse', to = 'human')
c3cf <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$C3_FB, from = 'mouse', to = 'human')
c3ar1mp <- ConvertGeneSpecies(Yap5sa_CM_Mac_CF_marker$C3ar1_Mac, from = 'mouse', to = 'human')
cm1 <- ConvertGeneID(cm1, species = 'human', from = 'symbol', to = 'id')$GENEID
cm2 <- ConvertGeneID(cm2, species = 'human', from = 'symbol', to = 'id')$GENEID
c3cf <- ConvertGeneID(c3cf, species = 'human', from = 'symbol', to = 'id')$GENEID
c3ar1mp <- ConvertGeneID(c3ar1mp, species = 'human', from = 'symbol', to = 'id')$GENEID

cm_mk <- list('CM1_mk' = cm1, 'CM2_mk' = cm2, 'C3CF_mk' = c3cf, 'C3AR1MP' = c3ar1mp)
saveRDS(cm_mk, 'analysis/part92.human_cm1_cm2_markers.list.rds')

hs_merge.srt <- AddModuleScore2(hs_merge.srt, features = cm_mk,
                                names = c('Score_aCM1', 'Score_aCM2', 'Score_C3CF', 'Score_C3AR1MP'), return_z = T)
n <- 6
p <- FeaturePlotST(hs_merge.srt, features = c('Score_aCM1', 'Score_aCM2', 'Score_C3CF', 'Score_C3AR1MP'),
                   minvals = rep(-1.5, n), maxvals = rep(2, n),
                   pt.sizes = c(0.42, 0.75, 0.6, 0.5, 0.5, 0.5), ncol = n) &
        theme(aspect.ratio = 1)
p[[6]] <- p[[6]] + RestoreLegend()
p[[12]] <- p[[12]] + RestoreLegend()
p[[18]] <- p[[18]] + RestoreLegend()
p[[24]] <- p[[24]] + RestoreLegend()
PlotPDF('08.3.feat_st.human_mi_cm1_cm2_score', 17, 12)
p
dev.off()
saveRDS(hs_merge.srt@meta.data, 'analysis/part92.Kramann_cm1_cm2_score_for_display.srt_meta.rds')


## For scoring
# dir <- '/Volumes/shire/data/visium/2022_Nature_RKramann/matrix_public/'
dir <- '~/Downloads/'
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
hs_merge.srt2 <- AddModuleScore2(hs_merge.srt2, features = cm_mk,
                                 names = c('Score_aCM1', 'Score_aCM2', 'Score_C3CF', 'Score_C3AR1MP'), return_z = T)
data <- hs_merge.srt2@meta.data[, c('Group', 'Sample', 'Score_aCM1', 'Score_aCM2', 'Score_C3CF', 'Score_C3AR1MP')] |>
        group_by(Sample) |>
        mutate(Score_aCM1 = mean(Score_aCM1),
               Score_aCM2 = mean(Score_aCM2),
               Score_C3CF = mean(Score_C3CF),
               Score_C3AR1MP = mean(Score_C3AR1MP))
colnames(data)[3:6] <- c('aCM1', 'aCM2', 'C3CF', 'C3AR1MP')
data <- data[! duplicated(data$Sample),]
data <- reshape2::melt(data)
p <- ggplot(data, aes(x = Group, y = value, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(size = 1.5, color = 'black', ) +
        labs(title = 'aCM1 and aCM2 Signature Scores', y = 'Z-score', x = '', fill = 'Human MI Sample',
             caption = 'aCM1 ANOVA p < 1e-200.001\naCM2 ANOVA p < 0.024') +
        scale_fill_manual(values = Color_Zone) +
        scale_y_continuous(limits = c(-1.8, 1.8)) +
        theme_classic() +
        theme(aspect.ratio = 3, axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~variable)
pval1 <- aov(data = data[data$variable == 'aCM1', ], formula = value~Group)
summary(pval1)
pval2 <- aov(data = data[data$variable == 'aCM2', ], formula = value~Group)
summary(pval2)
p
PlotPDF('08.4.box.human_mi_cm1_cm2_score_group_by_individual', 5, 10)
p
dev.off()

hs_merge.srt2 <- AddModuleScore2(hs_merge.srt2, features = zone_mk2, return_z = T, names = names(zone_mk))
hs_merge.srt2$Zone <- factor('Contol/RZ', levels = c('Contol/RZ', 'BZ', 'IZ'))
hs_merge.srt2$Zone[hs_merge.srt2$BZ2 > 0.5] <- 'BZ'
hs_merge.srt2$Zone[hs_merge.srt2$IZ > 1] <- 'IZ'
MappingHeatmap(hs_merge.srt2, que_var = 'Zone', ref_var = 'Group', percentage = F)
p1 <- BoxPlot(hs_merge.srt2, feature = 'Score_aCM1', group.by = 'Zone') +
        scale_y_continuous(limits = c(-2.6, 3)) +
        labs(title = 'aCM1 Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ')
p2 <- BoxPlot(hs_merge.srt2, feature = 'Score_aCM2', group.by = 'Zone') +
        scale_y_continuous(limits = c(-2.6, 3)) +
        labs(title = 'aCM2 Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ')
p3 <- BoxPlot(hs_merge.srt2, feature = 'Score_C3CF', group.by = 'Zone') +
        scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3+CF Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ')
p4 <- BoxPlot(hs_merge.srt2, feature = 'Score_C3AR1MP', group.by = 'Zone') +
        scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3AR1+MP Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ')

p <- p1 + p2 + p3 + p4 &
        theme(aspect.ratio = 3) &
        scale_fill_manual(values = Color_Zone)
p
pval1 <- c(t.test(split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[2]])$p.value,
           t.test(split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[3]])$p.value,
           t.test(split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[2]],
                  split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[3]])$p.value)
pval2 <- c(t.test(split(hs_merge.srt2$Score_aCM2, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_aCM2, hs_merge.srt2$Zone)[[2]])$p.value,
           t.test(split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_aCM1, hs_merge.srt2$Zone)[[3]])$p.value,
           t.test(split(hs_merge.srt2$Score_aCM2, hs_merge.srt2$Zone)[[2]],
                  split(hs_merge.srt2$Score_aCM2, hs_merge.srt2$Zone)[[3]])$p.value)
pval3 <- c(t.test(split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[2]])$p.value,
           t.test(split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[3]])$p.value,
           t.test(split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[2]],
                  split(hs_merge.srt2$Score_C3CF, hs_merge.srt2$Zone)[[3]])$p.value)
pval4 <- c(t.test(split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[2]])$p.value,
           t.test(split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[1]],
                  split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[3]])$p.value,
           t.test(split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[2]],
                  split(hs_merge.srt2$Score_C3AR1MP, hs_merge.srt2$Zone)[[3]])$p.value)
all(c(pval1 < 1e-10, pval2 < 1e-10, pval3 < 1e-10, pval4 < 1e-10))
PlotPDF('08.5.box.human_mi_cm1_cm2_score_group_by_spots', 8, 8)
p
dev.off()

hs_merge.srt2$Coloc_CM1_C3CF_C3AR1MP <- GetColocalProb(hs_merge.srt2,
                                                       meta_features = c('Score_aCM1', 'Score_C3CF', 'Score_C3AR1MP')) >
        -log10(0.01)
hs_merge.srt2$Coloc_CM2_C3CF_C3AR1MP <- GetColocalProb(hs_merge.srt2,
                                                       meta_features = c('Score_aCM2', 'Score_C3CF', 'Score_C3AR1MP')) >
        -log10(0.01)
hs_merge.srt2$Coloc_CM1_C3CF <- GetColocalProb(hs_merge.srt2,
                                               meta_features = c('Score_aCM1', 'Score_C3CF')) >
        -log10(0.01)
hs_merge.srt2$Coloc_CM2_C3CF <- GetColocalProb(hs_merge.srt2,
                                               meta_features = c('Score_aCM2', 'Score_C3CF')) >
        -log10(0.01)
CountCellBarPlot(hs_merge.srt2, group.var = 'Zone', stack.var = 'Coloc_CM1_C3CF',
                 stack.color = c('grey90', 'red'))
CountCellBarPlot(hs_merge.srt2, group.var = 'Zone', stack.var = 'Coloc_CM2_C3CF',
                 stack.color = c('grey90', 'red'))
CountCellBarPlot(hs_merge.srt2, group.var = 'Zone', stack.var = 'Coloc_CM1_C3CF_C3AR1MP',
                 stack.color = c('grey90', 'red'))
CountCellBarPlot(hs_merge.srt2, group.var = 'Zone', stack.var = 'Coloc_CM2_C3CF_C3AR1MP',
                 stack.color = c('grey90', 'red'))

BoxPlot(hs_merge.srt2, feature = 'Coloc_CM1_C3CF', group.by = 'Zone') +
        geom_hline(yintercept = -log10(0.05)) +
        #scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3AR1+MP Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ') |
BoxPlot(hs_merge.srt2, feature = 'Coloc_CM2_C3CF', group.by = 'Zone') +
        geom_hline(yintercept = -log10(0.05)) +
        #scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3AR1+MP Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ') |
BoxPlot(hs_merge.srt2, feature = 'Coloc_CM1_C3CF_C3AR1MP', group.by = 'Zone') +
        geom_hline(yintercept = -log10(0.05)) +
        #scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3AR1+MP Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ') |
BoxPlot(hs_merge.srt2, feature = 'Coloc_CM2_C3CF_C3AR1MP', group.by = 'Zone') +
        geom_hline(yintercept = -log10(0.05)) +
        #scale_y_continuous(limits = c(-2, 4)) +
        labs(title = 'C3AR1+MP Signature', y = 'Z-score', x = '', fill = 'Human ST Spots',
             caption = 't test p values < 1e-20 \nfor both BZ vs RZ and BZ vs IZ')


human_mp.srt <- readRDS('/Volumes/shire/data/scrnaseq/2022_Nature_RKramann/matrix_public/Myeloid_snRNA_snATAC.Rds')
DimPlot2(human_mp.srt, group.by = 'annotation', reduction = 'umap_harmony_v2', label = F) /
        FeaturePlot2(human_mp.srt, features = 'C3AR1', reduction = 'umap_harmony_v2', max.cutoff = 'q99')

human_fb.srt <- readRDS('/Volumes/shire/data/scrnaseq/2022_Nature_RKramann/matrix_public/Fibroblast_snRNA_snATAC.Rds')
DimPlot2(human_fb.srt, group.by = 'annotation', reduction = 'umap_harmony_v2', label = F) /
        FeaturePlot2(human_fb.srt, features = 'C3', reduction = 'umap_harmony_v2', max.cutoff = 'q99')
DotPlot2(human_mp.srt, feature = 'C3AR1', group.by = 'region') +
DotPlot2(human_fb.srt, feature = 'C3', group.by = 'region')

hs_mi_st_p3.srt$Sample <- factor('P3')
hs_mi_st_p3.srt <- AddModuleScore2(hs_mi_st_p3.srt, features = cm_mk,
                                   names = c('Score_aCM1', 'Score_aCM2'), return_z = T)
FeaturePlotST(hs_mi_st_p3.srt, features = c('ENSG00000125730', 'ENSG00000171860', 'Score_aCM1', 'Score_aCM2'),
              minvals = c(0, 0, 0, 0), maxvals = c(2, 2, 2, 2), ncol = 2, pt.sizes = 0.2)


hs_merge.srt <- AddModuleScore2(hs_merge.srt, features = list('ENSG00000125730', 'ENSG00000171860'),
                                names = c('C3', 'C3AR1'), return_z = T)
p <- FeaturePlotST(hs_merge.srt, features = c('C3', 'C3AR1'),
                   minvals = rep(-1, n), maxvals = rep(3, n),
                   pt.sizes = c(0.42, 0.75, 0.6, 0.5, 0.5, 0.5), ncol = n) &
        theme(aspect.ratio = 1)
p[[6]] <- p[[6]] + RestoreLegend()
p[[12]] <- p[[12]] + RestoreLegend()
p
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #9. C3+ CF Percentage  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drop_cf.srt <- drop.srt[, drop.srt$Cell_type2 == 'CF' & !is.na(drop.srt$Sample)]
Table(drop_cf.srt$Sample)

x <- as.matrix(table(drop_cf.srt@assays$RNA@data['C3', ] > 0.2, drop_cf.srt$Sample))
total <- colSums(x)
x[1,] <- x[1, ]/total
x[2,] <- x[2, ]/total
rownames(x) <- c('C3-CF', 'C3+CF')
data <- as.data.frame(x)
data$Group <- str_split(data$Var2, pattern = ' ', simplify = T)[, 1]
p <- ggplot(data[data$Var1 == 'C3+CF', ], aes(y = Freq, x = Group, fill = Group)) +
        geom_boxplot() +
        geom_point(size = 1.5, color = 'black') +
        labs(title = 'C3+ CF', y = 'Fraction of CFs', x = '') +
        scale_y_continuous(limits = c(0, 0.5)) +
        theme_classic() +
        theme(aspect.ratio = 2) +
        RotatedAxis() +
        NoLegend()
p
PlotPDF('09.1.box.c3_cf_composition_across_samples', 2, 4)
p
dev.off()

wilcox.test(Freq~Group,  data = data[data$Var1 == 'C3+CF', ])$p.value
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #10. P1MI vs P8MI aCM2 score  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Table(olson_cm.srt@active.ident, olson_cm.srt$CM2_like2)

olson_cm.srt$CM2_like3 <- 'Other CMs'
olson_cm.srt$CM2_like3[olson_cm.srt$CM2_like2 == 'Yap5sa CM2-like' &
                               olson_cm.srt$CM_subtype == 'CM4'] <- 'Immature aCM2-like'
olson_cm.srt$CM2_like3[olson_cm.srt$CM2_like2 == 'Yap5sa CM2-like' &
                               olson_cm.srt$CM_subtype == 'CM5'] <- 'Injury-induced aCM2-like'

olson_cm.srt$CM2_like3 <- factor(olson_cm.srt$CM2_like3,
                                 levels = c('Other CMs', 'Immature aCM2-like', 'Injury-induced aCM2-like'))
olson_cm.srt$Sample_pub2 <- factor(str_split(olson_cm.srt$Sample_pub, ' - ', simplify = T)[, 1],
                                   levels = c('P1 Sham', 'P1 MI', 'P8 Sham', 'P8 MI'))

p3 <- DimPlot2(olson_cm.srt, cols = c('grey85', mycol_10[4:5]), raster = F, pt.size = 0.1,
               label = F, group.by = 'CM2_like3') +
        labs(x = '', y = '', title = 'YAP5SA aCM2-like', caption = '* aCM2 Z-score >= 1', color = '') +
        theme(title = element_text(hjust = 1)) +
        NoLegend()

p5 <- CountCellBarPlot(olson_cm.srt, group.var = 'Sample_pub2', stack.var = 'CM2_like3',
                       percentage = T, stack.color = c('grey85', mycol_10[4:5])) +
        labs(x = '', y = 'Fraction of nuclei', fill = '')
p1 <- DimPlot2(olson_cm.srt, cols = mycol_10, raster = T, pt.size = 1, label = T) +
        NoLegend() +
        labs(x = '', y = '', title = 'CM subtypes -- Cui et al') +
        theme(title = element_text(hjust = 1))
p2 <- FeaturePlot2(olson_cm.srt, features = 'Zscore_y5sa_CM2', pt.size = 1, raster = T,
                   min.cutoff = 'q5', max.cutoff = 'q95', order = F) +
        labs(x = '', y = '', col = 'Z score', title = 'Yap5sa CM2 score')
p4 <- VlnPlot2(olson_cm.srt, features = c('Ctgf', 'Rtn4' ,'Sorbs2', 'Acta2', 'Eef2', 'Ptrf'),
               ncol = 3, cols = mycol_10) &
        labs(y = 'Expression')
p <- wrap_plots(p1 + p2 + p3, p4 | p5, ncol = 1)
p
PlotPDF('10.1.umap_box.cm2_like_cells_in_olson_data', 15, 14)
p
dev.off()

## Label Transfer
drop_cm.srt <- NormalizeData(drop_cm.srt) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
olson_cm.srt <- FindVariableFeatures(olson_cm.srt, selection.method = "vst", nfeatures = 2000)
anchors <- FindTransferAnchors(reference = olson_cm.srt, query = drop_cm.srt, dims = 1:30, query.assay = 'RNA')
predictions <- TransferData(anchorset = anchors, refdata = olson_cm.srt$CM_subtype, dims = 1:30)
drop_cm.srt <- AddMetaData(drop_cm.srt, metadata = predictions)
Table(drop_cm.srt$predicted.id, drop_cm.srt$CM_State)

data <- data.frame(Predicted_identity = paste0('CM', 1:5), Fraction_CM2 = c(0, 0, 0, 1004, 193)/1197)
p <- ggplot(data) +
        geom_col(aes(x = Predicted_identity, y = Fraction_CM2)) +
        theme_Publication() +
        theme(aspect.ratio = 1, panel.grid.major = element_blank())
p
PlotPDF('10.2.bar.cm2_label_transfer_prediction', 4, 4)
p
dev.off()

## Marker heatmap
olson_mk <- FindAllMarkers(olson_cm.srt, only.pos = T)
olson_mk_hi <- olson_mk[olson_mk$p_val_adj <= 0.05 & olson_mk$avg_log2FC > 0.5, ]
Table(olson_mk_hi$cluster)
olson_gl <- split(olson_mk_hi$gene, olson_mk_hi$cluster)
Idents(drop_cm.srt) <- 'CM_State'
y5sa_mk <- FindAllMarkers(drop_cm.srt, only.pos = T)
y5sa_mk_hi <- y5sa_mk[y5sa_mk$p_val_adj <= 0.05 & y5sa_mk$avg_log2FC > 0.5, ]
Table(y5sa_mk_hi$cluster)
y5sa_gl <- split(y5sa_mk_hi$gene, y5sa_mk_hi$cluster)

mk_df <- olson_mk_hi[olson_mk_hi$gene %in% RemoveRiboMito(olson_mk_hi$gene, 'mouse'), ]
Idents(merged.srt) <- 'Cell_state'
levels(merged.srt) <- c(paste0('Neo_CM', 1:5), paste0('Y5SA_CM', 1:2))
sub.srt <- merged.srt[, unlist(DownsampleByMeta(merged.srt[, merged.srt$Cell_state != 'Y5SA_CM1'],
                                                meta_var = 'Cell_state', n = 500))]
sub.srt <- RunALRA(sub.srt, genes.use = U(mk_df$gene))
p <- MarkerHeatmap(sub.srt, marker.df = mk_df, disp.min = 1, top = 100, raster = TRUE, group.cols = mycol_10) +
        scale_fill_viridis_c() +
        theme(aspect.ratio = 1, axis.text.y = element_blank())
PlotPDF('10.3.heat.cm2_express_neo_cm_markers', 8, 8)
p
dev.off()

saveRDS(list(olson_mk_hi, y5sa_mk_hi), 'analysis/part92.olson_y5sa_cm_state.srt_mk.rds')

data <- rbind(colMedians(AverageExpression(sub.srt, features = olson_gl$CM1, group.by = 'Cell_state')$alra),
              colMedians(AverageExpression(sub.srt, features = olson_gl$CM2, group.by = 'Cell_state')$alra),
              colMedians(AverageExpression(sub.srt, features = olson_gl$CM3, group.by = 'Cell_state')$alra),
              colMedians(AverageExpression(sub.srt, features = olson_gl$CM4, group.by = 'Cell_state')$alra),
              colMedians(AverageExpression(sub.srt, features = olson_gl$CM5, group.by = 'Cell_state')$alra))
rownames(data) <- paste0('Neo CM', 1:5, ' Sign.')
colnames(data) <- c(paste0('Neo_CM', 1:5), paste0('Y5SA_CM', 2))
df <- melt(data)
p <- ggplot(df) +
        geom_tile(aes(x = Var2, y = Var1, fill = value), color = 'black') +
        scale_y_discrete(limits = rev) +
        scale_fill_viridis_c() +
        theme_Publication() +
        labs(x = 'Cells', y = 'Genes', fill = 'Median Expression') +
        RotatedAxis() +
        theme(aspect.ratio = 5/6, panel.grid.major = element_blank())
p
PlotPDF('10.3.2.heat.cm2_express_neo_cm_markers', 8, 8)
p
dev.off()

##  PCA
merged.srt <- merge(olson_cm.srt, drop_cm.srt)
merged.srt$Cell_state <- NA
merged.srt$Cell_state[Cells(drop_cm.srt)] <- paste0('Y5SA_', drop_cm.srt$CM_State)
merged.srt$Cell_state[Cells(olson_cm.srt)] <- paste0('Neo_', olson_cm.srt$CM_subtype)
merged.srt$Dataset <- NA
merged.srt$Dataset[Cells(drop_cm.srt)] <- 'Y5SA'
merged.srt$Dataset[Cells(olson_cm.srt)] <- 'Neo'

merged.srt@assays$RNA@var.features <- U(c(olson_mk_hi$gene, y5sa_mk_hi$gene))
merged.srt <- ScaleData(merged.srt) |> RunPCA()
merged.srt <- RunHarmony(merged.srt, group.by.vars = 'Dataset')
df <- as.data.frame(merged.srt@reductions$pca@cell.embeddings[, c(1, 2)])
df$Group[Cells(merged.srt)] <- as.vector(merged.srt$Cell_state)
df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(PC_1), PC_2 = mean(PC_2))
p <- ggplot() +
        #geom_point(data = df, mapping = aes(x = PC_1, y = PC_2, color = Group), alpha = 0) +
        geom_point(data = df2[df2$Group != 'Y5SA_CM1', ], mapping = aes(x = PC_1, y = PC_2, color = Group), size = 3) +
        scale_color_manual(values = mycol_10) +
        #scale_x_continuous(limits = c(-30, 30)) +
        labs(color = '', x = 'PC1', y = 'PC2') +
        theme_classic() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p
PlotPDF('10.4.pc.cm2_similarity_with_neocm4', 4, 4)
p
dev.off()

## Coord mapping
que.srt <- merged.srt[, merged.srt$Cell_state == 'Y5SA_CM2']
# que.srt <- que.srt[, sample(Cells(que.srt), size = 20)]
genes <-  intersect(U(c(olson_mk_hi$gene, y5sa_mk_hi$gene)), rownames(ref.srt))
ref.srt <- olson_cm.srt
ref.srt$CM_subtype <- factor(ref.srt$CM_subtype)
que.srt$Cell_state <- factor(que.srt$Cell_state)
result <- EmbeddingMappingCorr(ref_seurat = ref.srt, ref_assay = 'RNA',
                               que_seurat = que.srt, que_assay = 'RNA',
                               features = genes)
p <- EmbeddingMappingPlot(ref_seurat = ref.srt, ref_cols = mycol_30[seq(3, 30, 3)],
                          ref_group_by = 'CM_subtype', ref_pt_size = 0.1,
                          que_seurat = que.srt, que_cols = 'grey20',
                          que_group_by = 'Cell_state', que_pt_size = 0.1,
                          pt_alpha = 0.5,
                          ref_reduction = 'umap', k = 10, pcc_mat = result, jitter = 0.5)
p
PlotPDF('10.5.umap.cm2_coord_mapping', 8, 10)
p
dev.off()

## Count cells by sample
olson_cm.srt$Sample_pub2 <- factor(str_split(olson_cm.srt$Sample_pub, pattern = ' - ', simplify = T)[, 1],
                                   levels = c('P1 Sham', 'P1 MI', 'P8 Sham', 'P8 MI'))
olson_cm.srt$tmp <- factor(ifelse(olson_cm.srt$CM_subtype == 'CM4', yes = 'aCM2-like Cluster', no = 'Other CM Clusters'),
                           levels = c('Other CM Clusters', 'aCM2-like Cluster'))
p <- CountCellBarPlot(olson_cm.srt, group.var = 'Sample_pub2', stack.var = 'tmp', stack.color = c('grey75', 'red3'))
p
PlotPDF('10.6.bar.cm4_composition', 4, 4)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(result, 'analysis/part92.olson_cm_coord_mapping_result.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  #11. P1MI vs P8MI Non-CMs  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comp.srt <- readRDS('/Volumes/anor/project/2019_scrna_compendium/rdata/mouse_v7.0/integrated/STEP29.pub.final.srt.rds')
data <- comp.srt[, comp.srt$study == '2020_CellRep_EOlson']
Table(data$condition, data$age)
data$Group <- factor('P2MI', levels = c('P2Sham', 'P2MI', 'P8Sham', 'P8MI'))
data$Group[data$condition_group == 'Sham' & data$age %in% c('P002_(0w)', 'P004_(0w)')] <- 'P2Sham'
data$Group[data$condition_group == 'Sham' & data$age %in% c('P009_(1w)', 'P011_(1w)')] <- 'P8Sham'
data$Group[data$condition_group != 'Sham' & data$age %in% c('P002_(0w)', 'P004_(0w)')] <- 'P2MI'
data$Group[data$condition_group != 'Sham' & data$age %in% c('P009_(1w)', 'P011_(1w)')] <- 'P8MI'
Table(data$Group)

data <- RunUMAP(data, dims = 1:30, reduction = 'harmony')
data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony') |> FindClusters(resolution = seq(0.1, 0.5, 0.1))

DimPlot2(data, group.by = 'DCX_snn_res.0.2', split.by = 'Group', cols = mycol_10, ncol = 2, reduction = 'umap')
DimPlot2(data, group.by = 'Cell_type_fine_abv', cols = mycol_30, reduction = 'umap')
Idents(data) <- 'DCX_snn_res.0.2'
mk <- FindAllMarkers(data, logfc.threshold = 1, only.pos = T)

data$Cell_type <- NA
data$Cell_type[data$DCX_snn_res.0.2 %in% c(0)] <- 'EC'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(1)] <- 'FB'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(2)] <- 'MP'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(3)] <- 'Mural'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(4)] <- 'RBC'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(5)] <- 'CM'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(6)] <- 'NPhil'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(7)] <- 'Lym'
data$Cell_type[data$DCX_snn_res.0.2 %in% c(8)] <- 'EpiC'
Idents(data) <- 'Cell_type'
non_cm.srt <- DietSeurat(data, scale.data = F, dimreducs = names(data@reductions)) ## Save this

data$C3_Group <- 'C3_Neg'
data$C3_Group <- ifelse(data@assays$DCX@data['C3', ] > 0, yes = 'C3_Pos', no = 'C3_Neg')
data$C3ar1_Group <- 'C3ar1_Neg'
data$C3ar1_Group <- ifelse(data@assays$DCX@data['C3ar1', ] > 0, yes = 'C3ar1_Pos', no = 'C3ar1_Neg')


p1 <- DimPlot2(data, group.by = 'Cell_type', cols = mycol_10, reduction = 'umap', label = T) +
        labs(x = 'UMAP1', y = 'UMAP2', title = '')
sub.data <- data[, unlist(DownsampleByMeta(data, 'Group', down_to_min_group = T))]
p2 <- FeaturePlot2(sub.data, features = c('C3', 'C3ar1'),
                   split.by = 'Group', ncol = 1, reduction = 'umap', max.cutoff = 'q95')
p3 <- CountCellBarPlot(data[, data$Cell_type %in% c('EpiC', 'FB')], group.var = 'Group', stack.var = 'C3_Group',
                       percentage = T, stack.color = c('grey75', 'red3')) +
        coord_cartesian(ylim = c(0, 0.5)) +
        labs(y = 'Fraction of CFs')
p4 <- CountCellBarPlot(data[, data$Cell_type %in% c('MP', 'NPhil', 'Lym', 'RBC')],
                       group.var = 'Group', stack.var = 'C3ar1_Group',
                         percentage = T, stack.color = c('grey75', 'red3')) +
        coord_cartesian(ylim = c(0, 0.5)) +
        labs(y = 'Fraction of Immune Cells')
p <- (p1 + p3 + p4) / p2
p
PlotPDF('11.1.umap_bar.neo_mi_non_cm', 12, 10)
p
dev.off()

p5 <- CountCellBarPlot(data, group.var = 'Group', stack.var = 'C3_Group',
                       percentage = T, stack.color = c('grey75', 'red3')) +
        coord_cartesian(ylim = c(0, 0.5)) +
        labs(y = 'Fraction of All Cells')
p6 <- CountCellBarPlot(data,
                       group.var = 'Group', stack.var = 'C3ar1_Group',
                       percentage = T, stack.color = c('grey75', 'red3')) +
        coord_cartesian(ylim = c(0, 0.5)) +
        labs(y = 'Fraction of All Cells')
p <- p3 + p4 + p5 + p6
PlotPDF('11.2.bar.c3_c3ar1_norm_to_all_cells', 10, 10)
p
dev.off()

Idents(non_cm.srt) <- 'Cell_type'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(non_cm.srt, 'analysis/part92.olson_non_cm.seurat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#### Fig 1B-C  ####
gl <- c('Sln', 'Nppa', 'Stard10', 'Bmp10', 'Hamp', 'Myl2', 'Ptgds', 'Acta1',
        'Pln', 'Strit1' ,'Lrtm1', 'Hopx', 'Nppb', 'Tcap', 'Synpo2l', 'Gm31659')
p <- SpatialFeaturePlot(st_wh.srt,
            features = gl,
            min.cutoff = 'q5', max.cutoff = 'q95', image.alpha = 0, ncol = 4, stroke = 0)
p1 <- wrap_plots(p[[1]], p[[3]], p[[5]], p[[7]], p[[9]], p[[11]], p[[13]], p[[15]],
                 p[[2]], p[[4]], p[[6]], p[[8]], p[[10]], p[[12]], p[[14]], p[[16]],
                 p[[17]], p[[19]], p[[21]], p[[23]], p[[25]], p[[27]], p[[29]], p[[31]],
                 p[[18]], p[[20]], p[[22]], p[[24]], p[[26]], p[[28]], p[[30]], p[[32]],
                 ncol = 8) &
        theme(text = element_text(colour = 'black'))
p1
PlotPDF('30.1.st_feat.example_genes', 32, 18)
p1
dev.off()


#### Fig 1D  ####
tmp.srt <- RunUMAP(drop.srt, reduction = 'harmony', dims = 1:30, min.dist = 0.5)
p2.1 <- DimPlot2(tmp.srt, cols = mycol_10[3:10], group.by = 'Cell_type', pt.size = 0.1)
p2.2 <- DimPlot2(drop_cm.srt, cols = mycol_10[1:2], group.by = 'CM_State', reduction = 'newumap', pt.size = 1.5)
PlotPDF('30.2.umap.cell_type', 10, 5)
p2.1 + p2.2
dev.off()


#### Fig 1E  ####
markers <- readRDS('external/fransico_data/analysis/outputs/overall.marker.list.rds')
markers <- rbind(
        markers[markers$cluster == 'CM1', ],
        markers[markers$cluster == 'CM2', ],
        markers[markers$cluster == 'CF', ],
        markers[markers$cluster == 'MAC', ],
        markers[markers$cluster == 'EC_Pecam1+', ],
        markers[markers$cluster == 'EC_Flt1+', ],
        markers[markers$cluster == 'SMC', ]
)
Idents(drop.srt) <- 'Cell_type2'
p3 <- MarkerHeatmap(drop.srt, marker.df = markers, disp.min = 1)
p3
PlotPDF('30.3.heat.cell_type_marker', 10, 10)
p3
dev.off()


#### Fig 2A  ####
gl <- c(
        'Decon_CM1',
        'Decon_CM2',
        'Decon_SMC',
        'Decon_CF',
        'Decon_Mac',
        'Decon_EC',
        'Decon_EndoC'
        )
p <- SpatialFeaturePlot(st_wh.srt,
                        features = gl,
                        image.alpha = 0, ncol = 7, stroke = 0)
p4 <- wrap_plots(p[[1]], p[[3]], p[[5]], p[[7]], p[[9]], p[[11]], p[[13]],
                 p[[2]], p[[4]], p[[6]], p[[8]], p[[10]], p[[12]], p[[14]],
                 ncol = 7) &
        my.scale_colour_distiller(limits = c(0, 1), palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill') &
        theme(aspect.ratio = 1)
p4
PlotPDF('30.4.st_feat.devon_result', 28, 10)
p4
dev.off()


#### Fig 2B-C  ####
Idents(drop_cm.srt) <- 'CM_State'
mk <- FindMarkers(drop_cm.srt, ident.1 = 'CM2', ident.2 = 'CM1')
mk <- mk[mk$p_val_adj < 0.05, ]
mk$gene <- rownames(mk)
gl <- split(mk$gene, mk$avg_log2FC > 0)
names(gl) <- c('aCM1', 'aCM2')

data <- read.csv(gzfile('external/fransico_data/analysis/outputs/CM1_vs_CM2_DEG.csv.gz'))
data <- data[data$p_val_adj < 0.05 & abs(data$avg_log2FC) > 0.25, ]
data$direction <- 'aCM2_upregulated'
data$direction[data$avg_log2FC < 0] <- 'aCM2_downregulated'
colnames(data)[1] <- 'gene'
gl2 <- split(data$gene, data$direction)

enrich <- ModuleEnrichment(gl2, human_or_mouse = 'mouse')
GO_CM1 <- enrich$GO$GO_aCM2_downregulated
GO_CM1 <- GO_CM1[GO_CM1$p.adjust < 0.05, ]
GO_CM2 <- enrich$GO$GO_aCM2_upregulated
GO_CM2 <- GO_CM2[GO_CM2$p.adjust < 0.05, ]
WriteCSV(GO_CM1, 'PART92.aCM1_vs_aCM2_Up_GOBP')
WriteCSV(GO_CM2, 'PART92.aCM2_vs_aCM1_Up_GOBP')
GO_CM1 <- ReadCSV('PART92.aCM1_vs_aCM2_Up_GOBP')
GO_CM2 <- ReadCSV('PART92.aCM2_vs_aCM1_Up_GOBP')

gooi <- c('sarcomere organization',
          'regulation of actin filament-based process',
          'muscle cell differentiation',
          'protein localization to cell periphery',
          'actin filament organization',
          'regulation of supramolecular fiber organization',
          'regulation of actin cytoskeleton organization',
          'muscle organ development',
          'regulation of actin filament depolymerization')
go_plot <- GO_CM2[GO_CM2$Description %in% gooi, ]
go_plot$Description <- factor(go_plot$Description, levels = go_plot$Description[order(go_plot$p.adjust)])
p5.1 <- ggplot(go_plot, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[2], alpha = 0.3) +
        geom_text(aes(label = str_to_sentence(Description), x = 0.1, y = Description), hjust = 'left') +
        labs(x = '-Log10 adjusted p value', y = 'Top Enriched Terms', title = 'GO Biological Processes') +
        scale_y_discrete(limit = rev) &
        theme_classic() &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

gooi <- c('aerobic respiration',
          'oxidative phosphorylation',
          'ATP biosynthetic process',
          'electron transport chain',
          'ATP metabolic process',
          'fatty acid oxidation',
          'heart contraction',
          'cytoplasmic translation',
          'myofibril assembly')
go_plot <- GO_CM1[GO_CM1$Description %in% gooi, ]
go_plot$Description <- factor(go_plot$Description, levels = go_plot$Description[order(go_plot$p.adjust)])
p5.2 <- ggplot(go_plot, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.3) +
        geom_text(aes(label = str_to_sentence(Description), x = 1, y = Description), hjust = 'left') +
        labs(x = '-Log10 adjusted p value', y = 'Top Enriched Terms', title = 'GO Biological Processes') +
        scale_y_discrete(limit = rev) &
        theme_classic() &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p5.1 + p5.2

gene_labels <- c(
        'Ttn', 'Sorbs2', 'Flnc', 'Acta2', 'Tmsb4x',   'Lmcd1', 'Rock2', 'Ankrd23', 'Sptbn1', 'Itga7',
        'Cox7b', 'Atp5j', 'Cycs', 'Cox7c', 'Ndufa4',     'Ybx1', 'Myl2', 'Etfa', 'Uqcrb', 'Fabp3', 'Fth1'
)
rownames(data) <- data$gene
p5.3 <- MarkerVolcano(mk.df = data, label = T, label_genes = gene_labels, line = T,
                      min.log2FC = 0.25, max.overlaps = 20) +
        coord_flip() +
        scale_x_continuous(limits = c(-3, 3)) +
        NoLegend()
p5.3
PlotPDF('30.5.barplot.aCM2_vs_aCM1_enriched_bp', 12, 6)
p5.3 + p5.1 + p5.2
dev.off()




#### Fig 2D  ####
DefaultAssay(st.srt) <- 'SCT_Fran'
# st.srt <- AddModuleScore2(st.srt, features = list(c(markers$gene[markers$cluster == 'CM2'], 'Ccn2', 'Cavin1')),
#                           names = 'aCM2_score', return_z = T, assay = 'SCT_Fran')
# st_wh.srt$aCM2_score <- NA
# st_wh.srt@meta.data[Cells(st.srt), 'aCM2_score'] <- st.srt$aCM2_score
#
# px <- SpatialFeaturePlot(st_wh.srt[, st_wh.srt$Sample == 'Control'],
#                          features = 'aCM2_score', image.alpha = 0, stroke = 0,
#                          min.cutoff = -2, max.cutoff = 0) &
#         my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill', limits = c(-2, 0))
# py <- SpatialFeaturePlot(st_wh.srt[, st_wh.srt$Sample == 'YAP5SA'],
#                          features = 'aCM2_score', image.alpha = 0, stroke = 0,
#                          min.cutoff = -1, max.cutoff = 2) &
#         my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill', limits = c(-1, 2))

st_wh.srt$Acta2_exp <- FetchData(st_wh.srt, vars = 'Acta2', slot = 'data')
st_wh.srt$Tmsb4x_exp <- FetchData(st_wh.srt, vars = 'Tmsb4x', slot = 'data')
st_wh.srt$S100a11_exp <- FetchData(st_wh.srt, vars = 'S100a11', slot = 'data')
st_wh.srt$Eef2_exp <- FetchData(st_wh.srt, vars = 'Eef2', slot = 'data')
st_wh.srt$Mprip_exp <- FetchData(st_wh.srt, vars = 'Mprip', slot = 'data')
st_wh.srt$Lmcd1_exp <- FetchData(st_wh.srt, vars = 'Lmcd1', slot = 'data')
st_wh.srt$Ahnak_exp <- FetchData(st_wh.srt, vars = 'Ahnak', slot = 'data')

gl <- c(
        'Acta2_exp',
        'Tmsb4x_exp',
        'S100a11_exp',
        'Eef2_exp',
        'Mprip_exp',
        'Lmcd1_exp',
        'Ahnak_exp'
)
max_cutoff <- c()
min_cutoff <- c()
for(i in 1:L(gl)){
        test <- st_wh.srt@meta.data[, gl[i]]
        max_cutoff[i] <- SetQuantile(cutoff = 'q95', data = test)
        min_cutoff[i] <- SetQuantile(cutoff = 'q5', data = test)
        st_wh.srt@meta.data[, gl[i]][st_wh.srt@meta.data[, gl[i]] > max_cutoff[i]] <- max_cutoff[i]
        st_wh.srt@meta.data[, gl[i]][st_wh.srt@meta.data[, gl[i]] < min_cutoff[i]] <- min_cutoff[i]
}
st_wh.srt@meta.data[!Cells(st_wh.srt) %in% Cells(st.srt), gl] <- NA

p <- SpatialFeaturePlot(st_wh.srt, features = gl, image.alpha = 0, ncol = 7, stroke = 0)
p6 <- wrap_plots(#px[[1]],
                 p[[1]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[1], max_cutoff[1])),
                 p[[3]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[2], max_cutoff[2])),
                 p[[5]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[3], max_cutoff[3])),
                 p[[7]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[4], max_cutoff[4])),
                 p[[9]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[5], max_cutoff[5])),
                 p[[11]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                     limits = c(min_cutoff[6], max_cutoff[6])),
                 p[[13]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                     limits = c(min_cutoff[6], max_cutoff[6])),

                 #py[[2]],
                 p[[2]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[1], max_cutoff[1])),
                 p[[4]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[2], max_cutoff[2])),
                 p[[6]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[3], max_cutoff[3])),
                 p[[8]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                    limits = c(min_cutoff[4], max_cutoff[4])),
                 p[[10]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                     limits = c(min_cutoff[5], max_cutoff[5])),
                 p[[12]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                     limits = c(min_cutoff[6], max_cutoff[6])),
                 p[[14]] + my.scale_colour_distiller(palette = 'Spectral', na.value = 'grey80', aesthetics = 'fill',
                                                     limits = c(min_cutoff[6], max_cutoff[6])),
                 ncol = 7) &
        theme(aspect.ratio = 1)
p6
PlotPDF('30.6.st_feat.cm2_markers', 28, 10)
p6
dev.off()


#### Fig 2E  ####
gl <- c(
        'C3_exp',
        'Cfh_exp',
        'Clu_exp',
        'C1qa_exp',
        'C1qb_exp',
        'C1qc_exp',
        'C3ar1_exp',
        'C5ar1_exp'
)
st_wh.srt$C3_exp <- FetchData(st_wh.srt, vars = 'C3', slot = 'data')
st_wh.srt$Cfh_exp <- FetchData(st_wh.srt, vars = 'Cfh', slot = 'data')
st_wh.srt$Clu_exp <- FetchData(st_wh.srt, vars = 'Clu', slot = 'data')
st_wh.srt$C1qa_exp <- FetchData(st_wh.srt, vars = 'C1qa', slot = 'data')
st_wh.srt$C1qb_exp <- FetchData(st_wh.srt, vars = 'C1qb', slot = 'data')
st_wh.srt$C1qc_exp <- FetchData(st_wh.srt, vars = 'C1qc', slot = 'data')
st_wh.srt$C3ar1_exp <- FetchData(st_wh.srt, vars = 'C3ar1', slot = 'data')
st_wh.srt$C5ar1_exp <- FetchData(st_wh.srt, vars = 'C5ar1', slot = 'data')

st_wh.srt@meta.data[!Cells(st_wh.srt) %in% Cells(st.srt), gl] <- NA

p <- SpatialFeaturePlot(st_wh.srt, features = gl, image.alpha = 0, stroke = 0, min.cutoff = 0, max.cutoff = 2.5)
p7 <- wrap_plots(
                 p[[1]], p[[3]], p[[5]], p[[7]], p[[9]], p[[11]], p[[13]], p[[15]],
                 p[[2]], p[[4]], p[[6]], p[[8]], p[[10]], p[[12]], p[[14]], p[[16]], ncol = 8) &
        my.scale_colour_distiller(palette = 'Spectral', aesthetics = 'fill', na.value = 'grey80') &
        theme(aspect.ratio = 1)
px <- SpatialFeaturePlot(st_wh.srt, features = 'Clu_exp', image.alpha = 0, stroke = 0, min.cutoff = 1, max.cutoff = 4) &
        my.scale_colour_distiller(palette = 'Spectral', aesthetics = 'fill', limits = c(1, 4), na.value = 'grey80') &
        theme(aspect.ratio = 1)

p7[[3]] <- px[[1]]
p7[[11]] <- px[[2]]
p7
PlotPDF('30.7.st_feat.complement', 28, 10)
p7
dev.off()


#### Fig 2F  ####
gl <- c(
        'C3',
        'Cfh',
        'Clu',
        'C1qa',
        'C1qb',
        'C1qc',
        'C3ar1',
        'C5ar1'
)
tmp.srt <- drop.srt[, drop.srt$Cell_type2 %in% c('aCM1', 'aCM2', 'CF', 'MP') &
                            drop.srt$Sample %in% c('Control 1', 'Control 2', 'YAP5SA 1', 'YAP5SA 2', 'YAP5SA 3')]
tmp.srt$tmp <- paste(tmp.srt$Cell_type2, tmp.srt$Experiment)
p8 <- DotPlot2(tmp.srt, features = gl, group.by = 'tmp', col.min = 0, col.max = 2) +
        scale_colour_gradient(low =c("grey90"), high =c("red"))
p8
PlotPDF('30.8.dot.complement', 6, 4)
p8
dev.off()


#### Fig 4C  ####
gl <- c(
        'C3ar1',
        'C3',
        'Cx3cr1',
        'Cx3cl1',
        'Csf1r',
        'Csf1',
        'Cd47',
        'Cd36',
        'Thbs1'
)
tmp.srt <- drop.srt[, drop.srt$Cell_type2 %in% c('aCM1', 'aCM2', 'CF', 'MP') & !is.na(drop.srt$Sample)]
tmp.srt$Experiment[is.na(tmp.srt$Experiment)] <- 'Ctrl'
tmp.srt$tmp <- paste(tmp.srt$Cell_type2, tmp.srt$Experiment)
p9 <- DotPlot2(tmp.srt, features = gl, group.by = 'tmp', col.min = -1, col.max = 2, scale.by = 'radius') +
        scale_colour_gradient(low =c("grey90"), high =c("red")) +
        scale_x_discrete(limits = rev) +
        coord_flip()
p9
PlotPDF('30.9.dot.ligand_receptor', 6, 4)
p9
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Tables  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Table S1  ####
markers <- readRDS('external/fransico_data/analysis/outputs/overall.marker.list.rds')
markers$duplicated <- NULL
markers$cluster <- revalue(markers$cluster, replace = c('EC_Flt1+' = 'EC2',
                                                        'EC_Pecam1+' = 'EC1',
                                                        'MAC' = 'MP',
                                                        'CM1' = 'aCM1',
                                                        'CM2' = 'aCM2'))
markers$pct.1 <- NULL
markers$pct.2 <- NULL
markers$p_val <- NULL
markers <- markers[! duplicated(markers$gene), ]
WriteCSV(markers, 'Table_S1.Top_cell_type_markers_for_deconvolution')

####  Table S2  ####
data <- read.csv(paste0('~/Documents/Bioinformatics/project/2022_yap5sa_st_rli/meta/mouse_v0/',
                        '91_NCR_Revision_1/part91.cm1_cm2_transition_correlated_genes.csv'))
WriteCSV(data, 'Table_S2.aCM1_to_aCM2_trajectory_correlated_genes')

####  Table S3  ####
data <- read.csv(gzfile('external/fransico_data/analysis/outputs/CM1_vs_CM2_DEG.csv.gz'))
data <- data[data$p_val_adj < 0.05 & abs(data$avg_log2FC) > 0.25, ]
data$direction <- 'aCM2_upregulated'
data$direction[data$avg_log2FC < 0] <- 'aCM2_downregulated'
colnames(data)[1] <- 'gene'
WriteCSV(data, 'Table_S3.aCM2_vs_aCM1_differnetially_express_genes')

####  Table S4  ####
zone_mk <- readRDS('~/Documents/Bioinformatics/r/db/adult_mouse_post_mi_zone_marker_lisftover.list.rds')
zone_mk2 <- data.frame(Gene = c(zone_mk$RZ,
                                c(zone_mk$BZ1, zone_mk$BZ2),
                                zone_mk$IZ),
                       Zone = c(rep('Remote_zone', L(zone_mk$RZ)),
                                rep('Border_zone', L(c(zone_mk$BZ1, zone_mk$BZ2))),
                                rep('Ischemic_zone', L(zone_mk$IZ)))
)

WriteCSV(zone_mk2, 'Table_S4.MI_zone_markers_for_annotating_human_mi_ST_data')

####  Table S5  ####
## Previous Table S2

####  Table S6  ####
## Previous Table S3

####  Table S7  ####
st_wh_ctrl.srt <- readRDS('analysis/part92.ST_YAP5SA_wholeheart_ctrl_niche_annotated.seurat.rds')
st_wh_y5sa.srt <- readRDS('analysis/part92.ST_YAP5SA_wholeheart_y5sa_niche_annotated.seurat.rds')
Idents(st_wh_ctrl.srt) <- 'Niche_split'
ctrl_mk <- FindAllMarkers(st_wh_ctrl.srt, only.pos = T, assay = 'ST')
ctrl_mk <- ctrl_mk[ctrl_mk$p_val_adj < 0.05 ,]
Table(ctrl_mk$cluster)
ctrl_mk$slide <- 'Control'

Idents(st_wh_y5sa.srt) <- 'Niche_split'
y5sa_mk <- FindAllMarkers(st_wh_y5sa.srt, only.pos = T, assay = 'ST')
y5sa_mk <- y5sa_mk[y5sa_mk$p_val_adj < 0.05 ,]
Table(y5sa_mk$cluster)
y5sa_mk$slide <- 'YAP5SA'

data <- rbind(ctrl_mk, y5sa_mk)
data$p_val <- NULL
data$pct.1 <- NULL
data$pct.2 <- NULL
WriteCSV(data, 'Table_S7.Spatial_niche_marker_genes')

####  Table S8  ####
## Previous Table S4

####  Table S9  ####
data <- read.csv(paste0('~/Documents/Bioinformatics/project/2023_neoc3ko_rli/meta/mouse_v0/',
                        'PART91_NCVR_Revision/PART91.CF_DEG.C3KO_MI_vs_WT_MI.csv'))
data$p_val <- NULL
data$direction <- 'C3_null_upregulated'
data$direction[data$avg_log2FC < 0] <- 'C3_null_downregulated'
WriteCSV(data, 'Table_S9.C3_null_MI_vs_WT_MI_differentially_expressed_genes_in_CF_Cluster_1')

####  Table S10  ####
data <- read.csv(paste0('~/Documents/Bioinformatics/project/2023_neoc3ko_rli/meta/mouse_v0/',
                        'PART91_NCVR_Revision/PART91.CM_DEG.C3KO_MI_vs_WT_MI.csv'))
data$p_val <- NULL
WriteCSV(data, 'Table_S10.C3_null_MI_vs_WT_MI_differentially_expressed_genes_in_all_CM_clusters')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
