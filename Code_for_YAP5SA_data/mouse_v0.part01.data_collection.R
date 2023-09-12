####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate directories  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- '01_Data_Collection'
Project <- '2022_yap5sa_st_rli'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'mouse', Project, 'ithil')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 1: αMHC-mcm:YAP5SA ST data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ctrl <- readRDS('/Volumes/shire/data/visium/2020_Yap5sa/control_st_whole_heart.rds')
ctrl <- RenameAssays(ctrl, 'Spatial' = 'ST')
DefaultAssay(ctrl) <- 'ST'
names(ctrl@images) <- 'Control'
ctrl@images$Control@key <- 'Control_'
ctrl@images$Control@assay <- 'ST'

ctrl$Orig_name <- Cells(ctrl)
ctrl$Sample <- 'Control'
ctrl$Sample_id <- 'Control'
ctrl$Genotype <- 'WT'
ctrl$Genotype_long <- 'αMHC-mcm'
ctrl$Coord_x_slide <- ctrl@images$Control@coordinates$imagecol
ctrl$Coord_y_slide <- -ctrl@images$Control@coordinates$imagerow
ctrl <- RenameCells(ctrl, new.names = paste0(ctrl$Sample_id, '_', Cells(ctrl)))

yap <- readRDS('/Volumes/shire/data/visium/2020_Yap5sa/mutant_st_whole_heart.rds')
yap <- RenameAssays(yap, 'Spatial' = 'ST')
DefaultAssay(yap) <- 'ST'
names(yap@images) <- 'YAP5SA'
yap@images$YAP5SA@key <- 'YAP5SA_'
yap@images$YAP5SA@assay <- 'ST'

yap$Orig_name <- Cells(yap)
yap$Sample <- 'YAP5SA'
yap$Sample_id <- 'YAP5SA'
yap$Genotype <- 'YAP5SA'
yap$Genotype_long <- 'αMHC-mcm:YAP5SA'
yap$Coord_x_slide <- yap@images$YAP5SA@coordinates$imagecol
yap$Coord_y_slide <- -yap@images$YAP5SA@coordinates$imagerow
yap <- RenameCells(yap, new.names = paste0(yap$Sample_id, '_', Cells(yap)))

## SCT normalize each slide
srt.list <- list(ctrl, yap)
for(i in 1:2){
        srt.list[[i]] <- SCTransform(srt.list[[i]],
                                     assay = 'ST',
                                     method = "glmGamPoi",
                                     vars.to.regress = c('nCount_ST', 'nFeature_ST'),
                                     do.correct.umi = T,
                                     seed.use = 505,
                                     return.only.var.genes = F)
}
## Merge dataset
merged <- merge(srt.list[[1]], srt.list[[2]])
VariableFeatures(merged) <- U(c(VariableFeatures(srt.list[[1]], assay = 'SCT'),
                                VariableFeatures(srt.list[[2]], assay = 'SCT')))

## Edit metadata
merged$Study <- '2020_YAP5SA'
merged$Age <- 'Unknown'
merged$Age_group <- 'Adult'
merged$Sex <- 'Unknown'
merged$Condition <- 'Normal'
merged$Condition_group <- 'Normal'
merged$Method <- 'Visium'
merged$Tissue <- 'Whole_heart'
merged$orig.ident <- NULL
merged$nCount_Spatial <- NULL
merged$nFeature_Spatial <- NULL
merged$percent.mt <- NULL
merged$Sample <- factor(merged$Sample, levels = c('Control', 'YAP5SA'))

## Store spatial coordinates as a reduced dimension
merged <- RunPCA(merged, verbose = F, npcs = 2)
merged@reductions$spatial <- merged@reductions$pca ## create a dimension
merged@reductions$spatial@key <- 'Spatial_'
colnames(merged@reductions$spatial@cell.embeddings) <- c('Spatial_1', 'Spatial_2')
merged@reductions$spatial@cell.embeddings[,1] <- merged$Coord_x_slide
merged@reductions$spatial@cell.embeddings[,2] <- merged$Coord_y_slide
merged@misc$spot_scale <- c(1, 1)

p1 <- FeaturePlotST(merged, features = 'Dcn',
                    minvals = rep(0, 4),
                    maxvals = rep(5, 4),
                    pt.sizes = merged@misc$spot_scale,
                    ncol = 2)
p2 <- SpatialFeaturePlot(merged, features = 'Dcn', image.alpha = 0, ncol = 2, pt.size.factor = 1.5) & NoLegend()

PlotPDF('Part01.01.st_feat.2020_YAP5SA', 8, 8)
print(p1/p2)
dev.off()

# ## Label zones
# merged2 <- RunPCA(merged) |> RunUMAP(dims = 1:30)
# merged2 <- FindNeighbors(merged2) |> FindClusters(res = 0.05)
# FeaturePlot2(merged2, features = 'Myl7') +
#         DimPlot2(merged2, group.by = 'Sample') +
#         DimPlot2(merged2, group.by = 'SCT_snn_res.0.05')
# merged$Zone <- 'Ventricles'
# merged$Zone[merged2$SCT_snn_res.0.05 %in% c(1)] <- 'Atria'

## Consistency with Francisco's version
x <- read.csv(gzfile('external/fransico_data/analysis/outputs/st_integrated_metadata.csv.gz'))
ctrl_name <- x$X[x$genotype == 'Control']
ctrl_name <- paste0('Control_', str_split(ctrl_name, pattern = '_', simplify = T)[,1])
y5sa_name <- x$X[x$genotype != 'Control']
y5sa_name <- paste0('YAP5SA_', str_split(y5sa_name, pattern = '_', simplify = T)[,1])

O(x = ctrl_name, y = as.vector(Cells(merged)[merged$Sample == 'Control']), uniq = F)
O(x = y5sa_name, y = as.vector(Cells(merged)[merged$Sample == 'YAP5SA']), uniq = F)

merged$Zone <- 'Atrium'
merged$Zone[c(ctrl_name, y5sa_name)] <- 'Ventricle'

p <- SpatialDimPlot(merged, group.by = 'Zone')
p
PlotPDF('Part01.02.st_dim.2020_YAP5SA_francisco_zone', 8, 4)
print(p)
dev.off()

## Add Spotlight decon result from Francisco
x$name <- str_split(x$X, pattern = '_', simplify = T)[,1]
x$name[x$genotype == 'Control'] <- paste0('Control_', x$name[x$genotype == 'Control'])
x$name[x$genotype == 'Mutant'] <- paste0('YAP5SA_', x$name[x$genotype == 'Mutant'])
rownames(x) <- x$name
O(Cells(merged), as.vector(rownames(x)), uniq = F)

merged$Decon_CF <- x[Cells(merged), 'CF']
merged$Decon_CM1 <- x[Cells(merged), 'CM1']
merged$Decon_CM2 <- x[Cells(merged), 'CM2']
merged$Decon_EC <- x[Cells(merged), 'EC.Flt1.']
merged$Decon_EndoC <- x[Cells(merged), 'EC.Pecam1.']
merged$Decon_SMC <- x[Cells(merged), 'SMC']
merged$Decon_Mac <- x[Cells(merged), 'MAC']

p <- SpatialFeaturePlot(merged, stroke = 0, image.alpha = 0.5,
                        features = c('Decon_CF', 'Decon_CM1', 'Decon_CM2', 'Decon_EndoC',
                                     'Decon_EC', 'Decon_SMC', 'Decon_Mac'))
PlotPDF('Part01.03.st_feat.2020_YAP5SA_francisco_spotlight', 6, 21)
print(p)
dev.off()

## Add expression matrix from Francisco
merged_vent <- merged[, merged$Zone == 'Ventricle']
fran_st.srt <- readRDS('external/fransico_data/analysis/outputs/st_integrated.rds')

orig_name <- Cells(fran_st.srt)
orig_name[fran_st.srt$genotype == 'Control'] <- paste0('Control_',
                                                       str_split(orig_name[fran_st.srt$genotype == 'Control'],
                                                                 pattern = '_',
                                                                 simplify = T)[,1])
orig_name[fran_st.srt$genotype != 'Control'] <- paste0('YAP5SA_',
                                                       str_split(orig_name[fran_st.srt$genotype != 'Control'],
                                                                 pattern = '_',
                                                                 simplify = T)[,1])
identical(orig_name, as.vector(Cells(merged_vent)))
fran_st.srt <- RenameCells(fran_st.srt, new.names = orig_name)

identical(merged_vent@assays$ST@data, fran_st.srt@assays$Spatial@data)
merged_vent[['SCT_Fran']] <- fran_st.srt@assays$SCT

## Run imputation on SCT
merged_vent <- RunALRA(merged_vent, assay = 'SCT', genes.use = rownames(merged_vent))
DefaultAssay(merged_vent) <- 'SCT'

## Export meta and matrices for GEO
merged_vent$Orig_name <- NULL
merged_vent$Sample_id <- NULL
merged_vent$Study <- NULL
merged_vent$Age <- NULL
merged_vent$Sex <- NULL
merged_vent$Condition <- NULL
merged_vent$Condition_group <- NULL
merged_vent$nCount_SCT_Fran <- NULL
merged_vent$nFeature_SCT_Fran <- NULL
merged_vent$Genotype <- merged_vent$Genotype_long
merged_vent$Genotype_long <- NULL

ctrl <- merged_vent[, merged_vent$Sample == 'Control']
write.table(GetAssayData(ctrl, slot = 'count', assay = 'ST'),
            file = '~/Desktop/ST_control.raw_count.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
write.table(GetAssayData(ctrl, slot = 'data', assay = 'SCT_Fran'),
            file = '~/Desktop/ST_control.normalized_count.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
write.table(ctrl@meta.data,
            file = '~/Desktop/ST_control.meta_data.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
y5sa <- merged_vent[, merged_vent$Sample == 'YAP5SA']
write.table(GetAssayData(y5sa, slot = 'count', assay = 'ST'),
            file = '~/Desktop/ST_yap5sa.raw_count.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
write.table(GetAssayData(y5sa, slot = 'data', assay = 'SCT_Fran'),
            file = '~/Desktop/ST_yap5sa.normalized_count.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
write.table(y5sa@meta.data,
            file = '~/Desktop/ST_yap5sa.meta_data.csv',
            quote = F, sep = ',', row.names = T, col.names = T)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged, 'individual/part01.YAP5SA_st.srt.rds')
saveRDS(merged_vent, 'individual/part01.YAP5SA_st_ventricle.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 2: 2021_NatComm_EOlson ST data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## P1ShamD3
P1ShamD3 <- readRDS('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public/P1ShamD3/P1ShamD3_spatial.rds')
SpatialFeaturePlot(P1ShamD3, features = 'nFeature_Spatial')

DefaultAssay(P1ShamD3) <- 'Spatial'
P1ShamD3 <- DietSeurat(P1ShamD3, counts = T, data = T, scale.data = F, assays = 'Spatial')
P1ShamD3 <- RenameAssays(P1ShamD3, 'Spatial' = 'ST')
P1ShamD3@images$P1ShamD3 <- P1ShamD3@images$slice1
P1ShamD3@images$P1ShamD3@key <- 'P1ShamD3_'
P1ShamD3@images$P1ShamD3@assay <- 'ST'
P1ShamD3$Coord_x_slide <- -P1ShamD3@images$P1ShamD3@coordinates$imagerow
P1ShamD3$Coord_y_slide <- P1ShamD3@images$P1ShamD3@coordinates$imagecol
P1ShamD3@images$slice1 <- NULL
P1ShamD3@images$P1ShamD3@coordinates$imagerow <- P1ShamD3$Coord_x_slide
P1ShamD3@images$P1ShamD3@coordinates$imagecol <- P1ShamD3$Coord_y_slide
SpatialFeaturePlot(P1ShamD3, features = 'Plvap')

## P1MID3
P1MID3 <- readRDS('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public/P1MID3/P1MID3_spatial.rds')
SpatialFeaturePlot(P1MID3, features = 'nFeature_Spatial')
SpatialDimPlot(P1MID3, cells.highlight = Cells(P1MID3)[P1MID3$nFeature_Spatial<=2000])
P1MID3 <- P1MID3[, P1MID3$nFeature_Spatial>2000]

DefaultAssay(P1MID3) <- 'Spatial'
P1MID3 <- DietSeurat(P1MID3, counts = T, data = T, scale.data = F, assays = 'Spatial')
P1MID3 <- RenameAssays(P1MID3, 'Spatial' = 'ST')
P1MID3@images$P1MID3 <- P1MID3@images$slice1
P1MID3@images$P1MID3@key <- 'P1MID3_'
P1MID3@images$P1MID3@assay <- 'ST'
P1MID3$Coord_x_slide <- -P1MID3@images$P1MID3@coordinates$imagecol
P1MID3$Coord_y_slide <- -P1MID3@images$P1MID3@coordinates$imagerow
P1MID3@images$slice1 <- NULL

## P1ShamD7
P1ShamD7 <- readRDS('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public/P1ShamD7/P1ShamD7_spatial.rds')
SpatialFeaturePlot(P1ShamD7, features = 'nFeature_Spatial')
SpatialDimPlot(P1ShamD7, cells.highlight = Cells(P1ShamD7)[P1ShamD7$nFeature_Spatial<=1000])
P1ShamD7 <- P1ShamD7[, P1ShamD7$nFeature_Spatial>1000]

DefaultAssay(P1ShamD7) <- 'Spatial'
P1ShamD7 <- DietSeurat(P1ShamD7, counts = T, data = T, scale.data = F, assays = 'Spatial')
P1ShamD7 <- RenameAssays(P1ShamD7, 'Spatial' = 'ST')
P1ShamD7@images$P1ShamD7 <- P1ShamD7@images$slice1
P1ShamD7@images$P1ShamD7@key <- 'P1ShamD7_'
P1ShamD7@images$P1ShamD7@assay <- 'ST'
P1ShamD7$Coord_x_slide <- P1ShamD7@images$P1ShamD7@coordinates$imagecol
P1ShamD7$Coord_y_slide <- P1ShamD7@images$P1ShamD7@coordinates$imagerow
P1ShamD7@images$slice1 <- NULL

## P1MID7 rep2
count <- read.csv(paste0('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public_error/',
                         'GSM5268647_P1MID7_rep2/GSM5268647_P1MID7.rep2_counts.csv'), header = T)
count <- count[!duplicated(count$X), ]
rownames(count) <- count$X
count$X <- NULL
meta <- read.csv(paste0('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public_error/',
                        'GSM5268647_P1MID7_rep2/GSM5268647_P1MID7.rep2.metadata.csv'), header = T)
rownames(meta) <- str_replace(meta$X, '-', replacement = '.')
identical(rownames(meta), colnames(count))
P1MID7 <- CreateSeuratObject(counts = count, meta.data = meta, assay = 'ST')
P1MID7 <- P1MID7[, P1MID7$imagerow >= 100 & P1MID7$imagerow <= 350 & P1MID7$imagecol >= 200 & P1MID7$imagecol <= 500]
P1MID7$Coord_x_slide <- P1MID7$imagecol
P1MID7$Coord_y_slide <- P1MID7$imagerow
P1MID7@images$P1MID7 <- new(## add image info
        Class = 'VisiumV1',
        image = png::readPNG(paste0('/Volumes/shire/data/visium/2021_NatComm_EOlson/matrix_public_error',
                                    '/GSM5268647_P1MID7_rep2/GSM5268647_P1MID7.rep2.tissue_hires_image.png')),
        scale.factors = scalefactors(spot = 0.2, fiducial = 100, hires = 1, lowres = 2), # fake scale factors
        coordinates = data.frame('imagerow' = -P1MID7$imagerow, 'imagecol' = P1MID7$imagecol),
        spot.radius = 0.01, # fake radius size
        key = 'P1MID7_',
        assay = 'ST'
)

## Add metadata:
P1ShamD3$Sample <- 'P1ShamD3'
P1ShamD3$Sample_id <- 'D001S001'
P1ShamD3$Orig_name <- Cells(P1ShamD3)
P1ShamD3$Condition <- '03 dpSham'
P1ShamD3$Age <- 'P004'
P1ShamD3 <- RenameCells(P1ShamD3,
                        new.names = paste0(P1ShamD3$Sample_id, '_',
                                           str_replace_all(Cells(P1ShamD3), pattern = '\\.', replacement = '-')))
P1MID3$Sample <- 'P1MID3'
P1MID3$Sample_id <- 'D001S002'
P1MID3$Orig_name <- Cells(P1MID3)
P1MID3$Condition <- '03 dpMI'
P1MID3$Age <- 'P004'
P1MID3 <- RenameCells(P1MID3,
                      new.names = paste0(P1MID3$Sample_id, '_',
                                         str_replace_all(Cells(P1MID3), pattern = '\\.', replacement = '-')))
P1ShamD7$Sample <- 'P1ShamD7'
P1ShamD7$Sample_id <- 'D001S003'
P1ShamD7$Orig_name <- Cells(P1ShamD7)
P1ShamD7$Condition <- '07 dpSham'
P1ShamD7$Age <- 'P008'
P1ShamD7 <- RenameCells(P1ShamD7,
                        new.names = paste0(P1ShamD7$Sample_id, '_',
                                           str_replace_all(Cells(P1ShamD7), pattern = '\\.', replacement = '-')))
P1MID7$Sample <- 'P1MID7'
P1MID7$Sample_id <- 'D001S004'
P1MID7$Orig_name <- Cells(P1MID7)
P1MID7$Condition <- '07 dpMI'
P1MID7$Age <- 'P008'
P1MID7 <- RenameCells(P1MID7,
                      new.names = paste0(P1MID7$Sample_id, '_',
                                         str_replace_all(Cells(P1MID7), pattern = '\\.', replacement = '-')))

## SCT normalize each slide
srt.list <- list(P1ShamD3, P1MID3, P1ShamD7, P1MID7)
for(i in 1:4){
        srt.list[[i]] <- SCTransform(srt.list[[i]],
                                     assay = 'ST',
                                     method = "glmGamPoi",
                                     vars.to.regress = c('nCount_ST', 'nFeature_ST'),
                                     do.correct.umi = T,
                                     seed.use = 505,
                                     return.only.var.genes = F)
}

## Correct SCT matrix for P1MID7_rep2
srt.list[[4]]@assays$SCT@counts <- srt.list[[4]]@assays$ST@counts
srt.list[[4]]@assays$SCT@data <- log1p(srt.list[[4]]@assays$SCT@counts) ## according to SCTransform()
srt.list[[4]]$nCount_SCT <- srt.list[[4]]$nCount_ST
srt.list[[4]]$nFeature_SCT <- srt.list[[4]]$nFeature_ST
srt.list[[4]]$nCount_ST <- srt.list[[4]]$nCount_Spatial
srt.list[[4]]$nFeature_ST <- srt.list[[4]]$nFeature_Spatial

## Merge dataset
merged <- merge(srt.list[[1]], srt.list[2:4])
Idents(merged) <- 'Sample'
VlnPlot2(merged, features = c('nCount_Spatial', 'nCount_ST', 'nCount_SCT'))

VariableFeatures(merged) <- U(c(VariableFeatures(srt.list[[1]], assay = 'SCT'),
                                VariableFeatures(srt.list[[2]], assay = 'SCT'),
                                VariableFeatures(srt.list[[3]], assay = 'SCT'),
                                VariableFeatures(srt.list[[4]], assay = 'SCT')))
## Edit metadata
merged$Study <- '2021_NatComm_EOlson'
merged$Age_group <- 'Neonate'
merged$Method <- 'Visium'
merged$Tissue <- 'Ventricle'
merged$Sex <- 'Unknown'
merged$Condition_group <- 'MI'
merged$Genotype <- 'WT'
merged$Genotype_long <- 'WT'
merged$Sample <- factor(merged$Sample, levels = c('P1ShamD3', 'P1MID3', 'P1ShamD7', 'P1MID7'))
merged$orig.ident <- NULL
merged$X <- NULL
merged$nCount_Spatial <- NULL
merged$nFeature_Spatial <- NULL
merged$CM4 <- NULL
merged$CM1 <- NULL
merged$EndoEC <- NULL
merged$CM3 <- NULL
merged$FB <- NULL
merged$EPI <- NULL
merged$Macrophage <- NULL
merged$Pericyte.SMC <- NULL
merged$EC <- NULL
merged$CM2 <- NULL
merged$CM5 <- NULL
merged$imagerow <- NULL
merged$imagecol <- NULL
merged$SCT_snn_res.2 <- NULL
merged$seurat_clusters <- NULL

## Store spatial coordinates as a reduced dimension
merged <- RunPCA(merged, verbose = F, npcs = 2)
merged@reductions$spatial <- merged@reductions$pca ## create a dimension
merged@reductions$spatial@key <- 'Spatial_'
colnames(merged@reductions$spatial@cell.embeddings) <- c('Spatial_1', 'Spatial_2')
merged@reductions$spatial@cell.embeddings[,1] <- merged$Coord_x_slide
merged@reductions$spatial@cell.embeddings[,2] <- merged$Coord_y_slide
merged@misc$spot_scale <- c(2.48, 3.20, 1.68, 2.25)

p1 <- FeaturePlotST(merged, features = 'Plvap',
                    minvals = rep(0, 4),
                    maxvals = rep(5, 4),
                    pt.sizes = merged@misc$spot_scale,
                    ncol = 4)
p2 <- SpatialFeaturePlot(merged, features = 'Plvap', image.alpha = 0, ncol = 4, pt.size.factor = 4) & NoLegend()
p3 <- SpatialFeaturePlot(merged, features = 'Plvap', alpha = 0, ncol = 4, pt.size.factor = 4) & NoLegend()

PlotPDF('Part01.04.st_feat.2021_NatComm_EOlson', 16, 12)
print(p1/p2/p3)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged, 'individual/part01.2021_NatComm_EOlson.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 3: 2020_DevCell_EOlson snRNA-seq data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sys.setenv("VROOM_CONNECTION_SIZE" = 1024*1024)
all.df <- read_delim(gzfile(paste0('/Volumes/shire/data/scrnaseq/2020_DevCell_EOlson/matrix_public/',
                                   'GSE130699_P1_8_AllCM.txt.gz')))
## shift colnames
all.mtx <- as.data.frame(all.df)
all.mtx$last <- as.numeric(str_split(all.mtx[, ncol(all.mtx)], pattern = ' ', simplify = T)[, 2])
all.mtx[, ncol(all.mtx)-1] <- as.numeric(str_split(all.mtx[, ncol(all.mtx)-1], pattern = ' ', simplify = T)[, 1])
colnames(all.mtx) <- c('Gene', colnames(all.mtx))[1:ncol(all.mtx)]
H(all.mtx)
all.mtx[1:5, (ncol(all.mtx)-5):ncol(all.mtx)]
genes <- all.mtx$Gene
all.mtx <- as.matrix(all.mtx[, 2:ncol(all.mtx)])
rownames(all.mtx) <- genes

cm.srt <- CreateSeuratObject(counts = all.mtx, assay = 'RNA', min.cells = 0, min.features = 0)
cm.srt$Sample <- paste0(str_split(Cells(cm.srt), '_', simplify = T, n = 2)[, 1],
                        '_',
                        str_split(Cells(cm.srt), '_', simplify = T, n = 3)[, 2])
cm.srt$Sample <- factor(cm.srt$Sample, levels = c(
        'P1_1Sham','P1_1MI','P1_3Sham','P1_3MI','P8_1Sham','P8_1MI','P8_3Sham','P8_3MI'))
cm.srt <- PercentageFeatureSet(cm.srt, pattern = '^mt-', col.name = 'pct_mito')
VlnPlot2(cm.srt, features = c('nCount_RNA', 'nFeature_RNA', 'pct_mito'), group.by = 'Sample') ## data already filtered

cm.srt <- cm.srt |> NormalizeData()
cm.srt@assays$RNA@data <- cm.srt@assays$RNA@counts ## data already normalized

## Follow Seurat integration pipeline
cm.srt.list <- SplitObject(cm.srt, split.by = "Sample")
cm.srt.list <- lapply(X = cm.srt.list, FUN = function(x) {
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = cm.srt.list)
cm.srt.list <- lapply(X = cm.srt.list, FUN = function(x) {
        x <- ScaleData(x, features = features, vars.to.regress = c('nCount_RNA', 'nFeature_RNA'))
        x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = cm.srt.list, anchor.features = features, reduction = "rpca")
cm.srt2 <- IntegrateData(anchorset = anchors)
DefaultAssay(cm.srt2) <- 'integrated'

cm.srt3 <- ScaleData(cm.srt2, verbose = FALSE)
cm.srt3 <- RunPCA(cm.srt3, npcs = 30, verbose = FALSE)
cm.srt3 <- RunUMAP(cm.srt3, reduction = "pca", dims = 1:30)
cm.srt3 <- FindNeighbors(cm.srt3, reduction = "pca", dims = 1:30)
cm.srt3 <- FindClusters(cm.srt3, resolution = 0.5)

## Evaluate
markers_lvl1 <- c('Tnnt2',  'Nppa',           # CM
                  'Col1a1', 'Tcf21',          # CF
                  'Pecam1',  'Cdh5',          # EC
                  'Prox1',   'Flt4',          # LEC
                  'Acta2',   'Rgs5',          # SMC
                  'Pdgfrb', 'Vtn',            # Pericyte
                  'Upk3b',  'Wt1',            # Epicardium
                  'Npr3',   'H19',            # Endocardium
                  'Plp1',   'Nrn1',           # Schwann
                  'Nrsn1',  'Npy',            # Neuronal
                  'Tenm4',  'Car3',           # Adipocyte
                  'Ptprc',  'H2-D1',          # Immune
                  'Hba-a1', 'Hbb-bh1',        # RBC
                  'Mki67',  'Top2a')          # Mitotic

DefaultAssay(cm.srt3) <- 'RNA'
p1 <- FeaturePlot2(cm.srt3, features = markers_lvl1, raster = T)

PlotPDF('Part01.05.feat.main_markers', 18, 14)
print(p1)
dev.off()

cm_subtype_mk <- read_excel('/Volumes/shire/data/scrnaseq/2020_DevCell_EOlson/doc/CM1-5_markers.xlsx')
cm_subtype_mk <- cm_subtype_mk[cm_subtype_mk$p_val_adj<0.0001, ]
genelist <- split(cm_subtype_mk$gene, cm_subtype_mk$cluster)
cm.srt3 <- AddModuleScore2(cm.srt3, features = genelist, assay = 'RNA', return_z = T,
                           names = paste0('Zscore_olson_', names(genelist)))

p2 <- FeaturePlot2(cm.srt3, features = paste0('Zscore_olson_', names(genelist)), raster = T, pt.size = 1) +
        DimPlot2(cm.srt3, cols = mycol_20, raster = T, pt.size = 1, label = T) + NoLegend()
PlotPDF('Part01.06.feat.cui_cm_markers', 9, 6)
print(p2)
dev.off()

cm.srt4 <- cm.srt3[, cm.srt3$seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 10)]
DefaultAssay(cm.srt4) <- 'integrated'
cm.srt4 <- FindVariableFeatures(cm.srt4, selection.method = "vst", nfeatures = 2000)
cm.srt4 <- RunPCA(cm.srt4, npcs = 30, verbose = FALSE, features = VariableFeatures(cm.srt4))
cm.srt4 <- RunUMAP(cm.srt4, reduction = "pca", dims = 1:15)
cm.srt4 <- FindNeighbors(cm.srt4, reduction = "pca", dims = 1:15)
cm.srt4 <- FindClusters(cm.srt4, resolution = 0.8)

FeaturePlot2(cm.srt4, features = paste0('Zscore_olson_', names(genelist)), raster = T, pt.size = 1) +
        DimPlot2(cm.srt4, cols = mycol_20, raster = T, pt.size = 1, label = T) + NoLegend()

cm.srt4$CM_subtype <- factor('CM1', levels = paste0('CM', 1:5))
cm.srt4$CM_subtype[cm.srt4$integrated_snn_res.0.8 %in% c(10)] <- 'CM2'
cm.srt4$CM_subtype[cm.srt4$integrated_snn_res.0.8 %in% c(5)] <- 'CM3'
cm.srt4$CM_subtype[cm.srt4$integrated_snn_res.0.8 %in% c(11,6,8)] <- 'CM4'
cm.srt4$CM_subtype[cm.srt4$integrated_snn_res.0.8 %in% c(4)] <- 'CM5'
Table(cm.srt4$CM_subtype)

Idents(cm.srt4) <- 'CM_subtype'
p3 <- FeaturePlot2(cm.srt4, features = paste0('Zscore_olson_', names(genelist)), raster = T, pt.size = 1) +
        DimPlot2(cm.srt4, cols = mycol_10, raster = T, pt.size = 1, label = T) + NoLegend()
PlotPDF('Part01.07.feat.cui_cm_markers', 9, 6)
print(p3)
dev.off()

cm.srt4 <- AddModuleScore2(cm.srt4, features = Yap5sa_CM_Mac_CF_marker[1:2], assay = 'RNA', return_z = T,
                           names = paste0('Zscore_y5sa_', names(Yap5sa_CM_Mac_CF_marker)[1:2]))
p4 <- FeaturePlot2(cm.srt4, features = paste0('Zscore_y5sa_', names(Yap5sa_CM_Mac_CF_marker)[1:2]),
                   raster = T, pt.size = 1, min.cutoff = 'q5', max.cutoff = 'q95') +
        DimPlot2(cm.srt4, cols = mycol_10, raster = T, pt.size = 1, label = T) + NoLegend()
PlotPDF('Part01.08.feat.y5sa_cm_markers', 6, 6)
print(p4)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(cm.srt4, 'individual/part01.2020_DevCell_EOlson.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 5: Yap5sa Drop-seq data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('external/fransico_data/reference_scrnaseq/Reference_combined_compendium_&_yap5sa_CMs_v2.rds')
DimPlot2(srt, group.by = 'celltype_lvl2')
DimPlot2(srt, group.by = 'pub')

cm.srt <- srt[, srt$celltype_lvl2 == 'cardiomyocyte' & srt$experiment == '2021_Unpublished_JMartin__Ventricule_']
DefaultAssay(cm.srt) <- 'RNA'
cm.srt <- NormalizeData(cm.srt)
cm.srt <- AddModuleScore2(cm.srt, features = Yap5sa_CM_Mac_CF_marker[1:2], assay = 'RNA', return_z = T,
                          names = paste0('Zscore_y5sa_', names(Yap5sa_CM_Mac_CF_marker)[1:2]))
cm.srt <- AddModuleScore2(cm.srt, features = genelist, assay = 'RNA', return_z = T,
                          names = paste0('Zscore_olson_', names(genelist)))
cm.srt <- AddModuleScore2(cm.srt, features = Yap_target[1:2], assay = 'RNA', return_z = T,
                          names = c('Zscore_yap_cm_target1', 'Zscore_yap_cm_target2'))
FeaturePlot2(cm.srt, features = c('Zscore_yap_cm_target1', 'Zscore_yap_cm_target2'))
DimPlot2(cm.srt, group.by = 'celltype') + DimPlot2(cm.srt, group.by = 'Experiment')

cm.srt <- RunUMAP(cm.srt, dims = 1:8, reduction = 'harmony')
cm.srt <- FindNeighbors(cm.srt, reduction = "harmony", dims = 1:8)
cm.srt <- FindClusters(cm.srt, resolution = 0.8, graph.name = 'SCT_snn')

DimPlot2(cm.srt, group.by = 'SCT_snn_res.0.8', label = T)
VlnPlot2(cm.srt, group.by = 'SCT_snn_res.0.8', features = c('Zscore_yap_cm_target1', 'Zscore_yap_cm_target2'))

cm.srt$New_CM_State <- factor('CM1', levels = c('CM1', 'CM2'))
cm.srt$New_CM_State[cm.srt$SCT_snn_res.0.8 %in% c(3, 5, 9)] <- 'CM2'
Idents(cm.srt) <- 'New_CM_State'
DimPlot2(cm.srt) +
DimPlot2(cm.srt, group.by = 'Experiment')

cm.srt$CM_State <- droplevels(cm.srt$celltype)
cm.srt$CM_State <- factor(cm.srt$CM_State, levels = c('CM1', 'CM2'))
colnames(cm.srt@meta.data)
cm.srt@meta.data <- cm.srt@meta.data[, c('Experiment',
                                         "Zscore_y5sa_CM1",
                                         "Zscore_y5sa_CM2",
                                         "Zscore_olson_CM1",
                                         "Zscore_olson_CM2",
                                         "Zscore_olson_CM3",
                                         "Zscore_olson_CM4",
                                         "Zscore_olson_CM5",
                                         "New_CM_State",
                                         'CM_State')]
Idents(cm.srt) <- 'New_CM_State'
new_mk <- FindAllMarkers(cm.srt, return.thresh = 0.001, logfc.threshold = 0.5, only.pos = T)
new_mk_top <- new_mk |> group_by(cluster) |> top_n(n = 20, wt = avg_log2FC)
new_mk_top <- split(new_mk_top$gene, new_mk_top$cluster)
cm.srt@misc$markers <- list(new_mk, new_mk_top)

DimPlot2(cm.srt, group.by = 'CM_State')
cm.srt <- RunUMAP(cm.srt, dims = 1:30, reduction = 'harmony', reduction.name = 'oriumap', reduction.key = 'ORIUMAP',
                  min.dist = 0.1)
DimPlot2(cm.srt, reduction = 'oriumap', raster = T, pt.size = 1.5, group.by = 'CM_State') +
        labs(x = 'UMAP1', y = 'UMAP2')
cm.srt <- RunALRA(cm.srt, assay = 'RNA', setDefaultAssay = F)

srt <- srt[, is.na(srt$pub)]
DimPlot2(srt)
srt <- DietSeurat(srt, assays = c('SCT', 'RNA'), dimreducs = names(srt@reductions))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(cm.srt, 'individual/part01.y5sa_dropseq_cm.srt.rds')
saveRDS(srt, 'individual/part01.y5sa_dropseq_all.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 6: 2022_JCardiovascDev_ESmall ST data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LoadST <- function(name, sample){
        srt <- Load10X_Spatial(
                data.dir = paste0('/Volumes/shire/data/visium/2022_JCardiovascDev_ESmall/matrix_public/', name, '/outs/'),
                filename = "filtered_feature_bc_matrix.h5",
                assay = "ST",
                slice = sample,
                filter.matrix = T)
}
## Add metadata:
D03 <- LoadST('D03', 'P2ResectD3')
D03$Orig_name <- Cells(D03)
D03$Sample <- 'P2ResectD3'
D03$Sample_id <- 'D002S001'
D03$Age <- 'P005'
D03$Condition <- '03 dpR'
D03$Coord_x_slide <- D03@images$P2ResectD3@coordinates$imagecol
D03$Coord_y_slide <- -D03@images$P2ResectD3@coordinates$imagerow
D03 <- RenameCells(D03, new.names = paste0(D03$Sample_id, '_', Cells(D03)))

D07 <- LoadST('D07', 'P2ResectD7')
D07$Orig_name <- Cells(D07)
D07$Sample_id <- 'D002S002'
D07$Sample <- 'P2ResectD7'
D07$Age <- 'P009'
D07$Condition <- '07 dpR'
D07$Coord_x_slide <- D07@images$P2ResectD7@coordinates$imagecol
D07$Coord_y_slide <- -D07@images$P2ResectD7@coordinates$imagerow
D07 <- RenameCells(D07, new.names = paste0(D07$Sample_id, '_', Cells(D07)))

D14 <- LoadST('D14', 'P2ResectD14')
D14$Orig_name <- Cells(D14)
D14$Sample <- 'P2ResectD14'
D14$Sample_id <- 'D002S003'
D14$Age <- 'P016'
D14$Condition <- '14 dpR'
D14$Coord_x_slide <- D14@images$P2ResectD14@coordinates$imagecol
D14$Coord_y_slide <- -D14@images$P2ResectD14@coordinates$imagerow
D14 <- RenameCells(D14, new.names = paste0(D14$Sample_id, '_', Cells(D14)))

D21 <- LoadST('D21', 'P2ResectD21')
D21$Orig_name <- Cells(D21)
D21$Sample <- 'P2ResectD21'
D21$Sample_id <- 'D002S004'
D21$Age <- 'P023'
D21$Condition <- '21 dpR'
D21$Coord_x_slide <- D21@images$P2ResectD21@coordinates$imagecol
D21$Coord_y_slide <- -D21@images$P2ResectD21@coordinates$imagerow
D21 <- RenameCells(D21, new.names = paste0(D21$Sample_id, '_', Cells(D21)))

## SCT normalize each slide
srt.list <- list(D03, D07, D14, D21)
for(i in 1:4){
        srt.list[[i]] <- SCTransform(srt.list[[i]],
                                     assay = 'ST',
                                     method = "glmGamPoi",
                                     vars.to.regress = c('nCount_ST', 'nFeature_ST'),
                                     do.correct.umi = T,
                                     seed.use = 505,
                                     return.only.var.genes = F)
}
## Merge dataset
merged <- merge(srt.list[[1]], srt.list[2:4])
VariableFeatures(merged) <- U(c(VariableFeatures(srt.list[[1]], assay = 'SCT'),
                                VariableFeatures(srt.list[[2]], assay = 'SCT'),
                                VariableFeatures(srt.list[[3]], assay = 'SCT'),
                                VariableFeatures(srt.list[[4]], assay = 'SCT')))

## Edit metadata
merged$Study <- '2022_JCardiovascDev_ESmall'
merged$Age_group <- 'Neonate'
merged$Method <- 'Visium'
merged$Tissue <- 'Whole_heart'
merged$Sex <- 'Unknown'
merged$Condition_group <- 'Resection'
merged$Genotype <- 'WT'
merged$Genotype_long <- 'WT'
merged$orig.ident <- NULL
merged$Sample <- factor(merged$Sample, levels = c('P2ResectD3', 'P2ResectD7', 'P2ResectD14', 'P2ResectD21'))

## Store spatial coordinates as a reduced dimension
merged <- RunPCA(merged, verbose = F, npcs = 2)
merged@reductions$spatial <- merged@reductions$pca ## create a dimension
merged@reductions$spatial@key <- 'Spatial_'
colnames(merged@reductions$spatial@cell.embeddings) <- c('Spatial_1', 'Spatial_2')
merged@reductions$spatial@cell.embeddings[,1] <- merged$Coord_x_slide
merged@reductions$spatial@cell.embeddings[,2] <- merged$Coord_y_slide
merged@misc$spot_scale <- c(1.26, 1.02, 0.84, 0.66)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged, 'individual/part01.2022_JCardiovascDev_ESmall.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 7: AAV-YAP5SA CD45+ scRNA-seq   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')

sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA') & sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')]
Idents(sc.srt) <- 'Cell_state'

## Export meta and matrices for GEO
# sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
# sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
# sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
# sc.srt <- RenameCells(sc.srt, new.names = paste(sc.srt$sample_pub, sc.srt$orig.name, sep = '_'))
# sc.srt$orig.name <- NULL
# sc.srt$sample <- sc.srt$sample_pub
# write.table(GetAssayData(sc.srt, slot = 'count', assay = 'RNA'),
#             file = '~/Desktop/scRNA_integrated.raw_count.csv',
#             quote = F, sep = ',', row.names = T, col.names = T)
# write.table(GetAssayData(sc.srt, slot = 'data', assay = 'SCT'),
#             file = '~/Desktop/scRNA_integrated.normalized_count.csv',
#             quote = F, sep = ',', row.names = T, col.names = T)
# write.table(sc.srt@meta.data,
#             file = '~/Desktop/scRNA_integrated.meta_data.csv',
#             quote = F, sep = ',', row.names = T, col.names = T)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(sc.srt, 'individual/part01.aav_y5sa_cd45_immune.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset 8: AAV-YAP5SA All snMulti-seq (Not used)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA') & sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')]
Idents(sc.srt) <- 'Cell_state'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(sc.srt, 'individual/part01.aav_y5sa_cd45_immune.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

