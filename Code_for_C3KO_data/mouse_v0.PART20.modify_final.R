####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  C3KO Neonatal MI -- Collaboration with Rich G. Li
####  2023-05-22 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- 'PART20_Modify_Final'

Code_dir <- paste0('/Volumes/shire/project/2023_neoc3ko_rli/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject(Machine = 'Rivendell', Ver = Ver, Part = Step, Catagory = 'mouse',
                Project_dir = '2023_neoc3ko_rli', Data_drive = 'shire')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load objects and update metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
scvi <- readRDS('integrated/PART19.all_reductions.srt_dimreducs.rds')
identical(Cells(full.srt), rownames(scvi$scVI@cell.embeddings))
full.srt@reductions$scVI <- scvi$scVI

full.srt$group2 <- revalue(full.srt$group2, replace = c('wt_p2m' = 'WT MI',
                                                        'wt_p2s' = 'WT Sham',
                                                        'c3_p2m' = 'C3KO MI',
                                                        'c3_p2s' = 'C3KO Sham'))
full.srt$group2 <- factor(full.srt$group2, levels = c('WT Sham', 'WT MI', 'C3KO Sham', 'C3KO MI'))
full.srt$group1 <- paste(full.srt$genotype, full.srt$condition, full.srt$replicate, sep = ' ')
full.srt$group1 <- factor(full.srt$group1, levels = c(
        'WT 5dpSham Rep1',
        'WT 5dpSham Rep2',
        'WT 5dpMI Rep1',
        'WT 5dpMI Rep2',
        'C3KO 5dpSham Rep1',
        'C3KO 5dpSham Rep2',
        'C3KO 5dpMI Rep1',
        'C3KO 5dpMI Rep2'
))
full.srt$Cell_type <- revalue(full.srt$Cell_type, replace = c(
        'Cardiomyocyte' = 'CM',
        'Endocardium' = 'EndoC',
        'Epicardium' = 'CF',
        'Fibroblast' = 'CF',
        'T cell' = 'TC',
        'B cell' = 'BC',
        'Glial' = 'Glial',
        'Pericyte' = 'PC',
        'Adipocyte' = 'Adipo',
        'BEC' = 'EC',
        'Contamination' = 'Ambiguous'
))

full.srt$Cell_type <- factor(full.srt$Cell_type, levels = c(
        'CM', 'CF', 'EC', 'EndoC', 'LEC', 'PC', 'SMC', 'Myeloid', 'TC', 'BC', 'Glial', 'Adipo',
        'Ambiguous', 'Doublet'))

full.srt$Cell_state <- revalue(full.srt$Cell_state, replace = c(
        'EC1' = 'EpiC1',
        'EC2' = 'EpiC2',
        'EC3' = 'EpiC3',
        'EC4' = 'EpiC4',
        'Endocardium' = 'EndoC',
        'Contamination' = 'Ambiguous',
        'BEC' = 'EC',
        'Granu' = 'Mono',
        'MP1' = 'RMP',
        'MP3' = 'Prol. MP',
        'MP2' = 'Spp1+ MP',
        'T cell' = 'TC',
        'B cell' = 'BC'
))
Table(full.srt$Cell_state, full.srt$Cell_type)
full.srt$Cell_state <- as.vector(full.srt$Cell_state)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re Embed CM  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cm.srt <- full.srt[, full.srt$Cell_type == 'CM']
DimPlot2(cm.srt, reduction = 'sub_clean_umap', label = T, cols = mycol_10, group.by = 'Cell_state')

cm.srt$Cell_state <- revalue(cm.srt$Cell_state, replace = c(
        CM2 = 'Prol. CM',
        CM1 = 'Vent. CM1',
        CM3 = 'Vent. CM2',
        CM5 = 'Vent. CM3',
        CM7 = 'Atr. CM',
        CM6 = 'Ambiguous',
        CM4 = 'Doublet'
))
cm.srt$Cell_state <- factor(cm.srt$Cell_state,
                            levels = c('Vent. CM1', 'Vent. CM2', 'Vent. CM3',
                                       'Atr. CM', 'Prol. CM', 'Ambiguous', 'Doublet'))
cm.srt$Non_ambiguous[cm.srt$Cell_state %in% c('Ambiguous', 'Doublet')] <- F
cm.srt2 <- cm.srt[, cm.srt$Non_ambiguous]
cm.srt2 <- RunUMAP(cm.srt2, reduction = 'scVI', dims = 1:50,
                  reduction.name = 'sub_clean_umap', reduction.key = 'subcUMAP_')
DimPlot2(cm.srt2, reduction = 'sub_clean_umap', label = T, cols = mycol_10, group.by = 'Cell_state')

## Return meta and DR to full.srt
full.srt$Cell_state <- revalue(full.srt$Cell_state, replace = c(
        CM2 = 'Prol. CM',
        CM1 = 'Vent. CM1',
        CM3 = 'Vent. CM2',
        CM5 = 'Vent. CM3',
        CM7 = 'Atr. CM',
        CM6 = 'Ambiguous',
        CM4 = 'Doublet'
))
full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cm.srt), 1] <- NA
full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cm.srt), 2] <- NA
full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cm.srt2),] <-
        cm.srt2@reductions$sub_clean_umap@cell.embeddings

DimPlot2(full.srt[, full.srt$Cell_type == 'CM'], reduction = 'sub_umap',
         label = T, cols = mycol_10, group.by = 'Cell_state')/
        DimPlot2(full.srt[, full.srt$Cell_type == 'CM'], reduction = 'sub_clean_umap',
                 label = T, cols = mycol_10, group.by = 'Cell_state')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re Embed CF  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cf.srt <- full.srt[, full.srt$Cell_type == 'CF']
DimPlot2(cf.srt, reduction = 'sub_clean_umap', label = T, cols = mycol_10, group.by = 'Cell_state')

cf.srt <- RunUMAP(cf.srt, reduction = 'scVI', dims = 1:50, min.dist = 0.6,
                  reduction.name = 'sub_clean_umap', reduction.key = 'subcUMAP_')
DimPlot2(cf.srt, reduction = 'sub_clean_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
cf.srt$C3_Pos <- ifelse(cf.srt@assays$RNA@data['C3',] > 2.5, 'C3+', 'C3-')
FeaturePlot2(cf.srt, reduction = 'sub_clean_umap',
             features = 'C3', split.by = 'group2',
             raster = T, ncol = 2, pt.size = 1, min.cutoff = 'q5', max.cutoff = 'q95')


cf.srt <- FindNeighbors(cf.srt, reduction = 'scVI', dims = 1:50) |> FindClusters(resolution = 0.05)
DimPlot2(cf.srt, reduction = 'sub_clean_umap', label = T, cols = mycol_10)

cf.srt$Cell_state <- revalue(cf.srt$RNA_snn_res.0.05, replace = c(
        '0' = 'C3- CF',
        '1' = 'C3+ CF'
))
cf.srt$Cell_state <- factor(cf.srt$Cell_state, levels = c('C3+ CF', 'C3- CF'))

DimPlot2(cf.srt, reduction = 'sub_clean_umap', label = T, cols = mycol_10, group.by = 'Cell_state')

## Return meta and DR to full.srt

full.srt@meta.data[Cells(cf.srt), 'Cell_state'] <- as.vector(cf.srt$Cell_state)

full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cf.srt), 1] <- NA
full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cf.srt), 2] <- NA
full.srt@reductions$sub_clean_umap@cell.embeddings[Cells(cf.srt),] <-
        cf.srt@reductions$sub_clean_umap@cell.embeddings

DimPlot2(full.srt[, full.srt$Cell_type == 'CF'], reduction = 'sub_umap',
         label = T, cols = mycol_10, group.by = 'Cell_state')/
        DimPlot2(full.srt[, full.srt$Cell_type == 'CF'], reduction = 'sub_clean_umap',
                 label = T, cols = mycol_10, group.by = 'Cell_state')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-embed full data ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt$Cell_state <- factor(full.srt$Cell_state, levels = c(
        'Vent. CM1',
        'Vent. CM2',
        'Vent. CM3',
        'Prol. CM',
        'Atr. CM',
        'C3+ CF',
        'C3- CF',
        'EC',
        'EndoC',
        'LEC',
        'Pericyte',
        'SMC',
        'Mono',
        'RMP',
        'Prol. MP',
        'Spp1+ MP',
        'TC',
        'BC',
        'Glial',
        'Adipocyte',
        'Ambiguous',
        'Doublet'
        ))
full.srt$Cell_type[full.srt$Cell_state == 'Ambiguous'] <- 'Ambiguous'
full.srt$Cell_type[full.srt$Cell_state == 'Doublet'] <- 'Doublet'
full.srt$Non_ambiguous[full.srt$Cell_state %in% c('Ambiguous', 'Doublet')] <- F

clean.srt <- full.srt[, full.srt$Non_ambiguous]
clean.srt <- RunUMAP(clean.srt, reduction = 'scVI', dims = 1:50,
                     reduction.name = 'clean_umap', reduction.key = 'cUMAP_')
full.srt@reductions$clean_umap@cell.embeddings[, 1:2] <- NA
full.srt@reductions$clean_umap@cell.embeddings[Cells(clean.srt), 1] <-
        clean.srt@reductions$clean_umap@cell.embeddings[, 1]
full.srt@reductions$clean_umap@cell.embeddings[Cells(clean.srt), 2] <-
        clean.srt@reductions$clean_umap@cell.embeddings[, 2]
full.srt@reductions$clean_umap@key <- 'cUMAP_'
colnames(full.srt@reductions$clean_umap@cell.embeddings) <- c('cUMAP_1', 'cUMAP_2')

DimPlot2(full.srt, reduction = 'clean_umap', group.by = 'Cell_type', cols = mycol_20, raster = T) /
        DimPlot2(full.srt, reduction = 'full_umap', cols = mycol_20, raster = T)

full.srt$Cell_type <- revalue(full.srt$Cell_type, replace = c(
        'EndoC' = 'EC',
        'LEC' = 'EC',
        'Myeloid' = 'Mye',
        'TC' = 'Lym',
        'BC' = 'Lym'
))
Idents(full.srt) <- 'Cell_type'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Verify sex  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt <- AddModuleScore2(full.srt,
                            features = list(chrX_genes$chrX_only_genes, chrY_genes$chrY_only_genes),
                            names = c('ChrX.Score', 'ChrY.Score'),
                            return_z = T)
BoxPlot(full.srt, feature = 'ChrY.Score', group.by = 'group1')
full.srt$sex[full.srt$group1 %in% c('WT 5dpSham Rep2', 'WT 5dpMI Rep2', 'C3KO 5dpMI Rep2')] <- 'F'
full.srt$sex[! full.srt$group1 %in% c('WT 5dpSham Rep2', 'WT 5dpMI Rep2', 'C3KO 5dpMI Rep2')] <- 'M'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full.srt, 'integrated/PART20.annotated_mod.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
