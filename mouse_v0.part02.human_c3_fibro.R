####--------------------------------------------------------------------------------------------------------------------
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate directories  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '0'
Step <- '02_human_c3_fibro'
Project <- '2022_yap5sa_st_rli'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'mouse', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Try reprocess raw (from human compendium)   ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('/Volumes/fangorn/project/2019_scrna_compendium/rdata/human_v1.1/integrated/STEP05.merged.clean.srt.rds')
Idents(srt) <- 'study'
srt <- FindNeighbors(srt, reduction = 'harmony') |> FindClusters(res = 0.1)
DimPlot2(srt, reduction = 'hmn_umap', label = T)

fb <- srt[, srt$CBN_snn_res.0.1 %in% c(2,5)]
fb <- fb[, fb$study %in% c('2020_Circulation_PEllinor', '2020_Nature_STeichmann')]
DimPlot2(fb, reduction = 'hmn_umap', label = T, group.by = 'study')

fb <- ProcessSrt_std(fb, assay = 'CBN', do.umap = F)
fb <- ProcessSrt_hmn(fb, assay = 'CBN', var.toal = 0.7)
DimPlot2(fb, reduction = 'hmn_umap', label = T, group.by = 'study')
fb <- AddModuleScore2(fb, features = Yap5sa_CM_Mac_CF_marker, name = paste0('Score_', names(Yap5sa_CM_Mac_CF_marker)))

fb <- fb[, fb$tissue %in% c('Interventricular_septum', 'Left_ventricle', 'Left_ventricular_apex', 'Right_ventricle')]
p1 <- FeaturePlot2(fb, features = c('Score_C3_FB', 'C3'), reduction = 'hmn_umap',
                   max.cutoff = 'q99', min.cutoff = 'q10')
p2 <- FeaturePlot3(fb, features = c('Score_C3_FB', 'C3'), reduction = 'hmn_umap', adjust = 3)
p3 <- DimPlot2(fb, reduction = 'hmn_umap', group.by = 'study', cols = mycol_10)
p4 <- DimPlot2(fb, reduction = 'hmn_umap', group.by = 'tissue', cols = mycol_10)
p <- wrap_plots(p1[[1]], p1[[2]], p2[[1]], p2[[2]], p3, p4, ncol = 2)
PlotPDF('01.combine_compendium_cf', 10, 15)
p
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(fb, 'analysis/part02.human_merged_compendium.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Try downloaded r data (processed by original authors)   ####
####--------------------------------------------------------------------------------------------------------------------
ori_data1 <- readRDS('/Volumes/shire/data/scrnaseq/2020_Nature_STeichmann/matrix_public/raw.srt.rds')
ori_data2 <- readRDS('/Volumes/shire/data/scrnaseq/2020_Circulation_PEllinor/matrix_public/raw.v4.srt.rds')
fb2 <- merge(ori_data1[, ori_data1$cell_type == 'Fibroblast'],
             ori_data2[, ori_data2$CellType %in% c('01. Fibroblast I', '02. Fibroblast II', '14. Fibroblast III')])
fb2$sample[is.na(fb2$sample)] <- fb2$experiment[is.na(fb2$sample)]
fb2$Cell_states <- fb2$cell_states
fb2$Cell_states[is.na(fb2$Cell_states)] <- fb2$Subcluster[is.na(fb2$Cell_states)]
sort(table(fb2$sample))
fb2 <- fb2[, fb2$sample %in% names(table(fb2$sample))[table(fb2$sample)>=200]]
fb2 <- ProcessSrt_std(fb2, assay = 'RNA', do.umap = F)
fb2 <- ProcessSrt_hmn(fb2, assay = 'RNA', var.toal = 0.9, haromize.by = 'sample')
fb2 <- RunUMAP(fb2, reduction = 'harmony', dims = 1:50, reduction.name = 'hmn_umap', reduction.key = 'hmnumap_',
               n.neighbors = 20, min.dist = 0.5)
DimPlot2(fb2, reduction = 'hmn_umap', group.by = 'Cell_states', cols = mycol_20)
DimPlot2(fb2, reduction = 'hmn_umap', split.by = 'Cell_states', group.by = 'Cell_states', cols = mycol_20, ncol = 3)
FeaturePlot3(fb2, features = c('C3'), reduction = 'hmn_umap', adjust = 4)
FeaturePlot2(fb2, features = c('C3'), reduction = 'hmn_umap', max.cutoff = 'q99', min.cutoff = 'q10')
VlnPlot2(fb2, features = c('C3'), group.by = 'Cell_states')
DotPlot2(fb2, features = c('C3'), group.by = 'Cell_states')

fb2 <- AddModuleScore2(fb2, features = Yap5sa_CM_Mac_CF_marker, name = paste0('Score_', names(Yap5sa_CM_Mac_CF_marker)))

p1 <- FeaturePlot2(fb2, features = c('Score_C3_FB', 'C3'), reduction = 'hmn_umap', max.cutoff = 'q99', min.cutoff = 'q10')
p2 <- FeaturePlot3(fb2, features = c('Score_C3_FB', 'C3'), reduction = 'hmn_umap', adjust = 3)
p3 <- DimPlot2(fb2, reduction = 'hmn_umap', group.by = 'Cell_states', cols = mycol_20)
p4 <- VlnPlot2(fb2, features = c('C3', 'Score_C3_FB'), group.by = 'Cell_states')
p5 <- DotPlot2(fb2, features = c('C3', 'Score_C3_FB'), group.by = 'Cell_states')

p <- wrap_plots(p1[[1]], p1[[2]], p2[[1]], p2[[2]], p3, p4[[1]], p5, p4[[2]], ncol = 2)
PlotPDF('02.combine_uploaded_rdata_cf', 10, 15)
p
dev.off()

## Find DEGs
fb2$C3_high <- factor('Low', levels = c('High', 'Low'))
fb2$C3_high[fb2$Score_C3_FB > 0.25] <- 'High'
DimPlot2(fb2, group.by = 'C3_high')
Idents(fb2) <- 'C3_high'
deg <- FindMarkers(fb2, ident.1 = 'High', ident.2 = 'Low',
                   logfc.threshold = 0.1,
                   test.use = 'MAST')
deg.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.2], human_or_mouse = 'human'),
                 'Low' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC < -0.2], human_or_mouse = 'human'))
degtop.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.5], human_or_mouse = 'human'),
                    'Low' = rownames(deg)[deg$avg_log2FC < -0.25])
p <- MarkerVolcano(deg, label = T, min.log2FC = 0.2, label_genes = unlist(degtop.list), max.overlaps = 40)
PlotPDF('03.vol.c3_cf_deg', 10, 10)
p
dev.off()

enr <- ModuleEnrichment(deg.list, human_or_mouse = 'human')
tab1 <- enr$GO$GO_High
tab1 <- tab1[tab1$p.adjust <= 0.05, ]
tab2 <- enr$Reactome$Reactome_High
tab2 <- tab2[tab2$p.adjust <= 0.1, ]

WriteCSV(tab1, 'part02.1.human_c3_cf_upgene_go_enrich')
WriteCSV(tab2, 'part02.2.human_c3_cf_upgene_pathway_enrich')
####--------------------------------------------------------------------------------------------------------------------
saveRDS(fb2, 'analysis/part02.human_merged_orig.srt.rds')
saveRDS(deg, 'analysis/part02.human_c3_cf_deg.df.rds')
saveRDS(enr, 'analysis/part02.human_c3_cf_upgene_enrich.list.rds')
####--------------------------------------------------------------------------------------------------------------------
