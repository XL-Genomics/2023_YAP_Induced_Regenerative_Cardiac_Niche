####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate directories  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- '03_ligand_receptor'
Project <- '2022_yap5sa_st_rli'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'mouse', Project, 'ithil')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  YAP5SA ST C3ar1 Ligand Receptor  NicheNet Manual Filter ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st.srt <- readRDS('individual/part01.YAP5SA_st.srt.rds')

# sn.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/04.scmulti.wnn.srt.rds')
cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
Idents(cm.srt) <- 'Experiment'

sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
Idents(sc.srt) <- 'sample_pub'
mac.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac')]

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Find DEGs high in Yap5sa cm as possible receptors
sn_cm_deg <- FindMarkers(cm.srt, ident.1 = 'YAP5SA', ident.2 = 'Ctrl', logfc.threshold = 0.2, min.pct = 0.05)
sn_cm_deg[c('Igf1r', 'Tnfrsf12a', 'Itgb1'), ]

rgene <- rownames(sn_cm_deg)[sn_cm_deg$avg_log2FC >= 0.2 & sn_cm_deg$p_val < 0.01] # 3037
sn_cm_deg <- sn_cm_deg[order(sn_cm_deg$avg_log2FC, decreasing = T), ]
sn_cm_deg <- sn_cm_deg[sn_cm_deg$avg_log2FC >= 0.2 & sn_cm_deg$p_val < 0.01, ]
sn_cm_deg$gene <- rownames(sn_cm_deg)
WriteCSV(sn_cm_deg, title = 'part03.1.snmulti_cm_deg_yap5sa_vs_wt')

## Find DEGs high in Yap5sa mac as possible ligands
sc_mac_deg <- FindMarkers(mac.srt, ident.1 = 'YAP5SA', ident.2 = 'WT')
sc_mac_deg[c('Igf1', 'Tnfsf12', 'Adam15'), ]
lgene <- rownames(sc_mac_deg)[sc_mac_deg$avg_log2FC >= 0.25 & sc_mac_deg$p_val < 0.01] # 1239
sc_mac_deg <- sc_mac_deg[order(sc_mac_deg$avg_log2FC, decreasing = T), ]
sc_mac_deg <- sc_mac_deg[sc_mac_deg$avg_log2FC >= 0.25 & sc_mac_deg$p_val < 0.01, ]
sc_mac_deg$gene <- rownames(sc_mac_deg)
WriteCSV(sc_mac_deg, title = 'part03.1.cd45_sc_mac_deg_yap5sa_vs_wt')

## Prepare LT database
# library('nichenetr')
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
# rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
# ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)),
#                                                   !is.na(colnames(ligand_target_matrix))]
# ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
# saveRDS(ligand_target_matrix, 'analysis/part03.ligand_target_matrix.rds')
ligand_target_matrix <- readRDS('analysis/part03.ligand_target_matrix.rds')

## Prepare LR database
# lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
# lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
# head(lr_network)
# Table(lr_network$bonafide)
# lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(ligand),
#                                    to = convert_human_to_mouse_symbols(receptor)) %>% drop_na()
# saveRDS(lr_network, 'analysis/part03.lr_network.rds')
lr_network <- readRDS('analysis/part03.lr_network.rds')
lr_network <- lr_network[lr_network$bonafide, c('from', 'to')]

## Filter for candidate LR pairs and LR network
lr_network_cand <- lr_network[lr_network$from %in% lgene  &  lr_network$to %in% rgene, ] # 26 pairs
p1 <- FeaturePlot2(mac.srt, features = U(lr_network_cand$from), reduction = 'sub_hmn_umap')
p2 <- DotPlot2(mac.srt, features = U(lr_network_cand$from), group.by = 'Cell_state', cols = 'RdBu')
p3 <- DotPlot2(mac.srt, features = U(lr_network_cand$from), group.by = 'Cell_state', split.by = 'sample_pub', cols = 'RdBu')

PlotPDF('01.1.feat.ligand_expr_on_mac', 12, 9)
p1
dev.off()
PlotPDF('01.2.dot.ligand_expr_on_mac', 15, 10)
p2 | p3
dev.off()

lr_network_cand <- lr_network_cand[lr_network_cand$from %in% c('Tnfsf12', 'Igf1', 'Adam15'), ] # 5 pairs
p1 <- FeaturePlot2(cm.srt, features = U(lr_network_cand$to), reduction = 'umap')
p2 <- DotPlot2(cm.srt, features = U(lr_network_cand$to), group.by = 'CM_State', cols = 'RdBu')
p3 <- DotPlot2(cm.srt, features = U(lr_network_cand$to), group.by = 'CM_State', split.by = 'Experiment', cols = 'RdBu')

PlotPDF('01.3.feat.ligand_expr_on_cm', 9, 6)
p1
dev.off()
PlotPDF('01.4.dot.ligand_expr_on_cm', 10, 5)
p2 | p3
dev.off()

lr_network_cand <- lr_network_cand[lr_network_cand$to %in% c('Tnfrsf12a', 'Igf1r', 'Itgb1'), ] # 3 pairs

## Filter LR network
lt_network_cand <- ligand_target_matrix[, lr_network_cand$from]
lt_network_cand_bnr <- lt_network_cand
target_list <- list()
for(i in 1:ncol(lt_network_cand)){
        lt_network_cand_bnr[, i] <- scale(lt_network_cand[, i]) >= 4
        target_list[[i]] <- rownames(lt_network_cand)[lt_network_cand_bnr[, i] == 1]
}
names(target_list) <- colnames(lt_network_cand)
target_list <- target_list[! duplicated(target_list)]

cm.srt <- AddModuleScore2(cm.srt, features = target_list, names = paste0('Score_', names(target_list)))
p1 <- VlnPlot2(mac.srt, features = c('Tnfsf12', 'Igf1', 'Adam15'), group.by = 'sample_pub', pt.size = 0.1)
p2 <- VlnPlot2(cm.srt, features = c('Tnfrsf12a', 'Igf1r', 'Itgb1'), group.by = 'Experiment', pt.size = 0.1)
p3 <- VlnPlot2(cm.srt, features = paste0('Score_', names(target_list)), group.by = 'CM_State', pt.size = 0.1)
PlotPDF('01.5.vln.ligand_receptor_target_expr', 10, 10)
p1/p2/p3
dev.off()

saveRDS(list(rgene, lgene, lr_network_cand, target_list), 'analysis/part03.list_of_nichenet_results.rds')

## Calculate target gene expression scores on ST
## Compute CM1/CM2/Mac scores
st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker,
                          assay = 'SCT', names = paste0('Yap5sa_score_', c('CM1', 'CM2', 'M2', 'C3CF')))

## Compute LRT scores
st.srt <- AddModuleScore2(st.srt, features = as.list(U(lr_network_cand$from)),
                          assay = 'SCT', names = paste0('Ligand_score_', U(lr_network_cand$from)))
st.srt <- AddModuleScore2(st.srt, features = as.list(U(lr_network_cand$to)),
                          assay = 'SCT', names = paste0('Receptor_score_', U(lr_network_cand$to)))
st.srt <- AddModuleScore2(st.srt, features = target_list,
                          assay = 'SCT', names = paste0('Target_score_', names(target_list)))
saveRDS(st.srt, 'analysis/part03.yap5sa_st_with_scores.srt.rds')

p <- FeaturePlotST(st.srt, features = c(
        paste0('Yap5sa_score_', c('CM1', 'CM2', 'M2')),
        paste0('Target_score_', names(target_list))),
        ncol = 4, pt.sizes = st.srt@misc$spot_scale*0.7,
        minvals = c(-0.5, -1,  -1,  0,   -0.02, 0.1, 0),
        maxvals = c(1.5,  0.5, 0.5 , 0.2, 0.05,  0.3, 0.1))
PlotPDF('02.01.feat_st.target_score', 6, 18)
p
dev.off()

p <- FeaturePlotST(st.srt, features = U(lr_network_cand$to),
                   ncol = 2, pt.sizes = st.srt@misc$spot_scale*0.7,
                   minvals = rep(0, 5),
                   maxvals = rep(2, 5))
PlotPDF('02.02.feat_st.recptor_expr', 6, 9)
p
dev.off()

p <- FeaturePlotST(st.srt, features = U(lr_network_cand$from),
                   ncol = 2, pt.sizes = st.srt@misc$spot_scale*0.7,
                   minvals = rep(0, 5),
                   maxvals = rep(2, 5))
PlotPDF('02.03.feat_st.ligand_expr', 6, 9)
p
dev.off()

## Compute Mac-L colocalization
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_M2', paste0('Ligand_score_', lig)))
        y <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_M2', paste0('Ligand_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_M2_', lig), paste0('Colocal_M2_',  lig))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}
plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_M2_', lig)),
                                    title = c(paste0('M2 vs. ', lig)),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = st.srt@misc$spot_scale*0.5, ncol = 2) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.04.feat_st.m2_lig_colocal', 6, 9)
print(p)
dev.off()

## Compute CM-R-T colocalization
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_CM1',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Target_score_', lig)))
        y <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_CM2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Target_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM1_Target_', lig, '_', rec), paste0('Colocal_CM2_Target_',  lig, '_', rec))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}

plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_CM1_Target_', lig, '_', rec),
                                                 paste0('Colocal_CM2_Target_',  lig, '_', rec)),
                                    title = c(paste0('CM1 vs. ', rec, ' vs. target'), paste0('CM2 vs. ', rec, ' vs. target')),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = st.srt@misc$spot_scale*0.5, ncol = 2) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.05.feat_st.cm_r_t_colocal', 6, 18)
print(p)
dev.off()

for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_CM1',
                                                      'Yap5sa_score_M2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Ligand_score_', lig),
                                                      paste0('Target_score_', lig)))
        y <- GetColocalProb(st.srt, meta_features = c('Yap5sa_score_CM2',
                                                      'Yap5sa_score_M2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Ligand_score_', lig),
                                                      paste0('Target_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM1_M2_', lig, '_', rec, '_target'),
                          paste0('Colocal_CM2_M2_',  lig, '_', rec, '_target'))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}

plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_CM1_M2_', lig, '_', rec, '_target'),
                                                 paste0('Colocal_CM2_M2_',  lig, '_', rec, '_target')),
                                    title = c(paste0('CM1 vs. M2 ', lig, ' vs. ', rec, ' vs. Target'),
                                              paste0('CM2 vs. M2 ', lig, ' vs. ', rec, ' vs. Target')),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = st.srt@misc$spot_scale*0.5, ncol = 2) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.06.feat_st.cm_m2_l_r_t_colocal', 6, 18)
print(p)
dev.off()

## Plot for Olson data
eo.srt <- readRDS('individual/part01.2021_NatComm_EOlson.srt.rds')
eo.srt <- AddModuleScore2(eo.srt, features = as.list(U(lr_network_cand$from)),
                          assay = 'SCT', names = paste0('Ligand_score_', U(lr_network_cand$from)))
eo.srt <- AddModuleScore2(eo.srt, features = as.list(U(lr_network_cand$to)),
                          assay = 'SCT', names = paste0('Receptor_score_', U(lr_network_cand$to)))
eo.srt <- AddModuleScore2(eo.srt, features = target_list, assay = 'SCT', names = paste0('Target_score_',
                                                                                        names(target_list)))
eo.srt <- AddModuleScore2(eo.srt, features = Yap5sa_CM_Mac_CF_marker[1:3], assay = 'SCT',
                          names = paste0('Yap5sa_score_', c('CM1', 'CM2', 'M2')))

eo.srt@meta.data <- eo.srt@meta.data[, 1:31]

for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_CM1',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Target_score_', lig)))
        y <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_CM2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Target_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM1_Target_', rec), paste0('Colocal_CM2_Target_', rec))
        rownames(df) <- Cells(eo.srt)
        eo.srt <- AddMetaData(eo.srt, metadata = df)
}

plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = eo.srt,
                                    features = c(paste0('Colocal_CM1_Target_', rec),
                                                 paste0('Colocal_CM2_Target_', rec)),
                                    title = c(paste0('CM1 vs. ', rec, ' vs. Target'),
                                              paste0('CM2 vs. ', rec, ' vs. Target')),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = eo.srt@misc$spot_scale*0.9, ncol = 4) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.07.feat_st.olson_rt_cm_colocal', 12, 18)
print(p)
dev.off()

for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_M2',
                                                      paste0('Ligand_score_', lig)))
        y <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_M2',
                                                      paste0('Ligand_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_M2_Target_', lig), paste0('Colocal_M2_Target_', lig))
        rownames(df) <- Cells(eo.srt)
        eo.srt <- AddMetaData(eo.srt, metadata = df)
}

plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = eo.srt,
                                    features = paste0('Colocal_M2_Target_', lig),
                                    title = paste0('M2 vs. ', lig),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = eo.srt@misc$spot_scale*0.9, ncol = 4) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.08.feat_st.olson_l_m2_colocal', 12, 9)
print(p)
dev.off()

for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        x <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_CM1',
                                                      'Yap5sa_score_M2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Ligand_score_', lig),
                                                      paste0('Target_score_', lig)))
        y <- GetColocalProb(eo.srt, meta_features = c('Yap5sa_score_CM2',
                                                      'Yap5sa_score_M2',
                                                      paste0('Receptor_score_', rec),
                                                      paste0('Ligand_score_', lig),
                                                      paste0('Target_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM1_M2_Target_', lig, '_', rec), paste0('Colocal_CM2_M2_Target_',  lig, '_', rec))
        rownames(df) <- Cells(eo.srt)
        eo.srt <- AddMetaData(eo.srt, metadata = df)
}

plist <- list()
for(i in 1:nrow(lr_network_cand)){
        lig <- lr_network_cand$from[i]
        rec <- lr_network_cand$to[i]
        plist[[i]] <- FeaturePlotST(srt = eo.srt,
                                    features = c(paste0('Colocal_CM1_M2_Target_', lig, '_', rec),
                                                 paste0('Colocal_CM2_M2_Target_',  lig, '_', rec)),
                                    title = c(paste0('CM1 vs. M2 ', lig, ' vs. ', rec, ' vs. Target'),
                                              paste0('CM2 vs. M2 ', lig, ' vs. ', rec, ' vs. Target')),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1, 1, 1, 1),
                                    pt.sizes = eo.srt@misc$spot_scale*0.9, ncol = 4) &
                scale_color_distiller(palette = 'Blues',  values = c(0, 0.5, 1), direction = 0)
}
p <- wrap_plots(plist, ncol = 1)
PlotPDF('02.09.feat_st.olson_lrt_cm2_m2_colocal', 12, 18)
print(p)
dev.off()

p <- FeaturePlotST(eo.srt, features = c(U(lr_network_cand$from), U(lr_network_cand$to)),
                   ncol = 4, pt.sizes = eo.srt@misc$spot_scale*0.9,
                   minvals = c(0, 0, 0, 0, 0, 1),
                   maxvals = c(2, 2, 2, 2, 2, 4))
PlotPDF('02.10.feat_st.olson_ligand_receptor_expr', 12, 18)
p
dev.off()


## Summarize plots
sc.srt <- AddModuleScore2(sc.srt, features = Yap5sa_CM_Mac_CF_marker[1:3],
                          assay = 'SCT', names = paste0('Yap5sa_score_', c('CM1', 'CM2', 'M2')))
cm.srt <- AddModuleScore2(cm.srt, features = Yap5sa_CM_Mac_CF_marker[1:3],
                          assay = 'SCT', names = paste0('Yap5sa_score_', c('CM1', 'CM2', 'M2')))
cm.srt <- AddModuleScore2(cm.srt, features = target_list, nbin = 20,
                          assay = 'SCT', names = paste0('Target_score_', names(target_list)))

p1 <- DimPlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], group.by = 'Cell_state',
               reduction = 'sub_hmn_umap', cols = mycol_10) + labs(title = 'Myeloid (Cd45+ scRNA-seq)')
p2 <- FeaturePlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], features = 'Yap5sa_score_M2',
                   reduction = 'sub_hmn_umap', min.cutoff = 'q2', max.cutoff = 'q98')  +
        labs(title = 'C3ar1+ Mac score')
p3 <- DimPlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], group.by = 'sample_pub',
               reduction = 'sub_hmn_umap') + labs(title = 'Genotype')
p4 <- DimPlot2(cm.srt, group.by = 'CM_State', reduction = 'umap', cols = mycol_10) +
        labs(title = 'CMs (snRNA-seq)')
p5 <- FeaturePlot2(cm.srt, features = 'Yap5sa_score_CM2',
                   reduction = 'umap', min.cutoff = 'q2', max.cutoff = 'q98') +
        labs(title = 'CM2 score')
p6 <- DimPlot2(cm.srt, group.by = 'Experiment', reduction = 'umap') + labs(title = 'Genotype')

p <- wrap_plots(p1, p4, p2, p5, p3, p6, ncol = 2) & labs(x = '', y = '')

PlotPDF('03.1.combine.umaps', 9, 10)
p
dev.off()

DefaultAssay(mac.srt) <- 'RNA'
DefaultAssay(cm.srt) <- 'RNA'
cm.srt <- AddModuleScore2(cm.srt, features = target_list, nbin = 20,
                          assay = 'SCT', names = paste0('Target_score_cbn_', names(target_list)))
levels(cm.srt) <- c("Ctrl", "YAP5SA")

p1 <- VlnPlot2(mac.srt, features = c('Tnfsf12', 'Igf1', 'Adam15'), pt.size = 0.01, raster = F)
p2 <- VlnPlot2(cm.srt, features = c('Tnfrsf12a', 'Igf1r', 'Itgb1'), pt.size = 0.01, raster = F)
p3 <- VlnPlot2(cm.srt, features = paste0('Target_score_cbn_', names(target_list)), pt.size = 0.01, raster = F)
p4 <- DotPlot(cm.srt, features = paste0('Target_score_cbn_', names(target_list))) +
        scale_size_continuous(limits = c(25, 100))

PlotPDF('03.2.combine.vlns', 15, 15)
p1/p2/p3/p4
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Cell Chat - CM to Mac ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('CellChat')

## Build LR db
CellChatDB.my <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
CellChatDB.my$interaction[2022:2025, ] <- NA
CellChatDB.my$interaction$interaction_name[2022:2025] <- c('ADAM15-ITGB1', 'ADAM15-ITGB3', 'ADAM15-ITGA5', 'ADAM15-ITGA9')
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
CellChatDB.my$interaction$interaction_name_2[2022:2025] <- c('Adam15-Itgb1', 'Adam15-Itgb3', 'Adam15-Itga5', 'Adam15-Itga9')

c('Adam15', 'Itgb1', 'Itgb3', 'Itga5', 'Itga9') %in% CellChatDB$geneInfo$Symbol

DoCellChat <- function(srt, group.by = "celltype"){
        CellChatDB <- CellChatDB.my
        ####  Build CellChat Object  ####
        cch <- createCellChat(object = srt, group.by = group.by)
        cch <- addMeta(cch, meta = srt@meta.data)
        groupSize <- as.numeric(table(cch@idents))
        cch@DB <- CellChatDB
        ####  Preprocessing the expression data for cell-cell communication analysis  ####
        cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
        cch <- identifyOverExpressedGenes(cch)
        cch <- identifyOverExpressedInteractions(cch)
        cch <- projectData(cch, PPI.mouse)
        ####  Inference of cell-cell communication network  ####
        cch <- computeCommunProb(cch)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}

## Prepare seurat data
#sn.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/04.scmulti.wnn.srt.rds')
cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
Idents(cm.srt) <- 'Experiment'
sn_cm_deg <- read.table(paste0(Meta_dir, 'part03.1.snmulti_cm_deg_yap5sa_vs_wt.csv'), header = T, sep = ',')

sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
Idents(sc.srt) <- 'sample_pub'
mac.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac')]
sc_mac_deg <- read.table(paste0(Meta_dir, 'part03.1.cd45_sc_mac_deg_yap5sa_vs_wt.csv'), header = T, sep = ',')

lr_include <- CellChatDB.my$interaction$interaction_name[CellChatDB.my$interaction$ligand %in% sc_mac_deg$gene &
                                                                 CellChatDB.my$interaction$receptor %in% sn_cm_deg$gene]


DefaultAssay(cm.srt) <- 'alra'
cm.min.srt <- DietSeurat(cm.srt, assays = 'alra', dimreducs = 'oriumap')
cm.min.srt <- RenameAssays(cm.min.srt, 'alra'  = 'RNA')
cm.min.srt@reductions$oriumap@assay.used <- 'RNA'
cm.min.srt@reductions$oriumap@key <- 'UMAP_'
colnames(cm.min.srt@reductions$oriumap@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
names(cm.min.srt@reductions) <- 'UMAP'
cm.min.srt$Cell_type <- 'CM'
cm.min.srt$Cell_state <- cm.min.srt$CM_State
cm.min.srt$Sample <- revalue(cm.min.srt$Experiment, c('Ctrl' = 'Control'))

DefaultAssay(mac.srt) <- 'SCT'
mac.srt <- AddModuleScore2(mac.srt, features = list('Csf1r', 'C3ar1'), nbin = 20,
                           names = c('Score_Csf1r', 'Score_C3ar1'),
                           return_z = T)

mac.min.srt <- DietSeurat(mac.srt, assays = 'SCT', dimreducs = 'sub_hmn_umap')
mac.min.srt <- RenameAssays(mac.min.srt, 'SCT'  = 'RNA')
mac.min.srt@reductions$sub_hmn_umap@assay.used <- 'RNA'
mac.min.srt@reductions$sub_hmn_umap@key <- 'UMAP_'
colnames(mac.min.srt@reductions$sub_hmn_umap@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
names(mac.min.srt@reductions) <- 'UMAP'
mac.min.srt$Sample <- revalue(mac.min.srt$sample_pub, c('WT' = 'Control'))
mac.min.srt$Cell_state <- c('C3ar1- Mac', 'C3ar1+ Mac')[(mac.min.srt$Score_C3ar1>1) + 1]
# mac.min.srt <- mac.min.srt[, mac.min.srt$Cell_state %in% c('Mac1_Ccr2', 'Mac2_Csf1r')]
# mac.min.srt$Cell_state <- revalue(mac.min.srt$Cell_state, c('Mac2_Csf1r' = 'C3ar1+ Mac', 'Mac1_Ccr2' = 'C3ar1- Mac'))
Idents(mac.srt) <- 'Cell_state'

merged_alra.srt <- merge(cm.min.srt, mac.min.srt, merge.dr = 'UMAP')

ctrl_alra.srt <- merged_alra.srt[, merged_alra.srt$Sample == 'Control']
y5sa_alra.srt <- merged_alra.srt[, merged_alra.srt$Sample != 'Control']

## Run Cell Chat on imputed data
ctrl_alra.cch <- DoCellChat(ctrl_alra.srt, group.by = 'Cell_state')
y5sa_alra.cch <- DoCellChat(y5sa_alra.srt, group.by = 'Cell_state')

PlotPDF('04.1.chrod.cellchat_all_pathway_cm_imputed_ctrl', 6, 4)
netVisual_chord_gene(ctrl_alra.cch, sources.use = c(1, 2), targets.use = c(3, 4),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Control',
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
PlotPDF('04.2.chrod.cellchat_all_pathway_cm_imputed_y5sa', 6, 4)
netVisual_chord_gene(y5sa_alra.cch, sources.use = c(1, 2), targets.use = c(3, 4),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Yap5sa',
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()

PlotPDF('04.3.chrod.cellchat_cm2_imputed_ctrl', 6, 4)
netVisual_chord_gene(ctrl_alra.cch,
                     sources.use = c(1, 2), targets.use = c(4),
                     scale = T,
                     slot.name = "netP",
                     title.name = 'Control',
                     link.target.prop = F,
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
PlotPDF('04.4.chrod.cellchat_cm2_imputed_y5sa', 6, 4)
netVisual_chord_gene(y5sa_alra.cch,
                     sources.use = c(1, 2), targets.use = 4,
                     scale = T,
                     slot.name = "netP",
                     title.name = 'Yap5sa',
                     link.target.prop = F,
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(ctrl_alra.cch, 'analysis/part03.cellchat_cm_mac.ctrl_alra.cch.rds')
saveRDS(y5sa_alra.cch, 'analysis/part03.cellchat_cm_mac.y5sa_alra.cch.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  YAP5SA ST C3 FB Ligand Receptor  NicheNet Manual Filter ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st.srt <- readRDS('individual/part01.YAP5SA_st.srt.rds')

# sn.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/04.scmulti.wnn.srt.rds')
drop.srt <- readRDS('individual/part01.y5sa_dropseq_all.srt.rds')
Idents(drop.srt) <- 'Experiment'
cm.srt <- drop.srt[, drop.srt$celltype %in% c('CM1', 'CM2')]
fb.srt <- drop.srt[, drop.srt$celltype %in% c('CF')]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Find DEGs high in Yap5sa cm as possible receptors
sn_cm_deg <- FindMarkers(cm.srt, ident.1 = 'YAP5SA', ident.2 = 'Ctrl', logfc.threshold = 0.1, min.pct = 0.05)
sn_cm_deg[c('Igf1r', 'Tnfrsf12a', 'Itgb1'), ]

lgene <- rownames(sn_cm_deg)[sn_cm_deg$avg_log2FC >= 0.1 & sn_cm_deg$p_val < 0.01] # 1265
sn_cm_deg <- sn_cm_deg[order(sn_cm_deg$avg_log2FC, decreasing = T), ]
sn_cm_deg <- sn_cm_deg[sn_cm_deg$avg_log2FC >= 0.1 & sn_cm_deg$p_val < 0.01, ]
sn_cm_deg$gene <- rownames(sn_cm_deg)
# WriteCSV(sn_cm_deg, title = 'part03.1.snmulti_cm_deg_yap5sa_vs_wt')

## Find DEGs high in Yap5sa mac as possible ligands
sc_fb_deg <- FindMarkers(fb.srt, ident.1 = 'YAP5SA', ident.2 = 'Ctrl', logfc.threshold = 0.25, min.pct = 0.05, assay = 'RNA')
sc_fb_deg[c('C3'), ]
rgene <- rownames(sc_fb_deg)[sc_fb_deg$avg_log2FC >= 0.1 & sc_fb_deg$p_val < 0.01] # 489
sc_fb_deg <- sc_fb_deg[order(sc_fb_deg$avg_log2FC, decreasing = T), ]
sc_fb_deg <- sc_fb_deg[sc_fb_deg$avg_log2FC >= 0.1 & sc_fb_deg$p_val < 0.01, ]
sc_fb_deg$gene <- rownames(sc_fb_deg)
# WriteCSV(sc_fb_deg, title = 'part03.1.cd45_sc_mac_deg_yap5sa_vs_wt')

## Prepare LT database
# library('nichenetr')
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
# rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
# ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)),
#                                                   !is.na(colnames(ligand_target_matrix))]
# ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
# saveRDS(ligand_target_matrix, 'analysis/part03.ligand_target_matrix.rds')
ligand_target_matrix <- readRDS('analysis/part03.ligand_target_matrix.rds')

## Prepare LR database
# lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
# lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
# head(lr_network)
# Table(lr_network$bonafide)
# lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(ligand),
#                                    to = convert_human_to_mouse_symbols(receptor)) %>% drop_na()
# saveRDS(lr_network, 'analysis/part03.lr_network.rds')
lr_network <- readRDS('analysis/part03.lr_network.rds')
lr_network <- lr_network[lr_network$bonafide, c('from', 'to')]

## Filter for candidate LR pairs and LR network
lr_network_cand <- lr_network[lr_network$from %in% lgene  &  lr_network$to %in% rgene, ] # 6 pairs
p1 <- FeaturePlot2(fb.srt, features = U(lr_network_cand$from))
p2 <- VlnPlot2(fb.srt, features = U(lr_network_cand$from), assay = 'RNA')
p2

PlotPDF('01.1.feat.ligand_expr_on_mac', 12, 9)
p1
dev.off()
PlotPDF('01.2.dot.ligand_expr_on_mac', 15, 10)
p2 | p3
dev.off()

lr_network_cand <- lr_network_cand[lr_network_cand$from %in% c('Tnfsf12', 'Igf1', 'Adam15'), ] # 5 pairs

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Cell Chat - CM to FB ####
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

DoCellChat <- function(srt, group.by = "celltype"){
        CellChatDB <- CellChatDB.my
        ####  Build CellChat Object  ####
        cch <- createCellChat(object = srt, group.by = group.by)
        cch <- addMeta(cch, meta = srt@meta.data)
        groupSize <- as.numeric(table(cch@idents))
        cch@DB <- CellChatDB
        ####  Preprocessing the expression data for cell-cell communication analysis  ####
        cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
        cch <- identifyOverExpressedGenes(cch)
        cch <- identifyOverExpressedInteractions(cch)
        cch <- projectData(cch, PPI.mouse)
        ####  Inference of cell-cell communication network  ####
        cch <- computeCommunProb(cch)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}

## Prepare seurat data
#sn.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/04.scmulti.wnn.srt.rds')
cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
Idents(cm.srt) <- 'Experiment'
sn_cm_deg <- read.table(paste0(Meta_dir, 'part03.1.snmulti_cm_deg_yap5sa_vs_wt.csv'), header = T, sep = ',')

sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
Idents(sc.srt) <- 'sample_pub'
mac.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac')]
sc_mac_deg <- read.table(paste0(Meta_dir, 'part03.1.cd45_sc_mac_deg_yap5sa_vs_wt.csv'), header = T, sep = ',')

lr_include <- CellChatDB.my$interaction$interaction_name[CellChatDB.my$interaction$ligand %in% sc_mac_deg$gene &
                                                                 CellChatDB.my$interaction$receptor %in% sn_cm_deg$gene]

DefaultAssay(cm.srt) <- 'alra'
cm.min.srt <- DietSeurat(cm.srt, assays = 'alra', dimreducs = 'oriumap')
cm.min.srt <- RenameAssays(cm.min.srt, 'alra'  = 'RNA')
cm.min.srt@reductions$oriumap@assay.used <- 'RNA'
cm.min.srt@reductions$oriumap@key <- 'UMAP_'
colnames(cm.min.srt@reductions$oriumap@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
names(cm.min.srt@reductions) <- 'UMAP'
cm.min.srt$Cell_type <- 'CM'
cm.min.srt$Cell_state <- cm.min.srt$CM_State
cm.min.srt$Sample <- revalue(cm.min.srt$Experiment, c('Ctrl' = 'Control'))

DefaultAssay(mac.srt) <- 'SCT'
mac.srt <- AddModuleScore2(mac.srt, features = list('Csf1r', 'C3ar1'), nbin = 20,
                           names = c('Score_Csf1r', 'Score_C3ar1'),
                           return_z = T)

mac.min.srt <- DietSeurat(mac.srt, assays = 'SCT', dimreducs = 'sub_hmn_umap')
mac.min.srt <- RenameAssays(mac.min.srt, 'SCT'  = 'RNA')
mac.min.srt@reductions$sub_hmn_umap@assay.used <- 'RNA'
mac.min.srt@reductions$sub_hmn_umap@key <- 'UMAP_'
colnames(mac.min.srt@reductions$sub_hmn_umap@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
names(mac.min.srt@reductions) <- 'UMAP'
mac.min.srt$Sample <- revalue(mac.min.srt$sample_pub, c('WT' = 'Control'))
mac.min.srt$Cell_state <- c('C3ar1- Mac', 'C3ar1+ Mac')[(mac.min.srt$Score_C3ar1>1) + 1]
# mac.min.srt <- mac.min.srt[, mac.min.srt$Cell_state %in% c('Mac1_Ccr2', 'Mac2_Csf1r')]
# mac.min.srt$Cell_state <- revalue(mac.min.srt$Cell_state, c('Mac2_Csf1r' = 'C3ar1+ Mac', 'Mac1_Ccr2' = 'C3ar1- Mac'))
Idents(mac.srt) <- 'Cell_state'

merged_alra.srt <- merge(cm.min.srt, mac.min.srt, merge.dr = 'UMAP')

ctrl_alra.srt <- merged_alra.srt[, merged_alra.srt$Sample == 'Control']
y5sa_alra.srt <- merged_alra.srt[, merged_alra.srt$Sample != 'Control']

## Run Cell Chat on imputed data
ctrl_alra.cch <- DoCellChat(ctrl_alra.srt, group.by = 'Cell_state')
y5sa_alra.cch <- DoCellChat(y5sa_alra.srt, group.by = 'Cell_state')

PlotPDF('04.1.chrod.cellchat_all_pathway_cm_imputed_ctrl', 6, 4)
netVisual_chord_gene(ctrl_alra.cch, sources.use = c(1, 2), targets.use = c(3, 4),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Control',
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
PlotPDF('04.2.chrod.cellchat_all_pathway_cm_imputed_y5sa', 6, 4)
netVisual_chord_gene(y5sa_alra.cch, sources.use = c(1, 2), targets.use = c(3, 4),
                     slot.name = "netP",
                     scale = T,
                     title.name = 'Yap5sa',
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()

PlotPDF('04.3.chrod.cellchat_cm2_imputed_ctrl', 6, 4)
netVisual_chord_gene(ctrl_alra.cch,
                     sources.use = c(1, 2), targets.use = c(4),
                     scale = T,
                     slot.name = "netP",
                     title.name = 'Control',
                     link.target.prop = F,
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
PlotPDF('04.4.chrod.cellchat_cm2_imputed_y5sa', 6, 4)
netVisual_chord_gene(y5sa_alra.cch,
                     sources.use = c(1, 2), targets.use = 4,
                     scale = T,
                     slot.name = "netP",
                     title.name = 'Yap5sa',
                     link.target.prop = F,
                     pairLR.use = data.frame(interaction_name = lr_include),
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(ctrl_alra.cch, 'analysis/part03.cellchat_cm_mac.ctrl_alra.cch.rds')
saveRDS(y5sa_alra.cch, 'analysis/part03.cellchat_cm_mac.y5sa_alra.cch.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
