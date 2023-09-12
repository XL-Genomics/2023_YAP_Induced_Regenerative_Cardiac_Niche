####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  C3KO Neonatal MI -- Collaboration with Rich G. Li
####  2023-05-22 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- 'PART24_CellChat'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2023_neoc3ko_rli/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject(Machine = 'Rivendell', Ver = Ver, Part = Step, Catagory = 'mouse',
                Project_dir = '2023_neoc3ko_rli', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load objects and update metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_w_amb.srt <- readRDS('integrated/PART20.annotated_mod.srt.rds')
full_w_amb.srt <- full_w_amb.srt[, full_w_amb.srt$condition != '5dpSham']
full.srt <- full_w_amb.srt[, full_w_amb.srt$Non_ambiguous]
full.srt <- DropMetaLevels(full.srt)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Prepare CellChat database  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~---
suppressMessages(library('CellChat'))

CellChatDB.my <- CellChatDB.mouse
DoCellChat <- function(srt, group.by, assay, secreted_only = F, trim = 0.1){
    CellChatDB <- CellChatDB.my
    if(secreted_only){CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")}
    ####  Build CellChat Object  ####
    cch <- createCellChat(object = srt, group.by = group.by, assay = assay)
    cch <- addMeta(cch, meta = srt@meta.data)
    groupSize <- as.numeric(table(cch@idents))
    cch@DB <- CellChatDB
    ####  Preprocessing the expression data for cell-cell communication analysis  ####
    cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
    cch <- identifyOverExpressedGenes(cch)
    cch <- identifyOverExpressedInteractions(cch)
    cch <- projectData(cch, PPI.human)
    ####  Inference of cell-cell communication network  ####
    cch <- computeCommunProb(cch, type = "truncatedMean", trim = trim)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cch <- filterCommunication(cch)
    # CellChat computes the communication probability on signaling pathway level
    cch <- computeCommunProbPathway(cch)
    cch <- aggregateNet(cch)
    return(cch)
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  CellChat  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
DimPlot2(full.srt, reduction = 'clean_umap')

tmp.srt <- full.srt[, full.srt$Cell_type %in% c('CM', 'CF') & full.srt$group2 %in% c('WT MI', 'C3KO MI')]
Table(tmp.srt$Cell_type, tmp.srt$group1)
DimPlot2(tmp.srt, reduction = 'sub_clean_umap', split.by = 'Cell_type')
tmp.srt <- DropMetaLevels(tmp.srt)

wt.cch <- DoCellChat(tmp.srt[, tmp.srt$group2 == 'WT MI'],
                     group.by = 'Cell_state', assay = 'RNA', secreted_only = T, trim = 0.1)
ko.cch <- DoCellChat(tmp.srt[, tmp.srt$group2 == 'C3KO MI'],
                     group.by = 'Cell_state', assay = 'RNA', secreted_only = T, trim = 0.1)
full.cch <- mergeCellChat(list(WT = wt.cch,
                               KO = ko.cch),
                          add.names = c('WT', 'C3KO'))

saveRDS(wt.cch, 'analysis/PART24.wt_cm_cf_cell_state.cellchat.rds')
saveRDS(ko.cch, 'analysis/PART24.c3ko_cm_cf_cell_state.cellchat.rds')

####  Per-cell subtype pairwise interactions  ####
cst <-  levels(tmp.srt$Cell_state)
mtx <- wt.cch@net$count + ko.cch@net$count
for (i in 1:L(cst)) {
    sender <- cst[i]
    p_list <- list()
    for(j in 1:L(cst)){
        receiver <- cst[j]
        if(mtx[sender, receiver] > 0) {
            p_list[[paste0(sender, "_", receiver)]] <- netVisual_bubble(full.cch,
                                                                        comparison = c(1, 2),
                                                                        sources.use = sender,
                                                                        targets.use = receiver,
                                                                        angle.x = 45,
                                                                        remove.isolate = T) +
                labs(title = paste(sender, "to", receiver))
        }
    }
    p <- wrap_plots(p_list, ncol = 3)
    PlotPDF(paste0('1.', str_pad(string = i, width = 2, pad = '0'),
                   '.dot.cbn_cellchat_', sender, '_to_each_cellsubtype'), 12, 30)
    print(p)
    dev.off()
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Direct filter from DEGs  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
lr_network <- readRDS('../../../2022_yap5sa_st_rli/rdata/mouse_v0/analysis/part03.lr_network.rds')
lr_network <- lr_network[lr_network$bonafide, c('from', 'to')]


c3_cf <- tmp.srt[, tmp.srt$Cell_state == 'C3+ CF']
c3_cf <- DropMetaLevels(c3_cf)
Idents(c3_cf) <- 'group2'
c3cf_deg <- FindMarkers(c3_cf, ident.1 = 'C3KO MI', ident.2 = 'WT MI')
c3cf_deg <- c3cf_deg[ c3cf_deg$p_val_adj < 0.01, ]
c3cf_up <- rownames(c3cf_deg)[c3cf_deg$avg_log2FC > 0]
c3cf_dn <- rownames(c3cf_deg)[c3cf_deg$avg_log2FC < 0]

prol_cm <- tmp.srt[, tmp.srt$Cell_state == 'Prol. CM']
prol_cm <- DropMetaLevels(prol_cm)
Idents(prol_cm) <- 'group2'
prol_cm_deg <- FindMarkers(prol_cm, ident.1 = 'C3KO MI', ident.2 = 'WT MI')
prol_cm_deg <- prol_cm_deg[ prol_cm_deg$p_val_adj < 0.01, ]
prol_cm_up <- rownames(prol_cm_deg)[prol_cm_deg$avg_log2FC > 0]
prol_cm_dn <- rownames(prol_cm_deg)[prol_cm_deg$avg_log2FC < 0]



## Filter for candidate LR pairs and LR network
up_cand <- lr_network[lr_network$from %in% c3cf_up  &  lr_network$to %in% prol_cm_up, ] # 2 pairs
dn_cand <- lr_network[lr_network$from %in% c3cf_dn  &  lr_network$to %in% prol_cm_dn, ] # 1 pairs


cm <- tmp.srt[, tmp.srt$Cell_type == 'CM']
cm <- DropMetaLevels(cm)
Idents(cm) <- 'group2'
cm_deg <- FindMarkers(cm, ident.1 = 'C3KO MI', ident.2 = 'WT MI')
cm_deg <- cm_deg[ cm_deg$p_val_adj < 0.01, ]
cm_up <- rownames(cm_deg)[cm_deg$avg_log2FC > 0]
cm_dn <- rownames(cm_deg)[cm_deg$avg_log2FC < 0]

## Filter for candidate LR pairs and LR network
up_cand2 <- lr_network[lr_network$from %in% c3cf_up  &  lr_network$to %in% cm_up, ] # 1 pairs
dn_cand2 <- lr_network[lr_network$from %in% c3cf_dn  &  lr_network$to %in% cm_dn, ] # 1 pairs

PlotPDF('2.1.direct_filter_from_DEGs', 6, 6)
VlnPlot2(tmp.srt, features = c('Thbs1', 'Cd47', 'Fgf1', 'Fgfr2'), 
         group.by = 'Cell_state', split.by = 'group2', ncol = 2) + 
    RestoreLegend()
DotPlot2(tmp.srt, features = c('Thbs1', 'Cd47', 'Fgf1', 'Fgfr2'), 
         group.by = 'Cell_state', split.by = 'group2', cols = 'RdBu') 
dev.off()
