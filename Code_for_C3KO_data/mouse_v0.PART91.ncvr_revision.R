####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  C3KO Neonatal MI -- Collaboration with Rich G. Li
####  2023-05-22 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- 'PART91_NCVR_Revision'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2023_neoc3ko_rli/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject(Machine = 'Gondor', Ver = Ver, Part = Step, Catagory = 'mouse',
                Project_dir = '2023_neoc3ko_rli', Data_drive = 'bree')

library(DESeq2)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_w_amb.srt <- readRDS('integrated/PART20.annotated_mod.srt.rds')
full.srt <- full_w_amb.srt[, full_w_amb.srt$Non_ambiguous]
full.srt$Cell_state <- droplevels(full.srt$Cell_state)
full.srt$Cell_type <- droplevels(full.srt$Cell_type)
Color_group2 <- c(mycol_14[c(2, 9)], mycol_10[2], mycol_14[10])
cf.srt <- full.srt[, full.srt$Cell_type == 'CF']
cm.srt <- full.srt[, full.srt$Cell_type == 'CM']
drop.srt <- readRDS(paste0('../../../2022_yap5sa_st_rli/rdata/mouse_v0/external/fransico_data/',
                           'reference_scrnaseq/Reference_combined_compendium_&_yap5sa_CMs_v2.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Compute DEG  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ct <- c('CM', 'CF')
srt_list <- SplitObject(full.srt, 'Cell_type')[c('CM', 'CF')]

sc_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        Idents(srt_list[[i]]) <- 'group2'
        deg <- FindMarkers(srt_list[[i]], ident.1 = 'C3KO MI', ident.2 = 'WT MI')
        deg <- deg[deg$p_val_adj < 0.05, ]
        deg$gene <- rownames(deg)
        sc_deg_list[[i]] <- deg
        names(sc_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(sc_deg_list, 'analysis/PART91.sc_deg_c3komi_vs_wtmi.srt_mk.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

psb_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        srt_list[[i]] <- srt_list[[i]][, srt_list[[i]]$group2 %in% c('C3KO MI', 'WT MI')]
        psb <- AggregateExpression(srt_list[[i]], assays = 'RNA', slot = 'count', group.by = 'group1')$RNA
        meta <- srt_list[[i]]@meta.data[!duplicated(srt_list[[i]]$group1), ]
        rownames(meta) <- meta$group1
        meta <- meta[colnames(psb), ]
        dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group2)
        dds <- estimateSizeFactors(dds)
        dds <- DESeq(dds)
        res <- results(dds, alpha = 0.05, independentFiltering = T, contrast = c('group2', 'C3KO MI', 'WT MI'))
        res <- lfcShrink(dds, contrast = c('group2', 'C3KO MI', 'WT MI'), res = res, type = "ashr")
        deg_up <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 0]
        deg_dn <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < 0]
        deg_up <- deg_up[!is.na(deg_up)]
        deg_dn <- deg_dn[!is.na(deg_dn)]
        psb_deg_list[[i]] <- list('UP' = deg_up, 'DN' = deg_dn)
        names(psb_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(psb_deg_list, 'analysis/PART91.psb_deg_c3komi_vs_wtmi.list.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

final_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        sc_deg <- sc_deg_list[[i]]
        final_deg_list[[i]] <- sc_deg[(sc_deg$gene %in% psb_deg_list[[i]][['UP']] & sc_deg$avg_log2FC > 0) |
                                              (sc_deg$gene %in% psb_deg_list[[i]][['DN']] & sc_deg$avg_log2FC < 0), ]
        gl <- split(final_deg_list[[i]]$gene, final_deg_list[[i]]$avg_log2FC > 0)
        print(str(gl))
        names(final_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(final_deg_list, 'analysis/PART91.psb_deg_c3komi_vs_wtmi.srt_mk.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Plots  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot 1 ####
p1 <- DimPlot2(full.srt, group.by = 'Cell_type', reduction = 'clean_umap', cols = mycol_10, pt.size = 1, raster = T) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'Integrated snRNA-seq')
p1
PlotPDF('1.01.umap_global_cell_type', 5, 5)
p1
dev.off()


#### Plot 2 ####
tmp.srt <- full.srt[, unlist(DownsampleByMeta(full.srt, meta_var = 'group2', random = T, down_to_min_group = T))]
p2 <- DimPlot2(tmp.srt, group.by = 'Cell_type', reduction = 'clean_umap',
               cols = mycol_10, pt.size = 2, raster = T, split.by = 'group2', ncol = 4) &
        labs(x = 'UMAP1', y = 'UMAP2')
p2
PlotPDF('1.02.umap_global_cell_type_split_by_exprm', 12, 4)
p2
dev.off()


#### Plot 3 ####
p3 <- DimPlot2(cf.srt, group.by = 'Cell_state', reduction = 'sub_clean_umap',
               cols = mycol_10[2:1], pt.size = 0.5, raster = F) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'Cardiac Fibroblast')
p3
PlotPDF('1.03.umap_cf_cell_state', 5, 5)
p3
dev.off()


#### Plot 4 ####
p4 <- FeaturePlot2(cf.srt, features = 'C3', reduction = 'sub_clean_umap', pt.size = 0.5, raster = F) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'C3 Expression')
p4
PlotPDF('1.04.feat_c3_in_cf', 5, 5)
p4
dev.off()


#### Plot 5 ####
tmp.srt <- cf.srt[, unlist(DownsampleByMeta(cf.srt, meta_var = 'group2', random = T, down_to_min_group = T))]
p5 <- FeaturePlot2(tmp.srt, features = 'C3', reduction = 'sub_clean_umap', pt.size = 3, raster = T,
                   split.by = 'group2', ncol = 2) &
        labs(x = 'UMAP1', y = 'UMAP2')
p5
PlotPDF('1.05.feat_c3_cf_split_by_exprm', 12, 4)
p5
dev.off()


#### Plot 6 ####
c3cf.srt <- cf.srt[, cf.srt$Cell_state == 'C3+ CF' & cf.srt$group2 %in% c('WT MI', 'C3KO MI')]
Idents(c3cf.srt) <- 'group2'
deg <- FindMarkers(c3cf.srt, ident.1 = 'C3KO MI', ident.2 = 'WT MI')
deg <- deg[deg$p_val_adj < 0.05,]
deg <- deg[order(deg$avg_log2FC, decreasing = T),]
deg$gene <- rownames(deg)
p6 <- MarkerVolcano(deg, label_genes = c(head(deg$gene, 8), tail(deg$gene, 8))) +
        scale_x_continuous(limits = c(-5, 2)) +
        theme_classic() +
        theme_Publication() +
        theme(aspect.ratio = 1)
p6
PlotPDF('1.06.volcano_c3ko_mi_vs_wt_mi', 6, 6)
p6
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(deg, 'analysis/PART91.deg_cf_c3mi_vs_wtmi.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### Plot 7 ####
deg <- readRDS('analysis/PART91.deg_cf_c3mi_vs_wtmi.rds')
gl <- split(rownames(deg), deg$avg_log2FC > 0)
names(gl) <- c('Down in C3KO', 'Up in C3KO MI')
str(gl)
enrich <- ModuleEnrichment(module_list = gl, human_or_mouse = 'mouse')
go <- enrich$GO$`GO_Up in C3KO MI`
go <- go[go$p.adjust < 0.05, ]
WriteCSV(go, 'PART91.cf_c3mi_vs_wtmi_go_bp')
react <- enrich$Reactome$`Reactome_Up in C3KO MI`
react <- react[react$p.adjust < 0.05, ]
WriteCSV(react, 'PART91.cf_c3mi_vs_wtmi_reactome')
kegg <- enrich$KEGG$`KEGG_Up in C3KO MI`
kegg <- kegg[kegg$p.adjust < 0.05, ]
WriteCSV(kegg, 'PART91.cf_c3mi_vs_wtmi_kegg')

go_oi <- c('ATP biosynthetic process',
           'oxidative phosphorylation',
           'cellular respiration',
           'ATP metabolic process',
           'extracellular matrix organization',
           'actin-mediated cell contraction',
           'cardiac muscle contraction'
           )
go_plot <- go[go$Description %in% go_oi, ]
go_plot$Description <- factor(go_plot$Description, levels = go_plot$Description[order(go_plot$p.adjust)])
p7 <- ggplot(go_plot, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.3) +
        geom_text(aes(label = str_to_sentence(Description), x = 0.1, y = Description), hjust = 'left') +
        labs(x = '-Log10 adjusted p value', y = 'Top Enriched Terms', title = 'GO Biological Processes') +
        scale_y_discrete(limit = rev) &
        theme_classic() &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p7
PlotPDF('1.07.barplot.c3ko_up_enriched_bp', 6, 3)
p7
dev.off()


#### Plot 8 ####
p8 <- DimPlot2(cm.srt, group.by = 'Cell_state', reduction = 'sub_clean_umap',
               cols = mycol_10, pt.size = 2, raster = T) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'Cardiomyocyte')
p8
PlotPDF('1.08.umap_cm_cell_state', 5, 5)
p8
dev.off()


#### Plot 9 ####
cm.srt <- full.srt[, full.srt$Cell_type == 'CM']
Idents(cm.srt) <- 'Cell_state'
mk <- FindAllMarkers(cm.srt, only.pos = T)
p9 <- DotPlot2(cm.srt, features = c('Rbfox1', 'Mpped2', 'Scn5a',
                                    'Xirp2', 'Nppb', 'Ankrd1',
                                    'Ddc', '1700042O10Rik', 'Fignl1',
                                    'Top2a', 'Mki67', 'Pcna',
                                    'Ryr3', 'Myl4', 'Fgf12'), col.min = 0)
p9
PlotPDF('1.09.dot_cm_cell_state_marker', 5, 5)
p9
dev.off()


#### Plot 10 ####
library('DESeq2')
tmp.srt <- NormalizeData(full.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'genotype', 'condition', 'replicate')]
tmp.srt <- FindVariableFeatures(tmp.srt, nfeatures = 5e3)

psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'group1',
                           features = VariableFeatures(tmp.srt))$RNA
meta <- tmp.srt@meta.data[!duplicated(tmp.srt$group1), ]
rownames(meta) <- meta$group1
meta <- meta[colnames(psb), ]
meta$age <- as.numeric(as.vector(meta$age))
dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group1)
rld <- rlog(dds, blind = T)
data <- plotPCA(rld, intgroup = "sex", ntop = 5000) +
        scale_color_manual(values = c(mycol_10))
data <- data$data
df2 <- U(full.srt@meta.data[, c('sex', 'group1', 'genotype')])
df2 <- df2[order(df2$group1), ]
data$Sex <- df2$sex
data$Genotype <- df2$genotype
p10.1 <- ggplot(data) +
        geom_point(aes(x = PC1, y = PC2, color = Genotype, shape = Sex), size = 3) +
        scale_color_manual(values = c('red3', 'black')) +
        labs(color = '', x = '', y = '') +
        labs(x = 'PC1 30% variance', y = 'PC2 24% variance', color = 'Genotype') +
        theme_classic() +
        theme(aspect.ratio = 1,
              axis.line = element_line(size = 0.5, color = 'black'),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.background = element_blank(),
              panel.background = element_blank())
data <- PlotPCA_DESeq2(rld, intgroup = "sex", ntop = 5000, dims = c(2, 3)) +
        scale_color_manual(values = c(mycol_10))
data <- data$data
df2 <- U(full.srt@meta.data[, c('sex', 'group1', 'condition')])
df2 <- df2[order(df2$group1), ]
data$Sex <- df2$sex
data$Condition <- df2$condition
p10.2 <- ggplot(data) +
        geom_point(aes(x = PC1, y = PC2, color = Condition, shape = Sex), size = 3) +
        scale_color_manual(values = c('red3', 'black')) +
        labs(color = '', x = '', y = '') +
        labs(x = 'PC2 24% variance', y = 'PC3 18% variance', color = 'Condition') +
        theme_classic() +
        theme(aspect.ratio = 1,
              axis.line = element_line(size = 0.5, color = 'black'),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.background = element_blank(),
              panel.background = element_blank())
p10.1 | p10.2
PlotPDF('1.10.pc.all_cell_global_pc_by_sex_genotype', 5, 5)
p10.1 | p10.2
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p10$data, 'analysis/PART91.all_cell_pseudobulk_pca.df.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Plot 10 v2 ####
tmp.srt <- NormalizeData(full.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'genotype', 'condition', 'replicate')]
tmp.srt <- FindVariableFeatures(tmp.srt, nfeatures = 5e3)
psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'group1',
                           features = VariableFeatures(tmp.srt), return.seurat = T)
meta <- U(tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'genotype', 'condition')])
rownames(meta) <- meta$group1
psb <- AddMetaData(psb, metadata = meta[Cells(psb), ])
psb <- FindVariableFeatures(psb)
psb@assays$RNA@var.features <- tmp.srt@assays$RNA@var.features
psb <- RunPCA(psb, npcs = 3)
p10.1 <- DimPlot2(psb, reduction = 'pca', group.by ='sex', pt.size = 3) +
        labs(x = 'PC1', y = 'PC2', title = 'Sex')
p10.2 <- DimPlot2(psb, reduction = 'pca', group.by ='genotype', pt.size = 3) +
        labs(x = 'PC1', y = 'PC2', title = 'Genotype')
p10.3 <- DimPlot2(psb, reduction = 'pca', group.by ='condition', dims = 2:3, pt.size = 3) +
        labs(x = 'PC2', y = 'PC3', title = 'Condition')
p10 <- wrap_plots(p10.1, p10.2, p10.3, ncol = 3) &
        scale_color_manual(values = c('black', 'red'))
p10
PlotPDF('1.10.2.pc.all_cell_global_pc_by_sex_genotype', 12, 4)
p10
dev.off()


#### Plot 11 ####
tmp.srt <- cm.srt
Idents(tmp.srt) <- 'group2'
deg <- FindMarkers(tmp.srt, ident.1 = 'C3KO MI', ident.2 = 'WT MI')
deg <- deg[deg$p_val_adj < 0.01,]
deg <- deg[order(deg$avg_log2FC, decreasing = T),]
deg$gene <- rownames(deg)
deg$direction <- ifelse(deg$avg_log2FC > 0, yes = 'Up in C3KO MI', no = 'Up in WT MI')

deg <- final_deg_list$CM
deg <- deg[order(deg$avg_log2FC, decreasing = T),]
p11 <- MarkerVolcano(deg, label_genes = c('Nppa', 'Nppb', 'Stat3', 'Il6st', 'Ppara', 'Pdlim3', 'Ryr2', 'Cox7b', 'Meg3')) +
        #scale_x_continuous(limits = c(-4.5, 2)) +
        theme_classic() +
        theme_Publication() +
        theme(aspect.ratio = 1)
p11
PlotPDF('1.11.volcano_cm_c3ko_mi_vs_wt_mi', 6, 6)
p11
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(deg, 'analysis/PART91.deg_cm_c3mi_vs_wtmi.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### Plot 12 ####
p12 <- CountCellBarPlot(cm.srt[, cm.srt$Cell_state != 'Atr. CM'],
                        group.var = 'group2', stack.var = 'Cell_state', percentage = T, stack.color = mycol_10)
p12
PlotPDF('1.12.bar.vent_cm_composition', 5, 5)
p12
dev.off()


#### Plot 13 ####
tmp.srt <- AddModuleScore2(cm.srt, features = list(Neo_MI_CM_marker$CM5,
                                                   Yap5sa_CM_Mac_CF_marker$CM2,
                                                   Yap_target$Yap_target_yap5sa_union,
                                                   CM_stress_genes),
                           names = c('Olson_CM5', 'aCM2', 'Yap_target', 'CM_stress'))
p13 <- BoxPlot(tmp.srt, feature = 'Olson_CM5', group.by = 'group2', cols = Color_group2) +
        labs(y = 'Expression Z-score', title = 'Regenerative CM Signature') +
        scale_y_continuous(limits = c(-2.5, 2.5)) +
        theme(aspect.ratio = 1)
p13
PlotPDF('1.13.box.olson_cm5_score', 5, 5)
p13
dev.off()


#### Plot 14 ####
deg <- readRDS('analysis/PART91.deg_cm_c3mi_vs_wtmi.rds')
gl <- split(rownames(deg), deg$avg_log2FC > 0)
names(gl) <- c('Down in C3KO', 'Up in C3KO MI')
str(gl)
enrich <- ModuleEnrichment(module_list = gl, human_or_mouse = 'mouse')
go <- enrich$GO$`GO_Up in C3KO MI`
go <- go[go$p.adjust < 0.05, ]
WriteCSV(go, 'PART91.cm_c3mi_vs_wtmi_go_bp')
react <- enrich$Reactome$`Reactome_Up in C3KO MI`
react <- react[react$p.adjust < 0.05, ]
WriteCSV(react, 'PART91.cm_c3mi_vs_wtmi_reactome')
kegg <- enrich$KEGG$`KEGG_Up in C3KO MI`
kegg <- kegg[kegg$p.adjust < 0.05, ]
WriteCSV(kegg, 'PART91.cm_c3mi_vs_wtmi_kegg')

go_oi <- c('oxidative phosphorylation',
           'ATP biosynthetic process',
           'cardiac muscle contraction',
           'adherens junction organization',
           'muscle cell differentiation'
)
go_plot <- go[go$Description %in% go_oi, ]
go_plot$Description <- factor(go_plot$Description, levels = go_plot$Description[order(go_plot$p.adjust)])
p14.1 <- ggplot(go_plot, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.3) +
        geom_text(aes(label = str_to_sentence(Description), x = 0.1, y = Description), hjust = 'left') +
        labs(x = '-Log10 adjusted p value', y = 'Top Enriched Terms', title = 'GO Biological Processes') +
        scale_y_discrete(limit = rev) &
        theme_classic() &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

kegg_oi <- c('Oxidative phosphorylation - Mus musculus (house mouse)',
             'Diabetic cardiomyopathy - Mus musculus (house mouse)',
             'Cardiac muscle contraction - Mus musculus (house mouse)',
             'Adherens junction - Mus musculus (house mouse)',
             'Dilated cardiomyopathy - Mus musculus (house mouse)',
             'Hypertrophic cardiomyopathy - Mus musculus (house mouse)'
)
kegg_plot <- kegg[kegg$Description %in% kegg_oi, ]
kegg_plot$Description <- factor(kegg_plot$Description, levels = kegg_plot$Description[order(kegg_plot$p.adjust)])
p14.2 <- ggplot(kegg_plot, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.3) +
        geom_text(aes(label = str_to_sentence(Description), x = 0.1, y = Description), hjust = 'left') +
        labs(x = '-Log10 adjusted p value', y = 'Top Enriched Terms', title = 'KEGG Pathways') +
        scale_y_discrete(limit = rev) &
        theme_classic() &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p14.1 / p14.2
PlotPDF('1.14.barplot.c3ko_up_cm_enriched_bp_kegg', 6, 6)
p14.1 / p14.2
dev.off()

#### Plot 15 ####
tmp.srt <- NormalizeData(full.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'genotype', 'condition', 'replicate')]
tmp.srt <- FindVariableFeatures(tmp.srt, nfeatures = 5e3)

psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'group1',
                           features = VariableFeatures(tmp.srt))$RNA
meta <- tmp.srt@meta.data[!duplicated(tmp.srt$group1), ]
rownames(meta) <- meta$group1
meta <- meta[colnames(psb), ]
meta$age <- as.numeric(as.vector(meta$age))
dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group1)
rld <- rlog(dds, blind = T)

rv <- rowVars(assay(rld))
select <- order(rv, decreasing = T)[seq_len(min(5000, length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
scree_plot <- data.frame(percentVar[1:7])
scree_plot[, 2]<- c(1:7)

colnames(scree_plot) <- c("variance", 'component_number')
p15 <- ggplot(scree_plot) +
        geom_bar(aes(x = component_number, y = variance), stat = "identity", fill = 'grey75', color = 'black') +
        labs(x = 'Principal Component', t = 'Variance Explained') +
        scale_x_continuous(breaks = 1:7) +
        theme_classic()
p15
PlotPDF('1.15.bar.global_pc_scree_plot', 5, 5)
p15
dev.off()

#### Plot 16 ####
p16 <- FeaturePlot2(cm.srt,
             features = c('G2M.Score', 'Top2a', 'Mki67', 'S.Score', 'Pcna', 'Ccne2'), 
             reduction = 'sub_clean_umap', ncol = 3, min.cutoff = 0, max.cutoff = 3)
PlotPDF('1.16.feat.cm_prolif_markers', 12, 8)
p16
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Tables ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Table 1 ####
deg <- readRDS('analysis/PART91.deg_cf_c3mi_vs_wtmi.rds')
WriteCSV(deg, 'PART91.CF_DEG.C3KO_MI_vs_WT_MI')


#### Table 2 ####
deg <- readRDS('analysis/PART91.deg_cm_c3mi_vs_wtmi.rds')
WriteCSV(deg, 'PART91.CM_DEG.C3KO_MI_vs_WT_MI')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### For GEO Submission  #####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta <- full.srt@meta.data[, c(
        'nCount_RNA',
        'nFeature_RNA',
        'pct_mito_RNA',
        'orig.name',
        'study',
        'method',
        'platform',
        'data_process',
        'tissue',
        'genotype',
        'condition',
        'sex',
        'age',
        'replicate',
        'S.Score',
        'G2M.Score',
        'Doublet_SC',
        'Doublet_SC_score',
        'Non_ambiguous',
        'Cell_type',
        'Cell_state',
        'group1',
        'group2'
)]
meta$global_umap_x <- full.srt@reductions$clean_umap@cell.embeddings[, 1]
meta$global_umap_y <- full.srt@reductions$clean_umap@cell.embeddings[, 2]
meta$cell_type_umap_x <- full.srt@reductions$sub_clean_umap@cell.embeddings[, 1]
meta$cell_type_umap_y <- full.srt@reductions$sub_clean_umap@cell.embeddings[, 2]

umi_raw <- GetAssayData(full.srt, slot = 'count', assay = 'RNA')
umi_norm <- GetAssayData(full.srt, slot = 'data', assay = 'RNA')

WriteCSV(meta, 'PART91.GEO_metadata')
WriteCSV(as.data.frame(umi_raw), 'PART91.GEO_raw_umi')
WriteCSV(as.data.frame(umi_norm), 'PART91.GEO_norm_umi')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Test #####
cm.srt$tmp <- paste(cm.srt$group2, cm.srt$sex)
cm.srt$tmp[cm.srt$group1 == 'C3KO 5dpSham Rep2'] <- 'C3KO Sham M 2'
cm.srt$tmp <- factor(cm.srt$tmp, levels = c(
        'WT Sham M',
        'WT MI M',
        'C3KO Sham M',
        'C3KO Sham M 2',
        'C3KO MI M',
        'WT Sham F',
        'WT MI F',
        'C3KO MI F'
))
CountCellBarPlot(cm.srt, group.var = 'tmp', stack.var = 'Cell_state', stack.color = mycol_10) |
        CountCellBarPlot(cm.srt, group.var = 'group2', stack.var = 'Cell_state', stack.color = mycol_10)

VlnPlot2(cm.srt, features = 'G2M.Score', group.by = 'tmp') +
        scale_y_continuous(limits = c(1, 5))

cm.srt <- AddModuleScore2(cm.srt, features = list(g2m.genes), names = 'G2M.Score', return_z = T)
cm.srt$G2M_cell <- c('Non-G2M CM', 'G2M CM')[(as.numeric(cm.srt$G2M.Score>1) +1)]
cm.srt$G2M_cell <- factor(cm.srt$G2M_cell, levels = c('Non-G2M CM', 'G2M CM'))
CountCellBarPlot(cm.srt, group.var = 'tmp', stack.var = 'G2M_cell', stack.color = mycol_10) |
        CountCellBarPlot(cm.srt, group.var = 'group2', stack.var = 'G2M_cell', stack.color = mycol_10)

tmp.srt <- cm.srt[, cm.srt$condition == '5dpMI']
tmp.srt$group2 <- droplevels(tmp.srt$group2)
tmp.srt$Cell_state <- droplevels(tmp.srt$Cell_state)
ChiSqP(srt.obj = tmp.srt, var1 = 'group2', var2 = 'Cell_state')
tmp.srt <- cm.srt[, cm.srt$condition != '5dpMI']
tmp.srt$group2 <- droplevels(tmp.srt$group2)
tmp.srt$Cell_state <- droplevels(tmp.srt$Cell_state)
ChiSqP(srt.obj = tmp.srt, var1 = 'group2', var2 = 'Cell_state')
tmp.srt <- cm.srt[, cm.srt$genotype != 'C3KO']
tmp.srt$group2 <- droplevels(tmp.srt$group2)
tmp.srt$Cell_state <- droplevels(tmp.srt$Cell_state)
ChiSqP(srt.obj = tmp.srt, var1 = 'group2', var2 = 'Cell_state')


tmp.srt <- cm.srt[, cm.srt$Cell_state != 'Atr. CM']
VlnPlot2(tmp.srt, features = 'Nppa', group.by = 'group2', pt.size = 1)

Idents(tmp.srt) <- 'group2'
VlnPlot(tmp.srt, feature = 'CC_Score', group.by = 'group2')
tmp.srt$Prolif <- tmp.srt$CC_Score>1
CountCellBarPlot(tmp.srt, group.var = 'group2', stack.var = 'Prolif')
PlotPDF('03.9.scatter.mp2_cell_cycle_score', 6, 6)
p9
dev.off()


CountCellBarPlot(cm.srt[, cm.srt$Cell_state != 'Atr. CM'], group.var = 'group2', stack.var = 'Cell_state', percentage = T)
VlnPlot(cm.srt, feature = 'Top2a', group.by = 'group2')

mk <- FindAllMarkers(cm.srt, only.pos = T)
mk <- mk[mk$p_val_adj < 0.05, ]
p9 <- MarkerHeatmap(cm.srt, marker.df = mk, raster = F, disp.min = 0)
labs(x = 'UMAP1', y = 'UMAP2', title = 'Cardiomyocyte')
p9
PlotPDF('1.09.heat_cm_cell_state_marker', 5, 5)
p9
dev.off()

tmp.srt <- AddModuleScore2(cm.srt, features = Neo_MI_CM_marker, names = names(Neo_MI_CM_marker))
FeaturePlot2(tmp.srt, features = names(Neo_MI_CM_marker), reduction = 'sub_clean_umap', pt.size = 3,
             raster = T, split.by = 'group2')
VlnPlot(tmp.srt, feature = 'CM5', group.by = 'group2')

tmp.srt <- AddModuleScore2(cf.srt, features = split(deg$gene, deg$direction), names = c('Dn', 'Up'))
FeaturePlot2(tmp.srt, features = c('Dn', 'Up'), reduction = 'sub_clean_umap', pt.size = 3, raster = T)

VlnPlot(tmp.srt, feature = 'CM1', group.by = 'group2')
DotPlot2(tmp.srt, feature = 'CM1', group.by = 'group2')



tmp.srt <- AddModuleScore2(cm.srt, features = Yap5sa_CM_Mac_CF_marker[1:2], names = c('aCM1', 'aCM2'))
FeaturePlot2(tmp.srt, features = c('aCM1', 'aCM2'), reduction = 'sub_clean_umap',
             pt.size = 3, raster = T, split.by = 'group2')
VlnPlot(tmp.srt, feature = 'aCM2', group.by = 'group2')
BoxPlot(tmp.srt, feature = 'aCM2', group.by = 'group2')


cm_deg <- FindAllMarkers(drop.srt, logfc.threshold = 0.5, only.pos = T)
cm_deg <- cm_deg[cm_deg$p_val_adj < 0.01, ]
cm_deg <- cm_deg[cm_deg$cluster %in% c('CM1', 'CM2'),]
gl <- split(cm_deg$gene, cm_deg$cluster)[2:1]
common <- intersect(gl$CM1, gl$CM2)
gl$CM1 <- gl$CM1[! gl$CM1 %in% common]
gl$CM2 <- gl$CM2[! gl$CM2 %in% common]
str(gl)

tmp.srt <- AddModuleScore2(cm.srt, features = gl, names = c('aCM1', 'aCM2'))
VlnPlot(tmp.srt, feature = 'aCM2', group.by = 'group2')
BoxPlot(tmp.srt, feature = 'aCM2', group.by = 'group2')
FeaturePlot2(tmp.srt, features = c('aCM1', 'aCM2'), reduction = 'sub_clean_umap',
             pt.size = 3, raster = T, split.by = 'group2')

