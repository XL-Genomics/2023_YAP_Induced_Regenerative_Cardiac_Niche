####--------------------------------------------------------------------------------------------------------------------
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate directories  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '0'
Step <- '04_Cd45sc_M2_DEG'
Project <- '2022_yap5sa_st_rli'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'mouse', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  DEG   ####
####--------------------------------------------------------------------------------------------------------------------
sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
sc.srt <- AddModuleScore2(sc.srt,
                          features = Yap5sa_CM_Mac_CF_marker,
                          names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)), return_z = T)
sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
Idents(sc.srt) <- 'sample_pub'
mac.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac', 'Mono')]
mac.srt$Cell_state <- droplevels(mac.srt$Cell_state)

## Find DEGs
mac.srt$C3ar1_high <- factor('Low', levels = c('High', 'Low'))
mac.srt$C3ar1_high[mac.srt$Zscore_C3ar1_Mac > 1] <- 'High'
DimPlot2(mac.srt, group.by = 'C3ar1_high', reduction = 'sub_hmn_umap')
Idents(mac.srt) <- 'C3ar1_high'
deg <- FindMarkers(mac.srt, ident.1 = 'High', ident.2 = 'Low',
                   logfc.threshold = 0.1,
                   test.use = 'MAST')
deg.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.25], human_or_mouse = 'mouse'),
                 'Low' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC < -0.25], human_or_mouse = 'mouse'))
degtop.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.5 & -log10(deg$p_val_adj) >= 30], human_or_mouse = 'mouse'),
                    'Low' = rownames(deg)[deg$avg_log2FC < -0.5 & -log10(deg$p_val_adj) >= 30])
p <- MarkerVolcano(deg, label = T, min.log2FC = 0.25, label_genes = unlist(degtop.list))
PlotPDF('01.1.vol.c3ar1_mac_deg', 10, 10)
p
dev.off()

enr <- ModuleEnrichment(deg.list, human_or_mouse = 'mouse')
tab1 <- enr$GO$GO_High
tab1 <- tab1[tab1$p.adjust <= 0.01, ]
tab2 <- enr$Reactome$Reactome_High
tab2 <- tab2[tab2$p.adjust <= 0.1, ]

terms <- c('myoblast differentiation',
           'myoblast fusion',
           'skeletal muscle tissue regeneration',
           'muscle cell proliferation',
           'muscle cell migration',
           'cell junction disassembly',
           'actin filament organization')
df <- tab1[tab1$Description %in% c(terms), ]
df$Description <- factor(df$Description, levels = df$Description[order(df$p.adjust)])
p <- ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.2) +
        #geom_text(aes(label = gene, x = 5, y = Description), color = 'black', size = 3) +
        scale_y_discrete(limit = rev) +
        theme_classic()
PlotPDF('01.2.bar.c3ar1_mac_upgene_go', 5, 5)
print(p)
dev.off()

WriteCSV(tab1, 'part04.1.c3ar1_mac_upgene_go_enrich')
WriteCSV(tab2, 'part04.2.c3ar1_mac_upgene_pathway_enrich')
deg$gene <- rownames(deg)
deg <- deg[deg$p_val_adj < 0.05, ]
WriteCSV(deg[order(deg$avg_log2FC, decreasing = T),], 'part04.3.c3ar1_mac_deg')

####--------------------------------------------------------------------------------------------------------------------
saveRDS(mac.srt, 'analysis/part04.merged_orig.srt.rds')
saveRDS(deg, 'analysis/part04.c3ar1_mac_deg.df.rds')
saveRDS(enr, 'analysis/part04.c3ar1_mac_upgene_enrich.list.rds')
####--------------------------------------------------------------------------------------------------------------------
