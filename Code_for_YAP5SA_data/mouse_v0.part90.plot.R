####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Yap5sa Spatial Transcriptomics -- Collaboration with Rich G. Li
####  2022-03-11 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate directories  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '0'
Step <- '90_Plot'
Project <- '2022_yap5sa_st_rli'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'mouse', Project, 'ithil')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  01. CM2 colocalization with Mac and CF  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st.srt <- readRDS('individual/part01.YAP5SA_st_ventricle.srt.rds')
DefaultAssay(st.srt) <- 'SCT'

st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = F,
                          names = paste0('Score_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))
st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = T,
                          names = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))

st.srt <- AddModuleScore2(st.srt, features = Yap_target[1:2], assay = 'SCT', return_z = F,
                          names = c('Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'))
st.srt <- AddModuleScore2(st.srt, features = Yap_target[1:2], assay = 'SCT', return_z = T,
                          names = c('Zscore_Y5SA_CM_Target1', 'Zscore_Y5SA_CM_Target2'))
RidgePlot(st.srt, group.by = 'Sample',
          features = c('Score_Y5SA_CM1', 'Score_Y5SA_CM2', 'Score_Y5SA_C3ar1_Mac', 'Score_Y5SA_C3_FB',
                       'Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'),
          ncol = 1) +
        NoLegend() +
        theme(axis.title = element_blank(), axis.text.y = element_blank())
shapiro.test(sample(st.srt$Score_Y5SA_CM1, size = 5000))$p.value
shapiro.test(sample(st.srt$Score_Y5SA_CM2, size = 5000))$p.value
shapiro.test(sample(st.srt$Score_Y5SA_C3ar1_Mac, size = 5000))$p.value
shapiro.test(sample(st.srt$Score_Y5SA_C3_FB, size = 5000))$p.value
shapiro.test(sample(st.srt$Score_Y5SA_C3_FB, size = 5000))$p.value


## Calculate z score of C3 and C3ar1
DefaultAssay(st.srt) <- 'alra'
st.srt <- ScaleData(st.srt)
RidgePlot(st.srt, group.by = 'Sample', slot = 'scale.data',
          features = c('C3', 'C3ar1'),
          ncol = 1) +
        NoLegend() +
        theme(axis.title = element_blank(), axis.text.y = element_blank())
shapiro.test(sample(st.srt@assays$alra@scale.data['C3',], size = 5000))$p.value
shapiro.test(sample(st.srt@assays$alra@scale.data['C3ar1', ], size = 5000))$p.value
st.srt$Score_C3 <- st.srt@assays$alra@scale.data['C3',]
st.srt$Score_C3ar1 <- st.srt@assays$alra@scale.data['C3ar1',]
st.srt$Zscore_C3 <- scale(st.srt$Score_C3)
st.srt$Zscore_C3ar1 <- scale(st.srt$Score_C3ar1)


## Plot expression pattern
minvals <- c(rep(-1, 3), 0, 0)
maxvals <- c(rep(1.5, 3), 2, 2)
p <- FeaturePlotST(srt = st.srt,
                   features = c('Zscore_Y5SA_CM1', 'Zscore_Y5SA_CM2',
                                'Zscore_Y5SA_CM_Target1',
                                'C3', 'C3ar1'),
                   title = c('CM1 score', 'CM2 score', 'Yap target score', 'C3', 'C3ar1'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = st.srt@misc$spot_scale*0.09,
                   ncol = 2)
p[[6]] <- p[[6]] + RestoreLegend() + labs(color='Z score')
p[[10]] <- p[[10]] + RestoreLegend() + labs(color='Expression')
p
PlotPDF('01.1.feat_st.yap5sa_marker_zscore', 6, 10)
print(p)
dev.off()


## Get cm-cf-mac colocalization probability
st.srt$Colocal_CM1_C3ar1_C3 <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3', 'Zscore_C3ar1'),
                         return_logp = T)
st.srt$Colocal_CM2_C3ar1_C3 <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3', 'Zscore_C3ar1'),
                       return_logp = T)
p_cutoff <- 0.1
minvals <- rep(0, 4)
maxvals <- rep(-log10(p_cutoff), 4)
p <- FeaturePlotST(srt = st.srt,
                   features = c('Colocal_CM1_C3ar1_C3',
                                'Colocal_CM2_C3ar1_C3'),
                   title = c('CM1 + C3ar1 Hi + C3 Hi',
                             'CM2 + C3ar1 Hi + C3 Hi'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = st.srt@misc$spot_scale*0.25,
                   ncol = 2) &
        scale_color_distiller(palette = 'Reds', limits = c(0, -log10(p_cutoff)),
                              direction = 0, values = c(0, 1))
p[[2]] <- p[[2]] + RestoreLegend() + labs(color='-Log10(p)')
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('01.2.feat_st.yap5sa_marker_colocolization_on_y5sa_st', 6, 5)
print(p)
dev.off()

st.srt$Psig_CM1 <- factor(paste0('p<', p_cutoff), levels = c('Not sig.', paste0('p<', p_cutoff)))
st.srt$Psig_CM2 <- factor(paste0('p<', p_cutoff), levels = c('Not sig.', paste0('p<', p_cutoff)))
st.srt$Psig_CM1[st.srt$Colocal_CM1_C3ar1_C3 < -log10(p_cutoff)] <- 'Not sig.'
st.srt$Psig_CM2[st.srt$Colocal_CM2_C3ar1_C3 < -log10(p_cutoff)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st.srt$Psig_CM1, st.srt$Sample)),
            as.data.frame(Table(st.srt$Psig_CM2, st.srt$Sample)))
df$CM <- rep(c('CM1', 'CM2'), each = 4)
df$Fraction = df$Freq/rep(rep(Table(st.srt$Sample), each = 2), 2)
p <- ggplot(df[df$Var1 == paste0('p<', p_cutoff),]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~CM) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('01.3.bar.yap5sa_marker_colocolization_on_y5sa_st', 4, 5)
print(p)
dev.off()

st.srt$Colocal_CM1_C3ar1_C3_prob <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3', 'Zscore_C3ar1'),
                       return_logp = F, return_prob_mat = F)
st.srt$Colocal_CM2_C3ar1_C3_prob <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3', 'Zscore_C3ar1'),
                       return_logp = F, return_prob_mat = F)
p <- VlnPlot(st.srt, split.by = 'Sample', same.y.lims = T, y.max = 1, adjust = 2, pt.size = -1,
          features = c('Colocal_CM1_C3ar1_C3_prob', 'Colocal_CM2_C3ar1_C3_prob'), cols = mycol_10,
          ncol = 2)
p[[1]] <- p[[1]] +
        labs(x = '', y = 'Colocolization probability', title = 'CM1 + C3ar1 Hi + C3 Hi')
p[[2]] <- p[[2]] +
        labs(x = '', y = '', title = 'CM2 + C3ar1 Hi + C3 Hi') +
        RestoreLegend()
p <- p & theme(aspect.ratio = 2, axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
p
PlotPDF('01.4.vln.yap5sa_marker_colocolization_on_y5sa_st', 8, 5)
print(p)
dev.off()


## Get cm-cf colocalization probability
st.srt$Colocal_CM1_C3 <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3'),
                       return_logp = T)
st.srt$Colocal_CM2_C3 <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3'),
                       return_logp = T)
p_cutoff <- 0.1
minvals <- rep(0, 4)
maxvals <- rep(-log10(p_cutoff), 4)
p <- FeaturePlotST(srt = st.srt,
                   features = c('Colocal_CM1_C3',
                                'Colocal_CM2_C3'),
                   title = c('CM1 + C3 Hi',
                             'CM2 + C3 Hi'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = st.srt@misc$spot_scale*0.25,
                   ncol = 2) &
        scale_color_distiller(palette = 'Reds', limits = c(0, -log10(p_cutoff)),
                              direction = 0, values = c(0, 1))
p[[2]] <- p[[2]] + RestoreLegend() + labs(color='-Log10(p)')
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('01.5.feat_st.yap5sa_cm_marker_and_c3_colocolization', 6, 5)
print(p)
dev.off()

st.srt$Colocal_CM1_C3_prob <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3'),
                       return_logp = F, return_prob_mat = F)
st.srt$Colocal_CM2_C3_prob <-
        GetColocalProb(st.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3'),
                       return_logp = F, return_prob_mat = F)
p <- VlnPlot(st.srt, split.by = 'Sample', same.y.lims = T, y.max = 1, adjust = 2, pt.size = -1,
             features = c('Colocal_CM1_C3_prob', 'Colocal_CM2_C3_prob'), cols = mycol_10,
             ncol = 2)
p[[1]] <- p[[1]] +
        labs(x = '', y = 'Colocolization probability', title = 'CM1 + C3 Hi')
p[[2]] <- p[[2]] +
        labs(x = '', y = '', title = 'CM2 + C3 Hi') +
        RestoreLegend()
p <- p & theme(aspect.ratio = 2, axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
p
PlotPDF('01.6.vln.yap5sa_cm_marker_and_c3_colocolization', 8, 5)
print(p)
dev.off()


## Get cm-cf colocalization probability on control only
st_ctl.srt <- st.srt[, st.srt$Sample == 'Control']
st_ctl.srt@images$YAP5SA <- NULL
st_ctl.srt$Sample <- droplevels(st_ctl.srt$Sample)

st_ctl.srt$Colocal_CM1_C3 <-
        GetColocalProb(st_ctl.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3'),
                       return_logp = T)
st_ctl.srt$Colocal_CM2_C3 <-
        GetColocalProb(st_ctl.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3'),
                       return_logp = T)
p_cutoff <- 0.001
minvals <- rep(0, 4)
maxvals <- rep(-log10(p_cutoff), 4)
p <- FeaturePlotST(srt = st_ctl.srt,
                   features = c('Colocal_CM1_C3',
                                'Colocal_CM2_C3'),
                   title = c('CM1 + C3 Hi',
                             'CM2 + C3 Hi'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = st_ctl.srt@misc$spot_scale*0.25,
                   ncol = 1) &
        scale_color_distiller(palette = 'Reds', limits = c(0, -log10(p_cutoff)),
                              direction = 0, values = c(0, 1))
p[[2]] <- p[[2]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('01.7.feat_st.yap5sa_cm_marker_and_c3_colocolization_control_only', 4, 5)
print(p)
dev.off()

st_ctl.srt$Colocal_CM1_C3_prob <-
        GetColocalProb(st_ctl.srt, meta_features = c('Zscore_Y5SA_CM1', 'Zscore_C3'),
                       return_logp = F, return_prob_mat = F)
st_ctl.srt$Colocal_CM2_C3_prob <-
        GetColocalProb(st_ctl.srt, meta_features = c('Zscore_Y5SA_CM2', 'Zscore_C3'),
                       return_logp = F, return_prob_mat = F)
p <- VlnPlot(st_ctl.srt, split.by = 'Sample', same.y.lims = T, y.max = 1, adjust = 1, pt.size = -1,
             features = c('Colocal_CM1_C3_prob', 'Colocal_CM2_C3_prob'), cols = mycol_10,
             ncol = 2)
p[[1]] <- p[[1]] +
        labs(x = '', y = 'Colocolization probability', title = 'CM1 + C3 Hi')
p[[2]] <- p[[2]] +
        labs(x = '', y = '', title = 'CM2 + C3 Hi') +
        RestoreLegend()
p <- p & theme(aspect.ratio = 2, axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
p
PlotPDF('01.8.vln.yap5sa_cm_marker_and_c3_colocolization_control_only', 8, 5)
print(p)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  02. Olson CM2 score and colocalization  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st2.srt <- readRDS('individual/part01.2021_NatComm_EOlson.srt.rds')
st2.srt <- AddModuleScore2(st2.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = F,
                           names = paste0('Score_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))
st2.srt <- AddModuleScore2(st2.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = T,
                           names = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))
st2.srt <- AddModuleScore2(st2.srt, features = Yap_target[1:2], assay = 'SCT', return_z = F,
                           names = c('Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'))
st2.srt <- AddModuleScore2(st2.srt, features = Yap_target[1:2], assay = 'SCT', return_z = T,
                           names = c('Zscore_Y5SA_CM_Target1', 'Zscore_Y5SA_CM_Target2'))

st2.srt$tmp <- 'ALL'
RidgePlot(st2.srt, features = c("Score_Y5SA_CM1", "Score_Y5SA_CM2", "Score_Y5SA_C3ar1_Mac", "Score_Y5SA_C3_FB",
                                'Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'), group.by = 'tmp',
          ncol = 1)
shapiro.test(st2.srt$Score_Y5SA_CM1)$p.value
shapiro.test(st2.srt$Score_Y5SA_CM2)$p.value
shapiro.test(st2.srt$Score_Y5SA_C3ar1_Mac)$p.value
shapiro.test(st2.srt$Score_Y5SA_C3_FB)$p.value
shapiro.test(st2.srt$Score_Y5SA_C3_FB)$p.value

minvals <- rep(-1, 5)
maxvals <- rep(2, 5)
scale <- c(2.48, 4, 1.3, 2.25)
p <- FeaturePlotST(srt = st2.srt,
                   features = c(paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)),
                                'Zscore_Y5SA_CM_Target1'),
                   title = c('CM1 score', 'CM2 score', 'C3ar1+ Mac score', 'C3+ FB score', 'Yap target score'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.4,
                   ncol = 4)
p[[20]] <- p[[20]] + RestoreLegend() + labs(color='Z score')
p
PlotPDF('02.1.feat_st.yap5sa_marker_zscore', 12, 10)
print(p)
dev.off()

st2.srt$Colocal_CM1_C3ar1Mac_C3FB <-
        GetColocalProb(st2.srt, meta_features = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)[c(1, 3, 4)]),
                       return_logp = T)
st2.srt$Colocal_CM2_C3ar1Mac_C3FB <-
        GetColocalProb(st2.srt, meta_features = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)[c(2, 3, 4)]),
                       return_logp = T)

minvals <- rep(0, 4)
maxvals <- rep(3, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_C3ar1Mac_C3FB',
                                'Colocal_CM2_C3ar1Mac_C3FB'),
                   title = c('CM1, C3ar1+ Mac, C3+ FB',
                             'CM2, C3ar1+ Mac, C3+ FB'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.6,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 3), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')

PlotPDF('02.2.feat_st.yap5sa_marker_colocolization_on_olson_st', 12, 5)
print(p)
dev.off()


st2.srt$Psig_CM1 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM2 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM1[st2.srt$Colocal_CM1_C3ar1Mac_C3FB < -log10(0.1)] <- 'Not sig.'
st2.srt$Psig_CM2[st2.srt$Colocal_CM2_C3ar1Mac_C3FB < -log10(0.1)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2, st2.srt$Sample)))
df$CM <- rep(c('CM1', 'CM2'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.1',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~CM) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('02.3.bar.yap5sa_marker_colocolization_on_olson_st', 4, 5)
print(p)
dev.off()


st2.srt$Colocal_CM1_YapTarget <-
        GetColocalProb(st2.srt, meta_features = c('Score_Y5SA_CM1', 'Score_Y5SA_CM_Target1'),
                       return_logp = T)
st2.srt$Colocal_CM2_YapTarget <-
        GetColocalProb(st2.srt, meta_features = c('Score_Y5SA_CM2', 'Score_Y5SA_CM_Target1'),
                       return_logp = T)

minvals <- rep(0, 4)
maxvals <- rep(3, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_YapTarget',
                                'Colocal_CM2_YapTarget'),
                   title = c('CM1, YapHigh',
                             'CM2, YapHigh'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.6,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 3), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')

PlotPDF('02.4.feat_st.yap5sa_marker_yap_target_colocolization_on_olson_st', 12, 5)
print(p)
dev.off()

st2.srt$Psig_CM1 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM2 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM1[st2.srt$Colocal_CM1_YapTarget < -log10(0.1)] <- 'Not sig.'
st2.srt$Psig_CM2[st2.srt$Colocal_CM2_YapTarget < -log10(0.1)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2, st2.srt$Sample)))
df$CM <- rep(c('CM1', 'CM2'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.1',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~CM) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('02.5.bar.yap5sa_marker__yap_target_colocolization_on_olson_st', 4, 5)
print(p)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  03. Human C3+ Fibro  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fb2 <- readRDS('analysis/part02.human_merged_orig.srt.rds')
fb2$Dataset <- factor('Litvinukova et al', levels = c('Tucker et al', 'Litvinukova et al'))
fb2$Dataset[is.na(fb2$source)] <- 'Tucker et al'

p1 <- FeaturePlot3(fb2, features = c('C3'), reduction = 'hmn_umap', adjust = 1.5) +
        labs(x = '', y = '', title = 'C3 expression')
p2 <- FeaturePlot3(fb2, features = c('Score_C3_FB'), reduction = 'hmn_umap', adjust = 1.5) +
        labs(x = '', y = '', title = 'C3+ Fibroblast score')
p3 <- DimPlot2(fb2, reduction = 'hmn_umap', group.by = 'Dataset', cols = mycol_10) +
        labs(x = '', y = '', title = '', color = 'Dataset')

p <- wrap_plots(p1, p2, p3, ncol = 3)
PlotPDF('03.1.umap.human_c3_cf', 15, 5)
p
dev.off()

deg <- readRDS('analysis/part02.human_c3_cf_deg.df.rds')
degtop.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.5], human_or_mouse = 'human'),
                    'Low' = rownames(deg)[deg$avg_log2FC < -0.25])
goi <- list(Up = c('Cxcl14', 'Mt2', 'Fstl1', 'Tmsb4x', 'Cst3', 'Fth1', 'Ly6a', 'Mt1', 'Apoe', 'Serpina3n',
                   'Igfbp4', 'Gsn', 'Gpx3', 'Serping1', 'Dcn'),
            Down = c('Lhfp', 'Lama2', 'Gpc6', 'Mast4'))
goi$Up <- str_to_upper(goi$Up)
goi$Down <- str_to_upper(goi$Down)

df <- deg
df$gene <- rownames(df)
df$log_p <- -log10(df$p_val_adj)
df_label <- df[df$gene %in% unlist(goi), ]
p <- ggplot(df) +
        geom_point(aes(y = avg_log2FC, x = log_p), size = 0.8) +
        labs(x = '-LogP', y = 'Log2FC') +
        theme_classic() +
        theme(aspect.ratio = 1, panel.grid.major = element_line(colour = 'grey90'),
              axis.line = element_blank(), axis.ticks = element_blank()) +
        ggrepel::geom_text_repel(data = df_label, aes(label = gene, y = avg_log2FC, x = log_p),
                                 color = 'black', size = 3, fontface = "italic", nudge_x = 1,  nudge_y = 0.2,
                                 max.overlaps = 100, min.segment.length = 0, force = 3)

p
PlotPDF('03.2.vol.c3_cf_deg', 5, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  04. Olson snRNA-seq CM scoring  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cm.srt <- readRDS('individual/part01.2020_DevCell_EOlson.srt.rds')
cm.srt$Sample_pub <- revalue(cm.srt$Sample, replace = c(
        'P1_1Sham' = 'P1 Sham - 1 day',
        'P1_1MI' = 'P1 MI - 1 day',
        'P1_3Sham' = 'P1 Sham - 3 day',
        'P1_3MI' = 'P1 MI - 3 day',
        'P8_1Sham' = 'P8 Sham - 1 day',
        'P8_1MI' = 'P8 MI - 1 day',
        'P8_3Sham' = 'P8 Sham - 3 day',
        'P8_3MI'  = 'P8 MI - 3 day'
        ))
cm.srt$Sample_pub <- factor(cm.srt$Sample_pub, levels = c(
        'P1 Sham - 1 day',
        'P1 MI - 1 day',
        'P1 Sham - 3 day',
        'P1 MI - 3 day',
        'P8 Sham - 1 day',
        'P8 MI - 1 day',
        'P8 Sham - 3 day',
        'P8 MI - 3 day'
))
cm2.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')

p1 <- FeaturePlot2(cm.srt, features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'), raster = T, pt.size = 2)
p2 <- FeaturePlot2(cm2.srt, features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'), raster = T, pt.size = 3)
PlotPDF('04.1.feat.cm2_marker_screen', 15, 30)
p1/p2
dev.off()

#cm.srt$Zscore_y5sa_CM2_pos <- cm.srt$Zscore_y5sa_CM2-min(cm.srt$Zscore_y5sa_CM2)
cm.srt$CM2_like <- factor('Other CMs', levels = c('Yap5sa CM2-like', 'Other CMs'))
cm.srt$CM2_like[cm.srt$Zscore_y5sa_CM2 >= 1] <- 'Yap5sa CM2-like'
cm.srt$CM2_like2 <- factor(cm.srt$CM2_like, levels = c('Other CMs', 'Yap5sa CM2-like'))
DefaultAssay(cm.srt) <- 'RNA'

p1 <- DimPlot2(cm.srt, cols = mycol_10, raster = T, pt.size = 1, label = T) +
        NoLegend() +
        labs(x = '', y = '', title = 'CM subtypes -- Cui et al') +
        theme(title = element_text(hjust = 1))
p2 <- FeaturePlot2(cm.srt, features = 'Zscore_y5sa_CM2', pt.size = 1, raster = T,
                   min.cutoff = 'q5', max.cutoff = 'q95', order = F) +
        labs(x = '', y = '', col = 'Z score', title = 'Yap5sa CM2 score')
p3 <- DimPlot2(cm.srt, cols = c(mycol_10[2], 'grey85'), raster = T, pt.size = 1, label = F, group.by = 'CM2_like') +
        labs(x = '', y = '', title = 'Yap5sa CM2-like', caption = '*Yap5sa CM2 score >= 1', color = '') +
        theme(title = element_text(hjust = 1))

p4 <- VlnPlot2(cm.srt, features = c('Ctgf', 'Rtn4' ,'Sorbs2', 'Acta2', 'Eef2', 'Ptrf'),
               ncol = 3, cols = mycol_10) &
        labs(y = 'Expression')

p5 <- CountCellBarPlot(cm.srt, group.var = 'Sample_pub', stack.var = 'CM2_like2',
                       percentage = T, stack.color = c('grey85', mycol_10[2])) +
        labs(x = '', y = 'Fraction of nuclei', fill = '')

p <- wrap_plots(p1 + p2 + p3, p4 | p5, ncol = 1)
PlotPDF('04.2.combine.cm2_similarity_olson_sn', 15, 14)
p
dev.off()

p <- VlnPlot2(cm.srt,
              # features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'),
              features = c('Ctgf', 'Flnc' ,'Serpinh1', 'Acta2', 'Eef2', 'Rtn4'),
              group.by = 'CM2_like2',
              ncol = 3, cols = mycol_10) &
        labs(y = 'Expression') &
        theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank())
p[[6]] <- p[[6]] + RestoreLegend()
p
PlotPDF('04.3.vln.cm2_like_marker_example', 8, 4)
p
dev.off()

cm.srt$tmp <- revalue(cm.srt$CM_subtype, c(CM4 = 'CM4/5',
                                           CM5 = 'CM4/5'))
p <- VlnPlot2(cm.srt,
              # features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'),
              features = c('Ctgf', 'Flnc' ,'Ptrf', 'Acta2', 'Eef2', 'Ahnak'),
              group.by = 'tmp',
              ncol = 3, cols = mycol_10) &
        labs(y = 'Expression')

p[[1]] <- p[[1]] + theme(axis.text.x.bottom = element_blank())
p[[2]] <- p[[2]] + theme(axis.text.x.bottom = element_blank())
p[[3]] <- p[[3]] + theme(axis.text.x.bottom = element_blank())
p
PlotPDF('04.4.vln.cm2_like_marker_example_cm45_merge', 6, 4.5)
p
dev.off()

p1 <- DotPlot(cm.srt, features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'), col.min = 0, cols = c('blue', 'red')) +
        theme(aspect.ratio = 0.5)
p2 <- DotPlot(cm.srt, features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'), col.min = 0, cols = c('blue', 'red'),
              group.by = 'tmp') +
        theme(aspect.ratio = 0.4) +
        NoLegend()
p3 <- DotPlot(cm.srt, features = c(Yap5sa_CM_Mac_CF_marker$CM2, 'Ctgf', 'Ptrf'), col.min = 0, cols = c('blue', 'red'),
              group.by = 'CM2_like2') +
        theme(aspect.ratio = 0.2) +
        NoLegend()
p <- wrap_plots(p1, p2, p3, ncol = 1) &
        RotatedAxis() &
        labs(y = '')
p
PlotPDF('04.5.dot.cm2_like_marker_example_cm45_merge', 8, 10)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  05. YAP5SA ST C3ar1 Ligand Receptor   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get data
st.srt <- readRDS('individual/part01.YAP5SA_st_ventricle.srt.rds')

drop.srt <- readRDS('individual/part01.y5sa_dropseq_all.srt.rds')
cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
Idents(cm.srt) <- 'Experiment'

sc.srt <- readRDS('/Volumes/ithil/project/nsp13/rdata/mouse_v0/integrated/11.cd45_annotated.simple.srt.rds')
sc.srt <- sc.srt[, sc.srt$sample_pub %in% c('WT', 'YAP5SA')]
sc.srt <- AddModuleScore2(sc.srt,
                          features = Yap5sa_CM_Mac_CF_marker,
                          names = paste0('Zscore_', names(Yap5sa_CM_Mac_CF_marker)), return_z = T)
sc.srt$sample_pub <- droplevels(sc.srt$sample_pub)
Idents(sc.srt) <- 'sample_pub'
mac.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac')]

## Check L-R specificity
p <- DotPlot2(drop.srt, features = c('Igf1', 'Adam15', 'Tnfsf12'), split.by = 'Experiment', cols = 'RdYlBu') /
        DotPlot2(drop.srt, features = c('Igf1r', 'Itgb1', 'Tnfrsf12a'), split.by = 'Experiment', cols = 'RdYlBu')
PlotPDF('05.1.dot.lig_rec_specificity_in_dropseq', 6, 12)
p
dev.off()

## Check L-R expression in M2 and CM
mac.srt@active.ident <- revalue(mac.srt@active.ident, replace = c('WT' = 'Control'))
levels(cm.srt) <- c("Ctrl", "YAP5SA")
cm.srt@active.ident <- revalue(cm.srt@active.ident, replace = c('Ctrl' = 'Control'))
p1 <- VlnPlot2(mac.srt, features = c('Tnfsf12', 'Igf1', 'Adam15'), pt.size = -1, raster = F, assay = 'RNA')
p2 <- VlnPlot2(cm.srt, features = c('Tnfrsf12a', 'Igf1r', 'Itgb1'), pt.size = -1, raster = F, assay = 'alra')
PlotPDF('05.2.vln.lig_rec_specificity_in_cd45sc_cmdropseq', 6, 6)
p1/p2
dev.off()

## Ligand on Mac cluster and Receptor on CM
p1 <- DimPlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], group.by = 'Cell_state',
               raster = T, pt.size = 1.5,
               reduction = 'sub_hmn_umap', cols = mycol_10) +
        labs(title = 'Myeloid (Cd45+ scRNA-seq)', x = '', y = '')
p2 <- FeaturePlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], features = 'Zscore_C3ar1_Mac',
                   reduction = 'sub_hmn_umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T, pt.size = 1.5)  +
        labs(title = 'C3ar1+ Mac score', fill = 'Z score')
p3 <- DimPlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], group.by = 'sample_pub',
               reduction = 'sub_hmn_umap', raster = T, pt.size = 1.5) +
        labs(title = 'Genotype', x = '', y = '')
DefaultAssay(cm.srt) <- 'alra'
p4 <- DimPlot2(cm.srt, group.by = 'CM_State',
               raster = T, pt.size = 1.5,
               reduction = 'oriumap', cols = mycol_10) +
        labs(title = 'CM (Drop-seq)', x = '', y = '')
p5 <- FeaturePlot2(cm.srt, features = 'Zscore_y5sa_CM2',
                   reduction = 'oriumap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T, pt.size = 1.5)  +
        labs(title = 'CM2 score', fill = 'Z score')
p6 <- DimPlot2(cm.srt, group.by = 'Experiment',
               reduction = 'oriumap', raster = T, pt.size = 1.5) +
        labs(title = 'Genotype', x = '', y = '')

p7 <- FeaturePlot2(sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')], features = c('Tnfsf12', 'Igf1', 'Adam15'),
                   reduction = 'sub_hmn_umap', max.cutoff = 'q95', raster = T, pt.size = 1.5, ncol = 1)
p8 <- FeaturePlot2(cm.srt, features = c('Tnfrsf12a', 'Igf1r', 'Itgb1'),
                   reduction = 'oriumap', min.cutoff = 'q5', max.cutoff = 'q99', raster = T, pt.size = 1.5, ncol = 1)

PlotPDF('05.3.dim.cd45_m2', 10, 20)
wrap_plots(p1, p4, p2, p5, p3, p6, p7[[1]], p8[[1]], p7[[2]], p8[[2]], p7[[3]], p8[[3]], ncol = 2)
dev.off()


## Compute Mac-L colocalization
lr.df <- data.frame(from = c('Tnfsf12', 'Igf1', 'Adam15'),  to = c('Tnfrsf12a', 'Igf1r', 'Itgb1'))
st.srt <- AddModuleScore2(st.srt, features = Yap5sa_CM_Mac_CF_marker,
                          names = paste0('Score_', names(Yap5sa_CM_Mac_CF_marker)))
st.srt <- AddModuleScore2(st.srt, features = as.list(lr.df$from),
                          names = paste0('Ligand_score_', lr.df$from))
st.srt <- AddModuleScore2(st.srt, features = as.list(lr.df$to),
                          names = paste0('Receptor_score_', lr.df$to))
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Score_C3ar1_Mac', paste0('Ligand_score_', lig)))
        y <- GetColocalProb(st.srt, meta_features = c('Score_C3ar1_Mac', paste0('Ligand_score_', lig)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_M2_', lig), paste0('Colocal_M2_',  lig))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}
plist <- list()
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_M2_', lig)),
                                    title = c(paste0('M2 vs. ', lig)),
                                    minvals = c(0, 0, 0, 0), maxvals = c(3, 3, 3, 3),
                                    pt.sizes = st.srt@misc$spot_scale*0.35, ncol = 2) &
                scale_color_distiller(palette = 'Reds',  values = c(0, 0.5, 1), direction = 0, limits = c(0, 3))
}
plist[[3]][[2]] <- plist[[3]][[2]] + RestoreLegend() + labs(color='-Log10(p)')
p <- wrap_plots(plist, ncol = 1)
PlotPDF('05.4.feat_st.m2_lig_colocal', 7, 9)
print(p)
dev.off()


## Compute CM-R colocalization
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Score_CM2', paste0('Receptor_score_', rec)))
        y <- GetColocalProb(st.srt, meta_features = c('Score_CM2', paste0('Receptor_score_', rec)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM2_', rec), paste0('Colocal_CM2_',  rec))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}
plist <- list()
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_CM2_', rec)),
                                    title = c(paste0('CM2 vs. ', rec)),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1,1,1,1),
                                    pt.sizes = st.srt@misc$spot_scale*0.35, ncol = 2) &
                scale_color_distiller(palette = 'Reds',  values = c(0, 0.5, 1), direction = 0, limits = c(0, 1))
}
plist[[3]][[2]] <- plist[[3]][[2]] + RestoreLegend() + labs(color='-Log10(p)')
p <- wrap_plots(plist, ncol = 1)
PlotPDF('05.5.feat_st.cm2_rec_colocal', 7, 9)
print(p)
dev.off()

for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        x <- GetColocalProb(st.srt, meta_features = c('Score_CM1', paste0('Receptor_score_', rec)))
        y <- GetColocalProb(st.srt, meta_features = c('Score_CM1', paste0('Receptor_score_', rec)))
        df <- data.frame(x, y)
        colnames(df) <- c(paste0('Colocal_CM1_', rec), paste0('Colocal_CM1_',  rec))
        rownames(df) <- Cells(st.srt)
        st.srt <- AddMetaData(st.srt, metadata = df)
}
plist <- list()
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt,
                                    features = c(paste0('Colocal_CM1_', rec)),
                                    title = c(paste0('CM1 vs. ', rec)),
                                    minvals = c(0, 0, 0, 0), maxvals = c(1,1,1,1),
                                    pt.sizes = st.srt@misc$spot_scale*0.35, ncol = 2) &
                scale_color_distiller(palette = 'Reds',  values = c(0, 0.5, 1), direction = 0, limits = c(0, 1))
}
plist[[3]][[2]] <- plist[[3]][[2]] + RestoreLegend() + labs(color='-Log10(p)')
p <- wrap_plots(plist, ncol = 1)
PlotPDF('05.6.feat_st.cm1_rec_colocal', 7, 9)
print(p)
dev.off()

## Plot only yap5sa
st.srt2 <- st.srt[ , st.srt$Sample == 'YAP5SA']
st.srt2$Sample <- droplevels(st.srt2$Sample)
plist <- list()
for(i in 1:nrow(lr.df)){
        lig <- lr.df$from[i]
        rec <- lr.df$to[i]
        plist[[i]] <- FeaturePlotST(srt = st.srt2,
                                    features = c(paste0('Colocal_M2_', lig),
                                                 paste0('Colocal_CM1_', rec),
                                                 paste0('Colocal_CM2_', rec)),
                                    title = c(paste0('C3ar1 Mac + ', lig, ' Hi'),
                                              paste0('CM1 + ', rec, ' Hi'),
                                              paste0('CM2 + ', rec, ' Hi')),
                                    minvals = rep(0, 6), maxvals = rep(2, 6),
                                    pt.sizes = st.srt@misc$spot_scale*0.35, ncol = 1) &
                scale_color_distiller(palette = 'Reds',  values = c(0, 0.5, 1), direction = 0, limits = c(0, 2))
}
plist[[3]][[3]] <- plist[[3]][[3]] + RestoreLegend() + labs(color='-Log10(p)')
p <- wrap_plots(plist, ncol = 3)
PlotPDF('05.7.feat_st.m2_lig_cm_rec_colocal', 10, 9)
print(p)
dev.off()


mtx <- matrix(NA, 3, 3)
colnames(mtx) <- paste(lr.df$from, lr.df$to, sep = ":")
rownames(mtx) <- c('C3ar1 Mac + Ligand', 'CM1 + Receptor', 'CM2 + Receptor')
for(i in 1:3){
        mtx[1, i] <- Table(st.srt2[[paste0('Colocal_M2_', lr.df$from[i])]]>=1)[2]/ncol(st.srt2)
}
for(i in 1:3){
        mtx[2, i] <- Table(st.srt2[[paste0('Colocal_CM1_', lr.df$to[i])]]>=1)[2]/ncol(st.srt2)
        mtx[3, i] <- Table(st.srt2[[paste0('Colocal_CM2_', lr.df$to[i])]]>=1)[2]/ncol(st.srt2)
}
mtx[is.na(mtx)] <- 0
df <- reshape2::melt(mtx)

p <- ggplot(df) +
        geom_bar(aes(x = Var1, y = value, fill = Var2), stat = 'identity') +
        scale_fill_manual(values = mycol_10) +
        facet_wrap(~Var2) +
        theme_classic()+
        RotatedAxis() +
        NoLegend() +
        labs(y = 'Fraction of spots (p < 0.1)', x = '')
p
PlotPDF('05.8.bar.lig_rec_sig_fraction', 4.5, 5)
print(p)
dev.off()


sc.srt <- AddModuleScore2(sc.srt, features = list('Csf1r', 'C3ar1', Yap5sa_CM_Mac_CF_marker$C3ar1_Mac),
                           names = c('Score_Csf1r', 'Score_C3ar1', 'Score_C3ar1_Mac'),
                           return_z = T, nbin = 20)
sc.srt$C3ar1_high <- factor('C3ar1+', levels = c('C3ar1+', 'C3ar1-'))
sc.srt$C3ar1_high[sc.srt$Score_C3ar1_Mac <= 0.5] <- 'C3ar1-'
mye.srt <- sc.srt[, sc.srt$Cell_type %in% c('DC', 'Mac', 'Mono')]
p <- DimPlot2(mye.srt, group.by = 'C3ar1_high', reduction = 'sub_hmn_umap', cols = mycol_10[2:1]) +
        labs(x = '', y = '', title = 'C3ar1+ MPs')
p
PlotPDF('05.9.umap.C3ar1+MP', 6, 6)
print(p)
dev.off()

df <- CountCellBarPlot(sc.srt, group.var = 'sample_pub', stack.var = 'Cell_state', stack.color = mycol_30)$data
df <- df[grepl('Mac|Mono', df$Cell_state) ,]
df$pct <- df$count/as.vector(Table(sc.srt$sample_pub[sc.srt$Cell_type!='Fibro']))
df2 <- df[seq(1,nrow(df),2), ]
df2$fc <- NA
for(i in 1:8) {df2$fc[i] <- df$pct[i*2]/df$pct[i*2-1]}
p1 <- ggplot(df2)+
        geom_bar(aes(x = Cell_state, y = fc, fill = Cell_state ), stat = 'identity') +
        scale_fill_manual(values = mycol_10) +
        geom_hline(yintercept = 1, linetype = 'dashed', cols = 'red4') +
        labs(x = '', y = 'Fold-change of normalized cell number', title = 'Yap5sa vs. Control', fill = '') +
        theme_classic() +
        theme(axis.text.x.bottom = element_blank())
p1
PlotPDF('05.10.bar.all_mye_enrichment', 4.5, 5)
print(p1)
dev.off()

df <- CountCellBarPlot(sc.srt[,sc.srt$Cell_type!='Fibro'],
                       group.var = 'sample_pub', stack.var = 'C3ar1_high', stack.color = mycol_10)$data
df$pct <- df$count/as.vector(Table(sc.srt$sample_pub[sc.srt$Cell_type!='Fibro']))
df2 <- df[seq(1,nrow(df),2), ]
df2$fc <- NA
for(i in 1:2) {df2$fc[i] <- df$pct[i*2]/df$pct[i*2-1]}
p2 <- ggplot(df2)+
        geom_bar(aes(x = C3ar1_high, y = fc, fill = C3ar1_high ), stat = 'identity') +
        geom_hline(yintercept = 1, linetype = 'dashed', cols = 'red4') +
        scale_fill_manual(values = mycol_10[2:1]) +
        labs(x = '', y = 'Fold-change of normalized cell number', title = 'Yap5sa vs. Control', fill = '') +
        theme_classic() +
        theme(axis.text.x.bottom = element_blank())
p2
PlotPDF('05.11.bar.C3ar1+MP_enrichment', 4.5, 5)
print(p2)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  06. C3ar1 mac deg and go   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mac.srt <- readRDS('analysis/part04.merged_orig.srt.rds')

deg <- readRDS('analysis/part04.c3ar1_mac_deg.df.rds')
deg$gene <- rownames(deg)
WriteCSV(deg[abs(deg$avg_log2FC)>0.5 & -log10(deg$p_val_adj)>=30, ], 'part90.06.c3ar1_mac_all_deg')

degtop.list <- list('High' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC > 0.5 & -log10(deg$p_val_adj) >= 30],
                                            human_or_mouse = 'mouse'),
                    'Low' = RemoveRiboMito(rownames(deg)[deg$avg_log2FC < -0.5 & -log10(deg$p_val_adj) >= 30],
                                           human_or_mouse = 'mouse'))
goi <- c('C3ar1', 'C5ar1', 'Rbpj', 'Gdf15', 'Igf1', 'Adam15', 'Tnfsf12', 'F13a1',
         'Mertk', 'Mrc1', 'Il10', 'Maf', 'Cd163', 'Cd68', 'Cd36')
df <- deg
df$gene <- rownames(df)
df$log_p <- -log10(df$p_val_adj)
df_label <- df[df$gene %in% goi, ]
p <- ggplot(df[abs(deg$avg_log2FC) >= 0.25, ]) +
        geom_point(aes(y = avg_log2FC, x = log_p), size = 0.5) +
        labs(x = '-LogP', y = 'Log2FC') +
        theme_classic() +
        theme(panel.grid.major = element_line(colour = 'grey90'),
              axis.line = element_blank(), axis.ticks = element_blank()) +
        ggrepel::geom_text_repel(data = df_label, aes(label = gene, y = avg_log2FC, x = log_p),
                                 color = 'red', size = 3, fontface = "italic", nudge_x = 1,  nudge_y = 0.2,
                                 max.overlaps = 100, min.segment.length = 0, force = 3)

p
PlotPDF('06.1.vol.c3ar1_mac_deg', 5, 3)
p
dev.off()

enr <- readRDS('analysis/part04.c3ar1_mac_upgene_enrich.list.rds')
tab1 <- enr$GO$GO_High
tab1 <- tab1[tab1$p.adjust <= 0.01, ]

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
PlotPDF('06.2.bar.c3ar1_mac_upgene_go', 5, 5)
print(p)
dev.off()

lr_network <- readRDS('analysis/part03.lr_network.rds')
lr_network <- lr_network[lr_network$bonafide, c('from', 'to')]
ligand <- U(lr_network$from)
upgene <- rownames(deg)[deg$avg_log2FC > 0.25 & -log10(deg$p_val_adj) >= 2]
upligand <- intersect(upgene, ligand)
WriteCSV(deg[upligand,], 'part90.06.c3ar1_mac_all_upregulated_ligand')

p <- DotPlot2(mac.srt, upligand, cols = 'RdBu') +
        labs(title = 'Upregulated ligand genes (C3ar1 High Mac vs. C3ar1 Low Mac)')
PlotPDF('06.3.dot.c3ar1_mac_upgene_ligand', 10, 3)
print(p)
dev.off()

GenePCC <- function(
                matrix,
                method.use = 'pearson' ## c("pearson", "kendall", "spearman")
){
        num_genes <- dim(matrix)[1]
        pcc_mat <- matrix(NA, num_genes, 1)
        rownames(pcc_mat) <- rownames(matrix)
        colnames(pcc_mat) <- 'C3ar1'
        for(i in 1:num_genes){
                message('Computing #', i, ' of ', num_genes)
                #for(j in 1:num_genes){
                        pcc_mat[i,] <- cor(matrix[i,], matrix['C3ar1',], method = method.use)
                #}
        }
        return(pcc_mat)
}
input_mat <- as.matrix(mac.srt@assays$SCT@data[mac.srt@assays$SCT@var.features, mac.srt$sample_pub == 'YAP5SA'])
pcc_yap <- GenePCC(matrix = input_mat, method.use = 'pearson')
pcc_yap.df <- as.data.frame(pcc_yap)
pcc_yap.df$gene <- rownames(pcc_yap)
pcc_yap.df$Z <- as.vector(scale(pcc_yap.df$C3ar1))
pcc_yap.df <- pcc_yap.df[!is.na(pcc_yap.df$Z),]
pcc_yap.df <- pcc_yap.df[pcc_yap.df$Z>=1, ]
pcc_yap.df <- pcc_yap.df[order(pcc_yap.df$Z, decreasing = T),]
colnames(pcc_yap.df) <- c('Pearson_r_C3ar1', 'Gene', 'Z_score')
WriteCSV(pcc_yap.df, 'part90.06.c3ar1_correlated_gene')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  07. Small ST CM2 score and colocalization  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st2.srt <- readRDS('individual/part01.2022_JCardiovascDev_ESmall.srt.rds')

st2.srt <- AddModuleScore2(st2.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = F,
                           names = paste0('Score_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))
st2.srt <- AddModuleScore2(st2.srt, features = Yap5sa_CM_Mac_CF_marker, assay = 'SCT', return_z = T,
                           names = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)))
st2.srt <- AddModuleScore2(st2.srt, features = Yap_target[1:2], assay = 'SCT', return_z = F,
                           names = c('Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'))
st2.srt <- AddModuleScore2(st2.srt, features = Yap_target[1:2], assay = 'SCT', return_z = T,
                           names = c('Zscore_Y5SA_CM_Target1', 'Zscore_Y5SA_CM_Target2'))

st2.srt$tmp <- 'ALL'
RidgePlot(st2.srt, features = c("Score_Y5SA_CM1", "Score_Y5SA_CM2", "Score_Y5SA_C3ar1_Mac", "Score_Y5SA_C3_FB",
                                'Score_Y5SA_CM_Target1', 'Score_Y5SA_CM_Target2'), group.by = 'tmp',
          ncol = 1)
# shapiro.test(st2.srt$Score_Y5SA_CM1)$p.value
# shapiro.test(st2.srt$Score_Y5SA_CM2)$p.value
# shapiro.test(st2.srt$Score_Y5SA_C3ar1_Mac)$p.value
# shapiro.test(st2.srt$Score_Y5SA_C3_FB)$p.value
# shapiro.test(st2.srt$Score_Y5SA_C3_FB)$p.value

minvals <- rep(-1, 5)
maxvals <- rep(2, 5)
scale <- st2.srt@misc$spot_scale
p <- FeaturePlotST(srt = st2.srt,
                   features = c(paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)),
                                'Zscore_Y5SA_CM_Target1'),
                   title = c('CM1 score', 'CM2 score', 'C3ar1+ Mac score', 'C3+ FB score', 'Yap target score'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.2,
                   ncol = 4)
p[[16]] <- p[[16]] + RestoreLegend() + labs(color='Z score')
p[[20]] <- p[[20]] + RestoreLegend() + labs(color='Z score')
p
PlotPDF('07.1.feat_st.yap5sa_marker_and_yap_target_zscores', 12, 10)
print(p)
dev.off()

st2.srt$Colocal_CM1_C3ar1Mac_C3FB <-
        GetColocalProb(st2.srt, meta_features = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)[c(1, 3, 4)]),
                       return_logp = T)
st2.srt$Colocal_CM2_C3ar1Mac_C3FB <-
        GetColocalProb(st2.srt, meta_features = paste0('Zscore_Y5SA_', names(Yap5sa_CM_Mac_CF_marker)[c(2, 3, 4)]),
                       return_logp = T)

minvals <- rep(0, 4)
maxvals <- rep(1, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_C3ar1Mac_C3FB',
                                'Colocal_CM2_C3ar1Mac_C3FB'),
                   title = c('CM1, C3ar1+ Mac, C3+ FB',
                             'CM2, C3ar1+ Mac, C3+ FB'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.5,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 1), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')

PlotPDF('07.2.feat_st.yap5sa_marker_colocolization_on_small_st', 12, 5)
print(p)
dev.off()


st2.srt$Psig_CM1 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM2 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM1[st2.srt$Colocal_CM1_C3ar1Mac_C3FB < -log10(0.1)] <- 'Not sig.'
st2.srt$Psig_CM2[st2.srt$Colocal_CM2_C3ar1Mac_C3FB < -log10(0.1)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2, st2.srt$Sample)))
df$CM <- rep(c('CM1', 'CM2'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.1',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~CM) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('07.3.bar.yap5sa_marker_colocolization_on_olson_st', 4, 5)
print(p)
dev.off()


st2.srt$Colocal_CM1_YapTarget <-
        GetColocalProb(st2.srt, meta_features = c('Score_Y5SA_CM1', 'Score_Y5SA_CM_Target1'),
                       return_logp = T)
st2.srt$Colocal_CM2_YapTarget <-
        GetColocalProb(st2.srt, meta_features = c('Score_Y5SA_CM2', 'Score_Y5SA_CM_Target1'),
                       return_logp = T)

minvals <- rep(0, 4)
maxvals <- rep(1.5, 4)
p <- FeaturePlotST(srt = st2.srt,
                   features = c('Colocal_CM1_YapTarget',
                                'Colocal_CM2_YapTarget'),
                   title = c('CM1, YapHigh',
                             'CM2, YapHigh'),
                   minvals = minvals, maxvals = maxvals, pt.sizes = scale*0.5,
                   ncol = 4) &
        scale_color_distiller(palette = 'Reds', limits = c(0, 1.5), direction = 0)
p[[4]] <- p[[4]] + RestoreLegend() + labs(color='-Log10(p)')
p[[8]] <- p[[8]] + RestoreLegend() + labs(color='-Log10(p)')
p
PlotPDF('07.4.feat_st.yap5sa_marker_yap_target_colocolization_on_olson_st', 12, 5)
print(p)
dev.off()

st2.srt$Psig_CM1 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM2 <- factor('p<0.1', levels = c('Not sig.', 'p<0.1'))
st2.srt$Psig_CM1[st2.srt$Colocal_CM1_YapTarget < -log10(0.1)] <- 'Not sig.'
st2.srt$Psig_CM2[st2.srt$Colocal_CM2_YapTarget < -log10(0.1)] <- 'Not sig.'
df <- rbind(as.data.frame(Table(st2.srt$Psig_CM1, st2.srt$Sample)),
            as.data.frame(Table(st2.srt$Psig_CM2, st2.srt$Sample)))
df$CM <- rep(c('CM1', 'CM2'), each = 2*LU(st2.srt$Sample))
df$Fraction = df$Freq/rep(rep(as.vector(table(st2.srt$Sample)), each = 2), 2)
p <- ggplot(df[df$Var1=='p<0.1',]) +
        geom_bar(aes(x = Var2, y = Fraction, fill = Var1), stat = 'identity') +
        scale_fill_manual(values = c('firebrick3')) +
        facet_wrap(~CM) +
        theme_classic()+
        RotatedAxis() +
        labs(y = 'Fraction of spots', x = '', fill = 'Colocalization')
p
PlotPDF('07.5.bar.yap5sa_marker__yap_target_colocolization_on_olson_st', 4, 5)
print(p)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  08. CD45+ Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc.srt <- readRDS('individual/part01.aav_y5sa_cd45_immune.srt.rds')
p1 <- DimPlot2(sc.srt, reduction = 'sub_hmn_umap', group.by = 'sample_pub', cols = mycol_10, pt.size = 0.1) +
        labs(x = 'UMAP1', y = 'UMAP2', title = 'Myeloid Cells')
p1

PlotPDF('08.1.dim.cd45_myeloid_color_by_genotype', 5, 4)
print(p1)
dev.off()

sc.srt <- sc.srt[, sc.srt$Cell_type != 'DC']
sc.srt$sample_pub <- revalue(sc.srt$sample_pub, replace = c('WT' = 'Control'))
p2.1 <- CountCellBarPlot(sc.srt, group.var = 'sample_pub', stack.var = 'Cell_state',
                         percentage = T, stack.color = mycol_10)$plot +
        NoLegend() +
        labs(x = '', y = 'Fraction')
p2.2 <- CountCellBarPlot(sc.srt, group.var = 'sample_pub', stack.var = 'Cell_state',
                         percentage = F, stack.color = mycol_10)$plot +
        labs(x = '', y = 'Cell count' , fill  = 'Myeloid cell states')
p2.1 + p2.2
PlotPDF('08.2.bar.cd45_myeloid_composition_by_genotype', 8, 4)
print(p2.1 + p2.2)
dev.off()

sc.srt <- readRDS('individual/part01.aav_y5sa_cd45_immune.srt.rds')
genes <- read_excel('external/sciimmunol.abf7777_tables_s1_to_s9/sciimmunol.abf7777_table_s1.xlsx',
                    sheet = 4, col_names = T, skip = 2) ## These genes are from the E. Slava's Science Immun Paper
gl <- split(genes$gene, genes$cluster)
str(gl)

sc.srt <- AddModuleScore2(sc.srt, features = gl, names = paste0('Score_SE_', c('CCR2', 'MHC2', 'TLF')),
                          return_z = T, nbin = 20)
sc.srt <- AddModuleScore2(sc.srt, features = list(Yap5sa_CM_Mac_CF_marker$C3ar1_Mac),
                          names = 'Score_C3ar1Mac',
                          return_z = T, nbin = 20)
FeatureScatter(sc.srt, 'Score_SE_TLF', 'Score_SE_CCR2', cols = mycol_10, group.by = 'SE_cat')

p <- FeaturePlot2(sc.srt, features = c(paste0('Score_SE_', c('CCR2', 'MHC2', 'TLF')), 'Score_C3ar1Mac'),
                  reduction = 'sub_hmn_umap', min.cutoff = -1, max.cutoff = 2.5)
p
PlotPDF('08.3.feat.epelman_score_cd45_sc_mye', 8, 8)
p
dev.off()

# sc.srt <- sc.srt[, sc.srt$Cell_type %in% c('Mac', 'Mono')]
sc.srt$C3ar1_hi <- 'C3ar1_lo'
sc.srt$C3ar1_hi[sc.srt$Score_C3ar1Mac >= 0.5] <- 'C3ar1_hi'
sc.srt$SE_cat <- NA
sc.srt$SE_cat[sc.srt$Score_SE_CCR2 >= 1] <- 'CCR2'
sc.srt$SE_cat[sc.srt$Score_SE_MHC2 >= 1] <- 'MHC2'
sc.srt$SE_cat[sc.srt$Score_SE_TLF >= 1] <- 'TLF'
DimPlot2(sc.srt, group.by = c('SE_cat', 'C3ar1_hi'), reduction = 'sub_hmn_umap')
p1 <- CountCellBarPlot(sc.srt, group.var = 'SE_cat', stack.var = 'C3ar1_hi', percentage = T)$plot
p2 <- CountCellBarPlot(sc.srt[,!is.na(sc.srt$SE_cat)], group.var = 'C3ar1_hi', stack.var = 'SE_cat', percentage = T)$plot
p1|p2
PlotPDF('08.4.bar.epelman_vs_c3ar1', 6, 3)
p1|p2
dev.off()
total <- p2$data$Count[p2$data$GroupVar == 'C3ar1_hi']
round(total*100/sum(total), 2)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  09. Drop-seq CM Clustering  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cm.srt <- readRDS('external/fransico_data/reference_scrnaseq/Reference_combined_compendium_&_yap5sa_CMs_v2.rds')
DimPlot2(cm.srt, group.by = 'celltype_lvl2')
DimPlot2(cm.srt, group.by = 'pub')

cm.srt <- cm.srt[, cm.srt$celltype_lvl2 == 'cardiomyocyte']


DimPlot2(cm.srt, group.by = 'celltype', reduction = 'umap')
FeaturePlot2(cm.srt, features = c('Igf1r', 'Itgb1', 'Tnfrsf12a'), min.cutoff = 'q1', max.cutoff = 'q99')

cm.srt <- RunUMAP(cm.srt, dims = 1:50, reduction = 'harmony', n.neighbors = 200, min.dist = 2)
DimPlot2(cm.srt, group.by = 'celltype', reduction = 'umap')

new_cm.srt <- readRDS('individual/part01.y5sa_dropseq_cm.srt.rds')
identical(Cells(new_cm.srt), Cells(cm.srt))
new_cm.srt@reductions$umap2 <- new_cm.srt@reductions$umap
new_cm.srt@reductions$umap2@cell.embeddings <- cm.srt@reductions$umap@cell.embeddings
colnames(new_cm.srt@reductions$umap2@cell.embeddings) <- c('UMAP2_1', 'UMAP2_2')
new_cm.srt@reductions$umap2@key <- 'UMAP2_'
DefaultAssay(new_cm.srt) <- 'alra'

p1 <- FeaturePlot2(new_cm.srt, features = c('Tnfrsf12a'),
             reduction = 'umap2', min.cutoff = 0, max.cutoff = 3, raster = T, pt.size = 1.5, ncol = 1)
p2 <- FeaturePlot2(new_cm.srt, features = c('Igf1r'),
                   reduction = 'umap2', min.cutoff = 0, max.cutoff = 2, raster = T, pt.size = 1.5, ncol = 1)
p3 <- FeaturePlot2(new_cm.srt, features = c('Itgb1'),
                   reduction = 'umap2', min.cutoff = 1.5, max.cutoff = 3, raster = T, pt.size = 1.5, ncol = 1)
p4 <- DimPlot2(new_cm.srt, group.by = 'CM_State', reduction = 'umap2', raster = T, pt.size = 1.5, cols = mycol_10[2:1]) +
        labs(x = '', y = '', title = '')
p <- wrap_plots(p4, p1, p2, p3, ncol = 2)
p
PlotPDF('09.1.umap.cm_reembeded', 6, 6)
print(p)
dev.off()

WriteCSV(data.frame(new_cm.srt@reductions$umap2@cell.embeddings), title = '09.cm_umap_coordinates')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
