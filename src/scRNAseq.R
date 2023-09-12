##  General functions and variable initialization for
##  scRNA-seq

####--------------------------------------------------------------------------------------------------------------------
suppressMessages(library('Seurat'))
suppressMessages(library('SeuratWrappers'))
suppressMessages(library('SeuratData'))
suppressMessages(library('SeuratDisk'))
suppressMessages(library('harmony'))
# suppressMessages(library('monocle3'))
suppressMessages(library('plotly'))
suppressMessages(library('cluster'))
suppressMessages(library('clusterProfiler'))
suppressMessages(library('ReactomePA'))
suppressMessages(library('org.Hs.eg.db'))
suppressMessages(library('org.Mm.eg.db'))
suppressMessages(library('future')) ## for parallel processing
# suppressMessages(library('parallelDist')) ## for faster dist computation
## Below are SeuratWrappers functions
suppressMessages(library('schex'))
suppressMessages(library('Nebulosa')) ## kernel density estimation of expression pattern
suppressMessages(library('tricycle')) ## Circular visualization of cell cycle
suppressMessages(library('CIPR')) ## Reference based annotation of clusters
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
# cc.genes.mouse <- readRDS('~/Documents/Bioinformatics/r/src/cc_genes_mouse.rds')
s.genes <- c(
        'Gmnn', 'Rad51', 'Cdca7', 'Prim1', 'Mcm7', 'Dscc1', 'Fen1', 'Cenpu', 'Slbp', 'Ccne2', 'Mcm4', 'Polr1b',
        'Rad51ap1', 'Tyms', 'Rrm2', 'Wdr76', 'Casp8ap2', 'Pcna', 'Usp1', 'Chaf1b', 'Hells', 'Uhrf1', 'Nasp',
        'Tipin', 'Clspn', 'Cdc45', 'Ung', 'Rrm1', 'Ubr7', 'Rfc2', 'Pola1', 'Blm', 'Mcm5', 'Dtl', 'E2f8',
        'Cdc6', 'Mrpl36', 'Mcm6', 'Exo1', 'Msh2', 'Gins2'
        )
g2m.genes <- c(
        'Tacc3', 'Cdk1', 'Smc4', 'Tmpo', 'Ckap2l', 'Cks2', 'Mki67', 'Nusap1', 'Ect2', 'Cks1b', 'Cdca3', 'Ckap5',
        'Cdc25c', 'Top2a', 'Kif11', 'Cdca2', 'Kif20b', 'Ube2c', 'Ncapd2', 'Tpx2', 'Dlgap5', 'Aurka', 'G2e3',
        'Pimreg', 'Lbr', 'Ttk', 'Cdca8', 'Tubb4b', 'Cks1brt', 'Cenpe', 'Anln', 'Cbx5', 'Ctcf', 'Ccnb2', 'Aurkb',
        'Anp32e', 'Ndc80', 'Kif2c', 'Hmgb2', 'Gas2l3', 'Cenpf', 'Gtse1', 'Rangap1', 'Hjurp', 'Nuf2', 'Bub1',
        'Ckap2', 'Cdc20', 'Hmmr', 'Psrc1', 'Nek2', 'Birc5', 'Cenpa', 'Kif23'
        )
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
MeanExp <- function(srt_obj, assay = 'RNA', slot = 'data', pattern = NULL, features = NULL){
        if(!is.null(pattern)){
                mean_exp <- colMeans(GetAssayData(srt_obj, slot = slot)[grep(pattern = pattern,
                                                                             rownames(all.srt), value = T), ])}
        if(!is.null(features)){
                mean_exp <- colMeans(GetAssayData(srt_obj, slot = slot)[features, ])}
        return(mean_exp)}

FindDimNumber <- function(srt_obj, var.toal = 0.95, reduction = 'pca'){
        if(is.null(srt_obj@reductions[[reduction]])){
                cat("Reduction", reduction, "not found!")
                return(NULL)
        }
        tmp.var <- (srt_obj@reductions[[reduction]]@stdev)^2
        var.cut <- var.toal*sum(tmp.var)
        dimNum = 0
        var.sum = 0
        while(var.sum < var.cut){
                dimNum = dimNum + 1
                var.sum <- var.sum + tmp.var[dimNum]
        }
        return(dimNum)}

PrintFeatureLoading <- function(seruat.obj, reduction = 'pca', top_n_components = 20, top_features = 20){
        feat.list <- list()
        x <- seruat.obj@reductions[[reduction]]
        for (i in 1:top_n_components){
                feat.list[[i]] <- (names(x@feature.loadings[, i][order(x@feature.loadings[, i], decreasing = T)][1:top_features]))}
        return(feat.list)
}

DimPlot2 <- function(srt_obj,  group.by = NULL, ...){
        p <- DimPlot(object = srt_obj, shuffle = T, group.by = group.by, ...) &
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
        return(p)}

DimPlotSplit <- function(srt_obj, split_by, split_order, ncol = round(length(split_order)^0.5),
                         reduction = 'umap',  cols.bg = 'grey80', cols.highlight = 'red2', pt.size = 0.1,
                         sizes.highlight = 0.1, ...){
        p_list <- list()
        for(i in seq_along(split_order)){
                p_list[[i]] <- DimPlot(object = srt_obj, reduction = reduction,
                                       cells.highlight = Cells(srt_obj)[srt_obj[[split_by]][, 1] == split_order[i]],
                                       cols = cols.bg, cols.highlight = cols.highlight[i],
                                       pt.size = pt.size, sizes.highlight = sizes.highlight, ...) +
                        labs(title = split_order[i]) +
                        theme(aspect.ratio = 1,
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank()) +
                        NoLegend()
        }
        p <- wrap_plots(p_list, ncol = ncol)
        return(p)}

DimPlot3D <- function(srt_obj, reduction = 'umap', group.by, ...){
        data_df <- data.frame(Embeddings(srt_obj, reduction = reduction))
        if(ncol(data_df) < 3){stop('3D embedding not fount!')}
        colnames(data_df) <- c('data_dim_1', 'data_dim_2', 'data_dim_3')
        data_df[, 'group.by'] <- srt_obj[[group.by]]
        p <- plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,  z = ~data_dim_3, color = ~group.by,
                     type = "scatter3d", mode = "markers", size = I(1.5),
                     ...)
        return(p)
}


FeaturePlot2 <- function(srt_obj, features, pt.size = 0.1, cols = mycol_RYB[2:19], order = T,
                         recolor_subplot = NULL, recolor = NULL,
                         ncol = ceiling(L(features)^0.5), reduction = tail(names(srt_obj@reductions), 1), ...){
        p <- FeaturePlot(srt_obj, features = features, pt.size = pt.size, order = order,
                         ncol = ncol, reduction = reduction, ...) &
                scale_colour_gradientn(colours = cols) &
                #plot_layout(guides = "collect") &
                theme(aspect.ratio = 1,
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank())
        if(!is.null(recolor_subplot)){
                p[[recolor_subplot]] <- p[[recolor_subplot]] +
                        scale_colour_gradientn(colours = recolor)
        }
        return(p)}
FeaturePlot3 <- function(srt_obj, features, pt.size = 0.1, cols = mycol_Spec[90:5], joint = F,
                         adjust = 1,
                         recolor_subplot = NULL, recolor = NULL,
                         ncol = ceiling(L(features)^0.5),
                         reduction = tail(names(srt_obj@reductions), 1),
                         rescale_neg = F, ...){
        tmp <- srt_obj
        if(rescale_neg){
                for(i in 1:L(features)){
                        tmp@meta.data[, features[i]] <- Range01(tmp@meta.data[, features[i]])
                }
        }
        plist <- plot_density(tmp, features = features, size = pt.size, combine = F,
                              reduction = reduction, joint = joint, adjust = adjust, method = 'wkde', ...)
        p <- wrap_plots(plist, ncol = ncol) &
                scale_colour_gradientn(colours = cols) &
                theme(aspect.ratio = 1,
                      axis.title = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank())
        if(!is.null(recolor_subplot)){
                p[[recolor_subplot]] <- p[[recolor_subplot]] + scale_colour_gradientn(colours = recolor)
        }
        return(p)
}
FeaturePlot2_Dark <- function(srt_obj, pt.size = 0.2, cols = mycol_GG, order = T, ...){
        p <- FeaturePlot(srt_obj, pt.size = pt.size, cols = cols, order = order, ...) +
                plot_layout(guides = "collect") &
                theme(aspect.ratio = 1,
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"),
                      legend.background = element_rect(fill = "black"),
                      legend.box.background = element_rect(fill = "black", size = 0),
                      legend.key = element_rect(fill = "black", size = 0),
                      strip.background = element_rect(fill = "grey50", colour = NA),
                      axis.line.x = element_line(colour = "white"),
                      axis.line.y = element_line(colour = "white"),
                      panel.grid = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      text = element_text(colour = "white")
                      )
        return(p)}

VlnPlot2 <- function(srt_obj, pt.size = -1, ...){
        p <- VlnPlot(srt_obj, pt.size = pt.size, ...) &
                theme(aspect.ratio = 1,
                      axis.title = element_blank(), axis.title.x.bottom = element_blank())
        return(p)
}

BoxPlot <- function(srt_obj, feature, group.by, cols = mycol_10, split.by = NULL, split.ncol = NULL, ...){
        p <- VlnPlot(srt_obj, pt.size = -1, features = feature, group.by = group.by, ...)
        df <- p$data
        df$BOXPLOT_TMP <- df[,1]
        if(!is.null(split.by)){
                df$BOXPLOT_SPLIT <- srt_obj@meta.data[, split.by]
        }
        # levels <- levels(srt_obj@meta.data[, feature])
        # if(!is.null(levels)){df$ident <- factor(df$ident, levels = levels)}
        p <- ggplot(df) +
                geom_boxplot(aes(x = ident, y = BOXPLOT_TMP, fill = ident), outlier.shape = NA) +
                scale_fill_manual(values = cols) &
                theme_classic() &
                theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) &
                labs(fill = '', y = feature)
        if(!is.null(split.by)){
                p <- p + facet_wrap(~BOXPLOT_SPLIT, ncol = split.ncol, drop = T)
        }
        return(p)
}


DotPlot2 <- function(srt_obj, features, cols = c('grey85', mycol_BuGr[1]), scale.by = 'size', ...){
        p <- DotPlot(srt_obj, features = features, cols = cols, scale.by = scale.by, ...) +
                RotatedAxis() +
                scale_y_discrete(limits = rev) +
                theme(axis.title = element_blank(),
                      panel.grid.major = element_line(colour = 'grey85'))
        return(p)
}

HexPlot2 <- function(srt_obj, features, slot = 'data', aggregation = "mean", colors = rev(mycol_Spec), ...){
        p <- plot_hexbin_gene(sce = srt_obj, gene = features, type = slot, action = aggregation, ...) +
                guides(fill = guide_colourbar(barwidth = 0.3, barheight = 3)) +
                scale_fill_gradientn(colours = colors) +
                labs(title = paste0(features), fill = 'Cell Prop.') +
                theme(aspect.ratio = 1,
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank())
        return(p)
}

HexPlot2Multi <- function(srt_obj, features, slot = 'data', aggregation = "mean", colors = mycol_RnBw, ncol = 1, ...){
        p.list <- list()
        for(i in 1:L(features)){
                p.list[[i]] <- HexPlot2(srt_obj = srt_obj, features = features[i],
                                        slot = 'data', aggregation = "mean", colors = colors)
        }
        p <- wrap_plots(p.list, ncol = ncol)
        return(p)}

HexPlot2_Dark <- function(srt_obj, features, slot = 'data', aggregation = "mean", ...){
        p <- plot_hexbin_gene(sce = srt_obj, gene = features, type = slot, action = aggregation, ...) +
                guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5)) +
                scale_fill_gradient(low = mycol_GG[1], high = mycol_GG[100]) +
                labs(title = paste0(features,'-Positive Bins'), fill = 'Cell Prop.') +
                theme(aspect.ratio = 1,
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"),
                      legend.background = element_rect(fill = "black"),
                      legend.box.background = element_rect(fill = "black", size = 0),
                      legend.key = element_rect(fill = "black", size = 0),
                      strip.background = element_rect(fill = "grey50", colour = NA),
                      plot.title = element_text(colour = "white"), plot.subtitle = element_text(colour = "white"),
                      legend.title = element_text(colour = "white"),
                      legend.text = element_text(colour = "white"), strip.text = element_text(colour = "white"),
                      axis.line.x = element_line(colour = "white"),  axis.line.y = element_line(colour = "white"),
                      panel.grid = element_blank(), panel.grid.minor = element_blank(),
                      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
        return(p)}

MonocleDimPlot3D <- function(cds, dimension = "UMAP", col_group_by = 'orig.ident', colors, pt.size = 1.5,
                             angle = c(0, 0, 0), plot_title = NULL, width = 300*4, height = 300*4, ...){
        S_matrix <- reducedDims(cds)[[dimension]]
        data_df <- data.frame(S_matrix[, c(1, 2, 3)])
        colnames(data_df) <- c("data_dim_1", "data_dim_2", "data_dim_3")
        data_df$sample_name <- row.names(data_df)
        data_df <- as.data.frame(cbind(data_df, colData(cds)))
        data_df$col_group_by <- data_df[, col_group_by]
        font <- list(family = "Helvetica", size = 30, color = 'black')
        x_style <- list(title = '', showticklabels = F, titlefont = font)
        y_style <- list(title = '', showticklabels = F, titlefont = font)
        z_style <- list(title = '', showticklabels = F, titlefont = font)
        legend_style = list(font = font)
        p <- plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,  z = ~data_dim_3, color = ~col_group_by,
                     type = "scatter3d", mode = "markers", size = I(pt.size), colors = colors)
        p %>%
                layout(legend = legend_style, showlegend = T,
                       scene = list(xaxis = x_style, yaxis = y_style, zaxis = z_style,
                                    camera = list(eye = list(x = angle[1], y = angle[2], z = angle[3])))) %>%
                config(toImageButtonOptions = list(format = "png",
                                                   filename = plot_title, width = width, height = height))
        return(p)}

MonocleFeaturePlot3D <- function(cds, dimension = "UMAP", color_scale = rev(mycol_RnBw), pt.size = 1.5, expr_data,
                             angle = c(0, 0, 0), plot_title = NULL, width = 300*4, height = 300*4, ...){
        S_matrix <- reducedDims(cds)[[dimension]]
        data_df <- data.frame("data_dim_1" = S_matrix[, 1],
                              "data_dim_2" = S_matrix[, 2],
                              "data_dim_3" = S_matrix[, 3],
                              "gene" = expr_data[, 1])
        data_df$sample_name <- row.names(data_df)
        font <- list(family = "Helvetica", size = 30, color = 'black')
        x_style <- list(title = '', showticklabels = F, titlefont = font)
        y_style <- list(title = '', showticklabels = F, titlefont = font)
        z_style <- list(title = '', showticklabels = F, titlefont = font)
        legend_style = list(font = font)
        p <- plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,  z = ~data_dim_3, color = ~gene,
                     type = "scatter3d", mode = "markers", size = I(pt.size), colors = color_scale)
        p %>%
                layout(legend = legend_style, showlegend = T,
                       scene = list(xaxis = x_style, yaxis = y_style, zaxis = z_style,
                                    camera = list(eye = list(x = angle[1], y = angle[2], z = angle[3])))) %>%
                config(toImageButtonOptions = list(format = "png",
                                                   filename = plot_title, width = width, height = height))
        return(p)}

ChiSqP <- function(srt.obj, var1 = 'genotype', var2 = 'seurat_cluster'){
        celltypes = names(table(srt.obj[[var2]]))
        group1 = names(table(srt.obj[[var1]]))[1]
        group2 = names(table(srt.obj[[var1]]))[2]
        total_cell_num <- table(droplevels(srt.obj[[var1]]))
        count_mtx <- matrix(NA, L(celltypes), 2)
        rownames(count_mtx) <- celltypes
        colnames(count_mtx) <- names(table(srt.obj[[var1]]))

        n <- 0
        for(i in celltypes){
                n <- n+1
                count_mtx[n,] <- table(srt.obj[[var1]][srt.obj[[var2]] == i])[colnames(count_mtx)]
        }
        p_list <- c()
        for(celltype in celltypes){
                test_mtx <- rbind(c(sum(count_mtx[celltype, group1]), sum(count_mtx[celltype, group2])),
                                  c(sum(total_cell_num[group1])-sum(count_mtx[celltype, group1]),
                                    sum(total_cell_num[group2])-sum(count_mtx[celltype, group2])))
                p_list <- c(p_list, chisq.test(test_mtx)$p.value)
        }
        count_mtx <- cbind(count_mtx, 'pvalue' = p_list)
        return(as.data.frame(count_mtx))}

GetSilhouetteScore <- function(srt.obj, reduction = 'pca', dims, cluster_meta_col){
        ## srt.obj must be pre-processed and have nn graph computed
        if(ncol(srt.obj) <= 50e3){
                cells <- colnames(srt.obj)
        } else {
                cells <- sample(ncol(srt.obj), size = 50e3)
                message('Too many cells... Subsampling for 50k cells')
        }
        dist <- parDist(Embeddings(object = srt.obj, reduction = reduction)[cells, dims])
        clusters <- srt.obj@meta.data[cells, cluster_meta_col]
        sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist)
        sil_score.df <- data.frame('score' = sil[, 3], 'cluster' = clusters)
        sil_score.df <- tibble::as_tibble(sil_score.df[, ]) %>%
                dplyr::group_by(cluster) %>%
                dplyr::summarise("median_sil_per_clust" = median(score)) %>%
                tibble::add_column("res" = cluster_meta_col)
        return(sil_score.df)}

GetOptimalRes <- function(srt_obj, test_res, assay = 'RNA', reduction = 'pca', dims, ...){
        sil_score.df <- NULL
        for(i in 1:L(test_res)){
                message('Testing Resolution ', test_res[i], ' ...')
                srt_obj <- FindClusters(srt_obj, resolution = test_res[i])
                sil_score <- GetSilhouetteScore(srt_obj,
                                                reduction = reduction,
                                                dims = dims,
                                                cluster_meta_col = paste0(assay, '_snn_res.', test_res[i]))
                sil_score.df <- rbind(sil_score.df, sil_score)
                gc()
        }
        res_score.df <- sil_score.df %>%
                group_by(res) %>%
                summarise("median_sil_per_res" = median(median_sil_per_clust))
        optimal_res <- res_score.df$res[which.max(res_score.df$median_sil_per_res)]
        p <- ggplot(sil_score.df) +
                geom_boxplot(aes(x = res, y = median_sil_per_clust), fill = mycol_10[1:L(test_res)]) +
                geom_jitter(aes(x = res, y = median_sil_per_clust), width = 0.2) +
                WhiteBackground() +
                RotatedAxis() +
                labs(x = 'Resolution', y = 'Median Silhouette Scores', title = paste0('Optimal Resolution: ', optimal_res)) +
                theme(aspect.ratio = 2, axis.line = element_line())
        return(list(sil_score.df = sil_score.df,
                    res_score.df = res_score.df,
                    box.ggplot = p))}

DynamicQualityFilter <- function(srt_obj, filter_data, upper_c = 2, lower_c = 2){
        cutoff <- c()
        cells_toss <- c()
        if(is.null(upper_c)){
                cutoff[1] <- Inf
        } else{
                sd <- sd(srt_obj@meta.data[, filter_data])
                m <- mean(srt_obj@meta.data[, filter_data])
                cutoff[1] <- upper_c*sd + m ## Upper cutoff
        }
        if(is.null(lower_c)){
                cutoff[2] <- -Inf
        } else{
                sd <- sd(log10(srt_obj@meta.data[, filter_data]+1))
                m <- mean(log10(srt_obj@meta.data[, filter_data]+1))
                cutoff[2] <- 10^(m - lower_c*sd) ## Lower cutoff, built on log-transformed distribution
        }
        cells_toss <- Cells(srt_obj)[srt_obj@meta.data[, filter_data] > cutoff[1] |
                                             srt_obj@meta.data[, filter_data] < cutoff[2]]
        return(list('hi_cut' = cutoff[1],
                    'lo_cut' = cutoff[2],
                    'cells_toss' = cells_toss))}

ModuleEnrichment <- function(module_list, human_or_mouse){
        if(human_or_mouse == 'human'){
                OrgDb <- 'org.Hs.eg.db'
                organism <- 'human'
                kegg_org <- 'hsa'
        } else if (human_or_mouse == 'mouse') {
                OrgDb <- 'org.Mm.eg.db'
                organism <- 'mouse'
                kegg_org <- 'mmu'
        } else (stop('Wrong Species...'))
        module_geneID <- lapply(module_list,
                                function(gr) as.numeric(bitr(gr,
                                                             fromType = 'SYMBOL',
                                                             toType = 'ENTREZID',
                                                             OrgDb = OrgDb)$ENTREZID))
        modules.go <- list()
        modules.go_mf <- list()
        modules.go_cc <- list()
        modules.reactome <- list()
        modules.kegg <- list()
        for(i in 1:L(module_list)){
                modules.go[[i]] <- enrichGO(module_geneID[[i]], OrgDb = OrgDb, ont = 'BP', readable = T)@result
                modules.go_mf[[i]] <- enrichGO(module_geneID[[i]], OrgDb = OrgDb, ont = 'MF', readable = T)@result
                modules.go_cc[[i]] <- enrichGO(module_geneID[[i]], OrgDb = OrgDb, ont = 'CC', readable = T)@result
                modules.reactome[[i]] <- enrichPathway(module_geneID[[i]], organism = organism, readable = T)@result
                kegg <- enrichKEGG(module_geneID[[i]], organism = organism)
                modules.kegg[[i]] <-  setReadable(kegg, OrgDb = OrgDb, keyType = 'ENTREZID')@result
        }
        names(modules.go) <- paste0('GO_', names(module_list))
        names(modules.go_mf) <- paste0('GOMF_', names(module_list))
        names(modules.go_cc) <- paste0('GOCC_', names(module_list))
        names(modules.reactome) <- paste0('Reactome_', names(module_list))
        names(modules.kegg) <- paste0('KEGG_', names(module_list))
        return(list('GO' = modules.go, 'GOMF' = modules.go_mf, 'GOCC' = modules.go_cc, 
                    'Reactome' = modules.reactome, 'KEGG' = modules.kegg))
}

RemoveRiboMito <- function(genes, human_or_mouse, ribo = T, mito = T){
        if(human_or_mouse == 'human'){
                ribo_surfix <- '^RPL[[:digit:]]+'
                ribo_surfix2 <- '^RPS[[:digit:]]+'
                ribo_surfix3 <- '^RP[[:digit:]]+-'
                mito_surfix <- '^MT-'
        } else if (human_or_mouse == 'mouse') {
                ribo_surfix <- '^Rpl[[:digit:]]+'
                ribo_surfix2 <- '^Rps[[:digit:]]+'
                ribo_surfix3 <- '^Rp[[:digit:]]+-'
                mito_surfix <- '^mt-'
        } else (stop('Wrong Species...'))
        if(ribo & mito){
                genes <- grep(pattern = ribo_surfix, genes, value = T, invert = T)
                genes <- grep(pattern = ribo_surfix2, genes, value = T, invert = T)
                genes <- grep(pattern = ribo_surfix3, genes, value = T, invert = T)
                genes <- grep(pattern = mito_surfix, genes, value = T, invert = T)
        } else if(ribo){
                genes <- grep(pattern = ribo_surfix, genes, value = T, invert = T)
                genes <- grep(pattern = ribo_surfix2, genes, value = T, invert = T)
                genes <- grep(pattern = ribo_surfix3, genes, value = T, invert = T)
        } else if(mito){
                genes <- grep(pattern = mito_surfix, genes, value = T, invert = T)
        } else {stop('One of ribo or mito option must be TRUE...')}
        return(genes)
}
CountCellBarPlot <- function(srt_obj, group.var = NULL, stack.var = NULL,
                             stack.color = mycol_20, percentage = T, width = 0.7, ...){
        if(is.null(group.var)){
                srt_obj$CountCellBarPlot <- 1
                group.var <- 'CountCellBarPlot'
        }
        if(is.null(stack.var)){
                srt_obj$CountCellBarPlot <- 1
                stack.var <- 'CountCellBarPlot'
        }
        df <- srt_obj[[c(group.var, stack.var)]]
        df <- df %>% group_by(across(c(stack.var, group.var))) %>% summarise(count = n())
        colnames(df) <- c('StackVar', 'GroupVar', 'Count')
        if(percentage){
                p <- ggplot(df, aes(fill = StackVar, y = Count, x = GroupVar)) +
                        geom_bar(position = "fill", stat = "identity", width = width) +
                        scale_fill_manual(values = stack.color) +
                        labs(y = 'Fraction')
        } else {
                p <- ggplot(df, aes(fill = StackVar, y = Count, x = GroupVar)) +
                        geom_bar(position = "stack", stat = "identity", width = width) +
                        scale_fill_manual(values = stack.color) +
                        labs(y = 'Cell Count')
        }
        p <- p +
                labs(x = '', fill = '') +
                theme_minimal() +
                theme(
                        panel.grid.major.x =  element_blank(),
                        panel.grid.minor.x = element_blank()) +
                RotatedAxis()
        return(p)
}

SeuratToH5 <- function(srt, out_folder, embedding_label = NULL, save.meta = F, ...){
        dir.create(file.path(out_folder), showWarnings = FALSE)
        SaveH5Seurat(object = srt, filename = paste0(out_folder, '/out.h5seurat'))
        Convert(paste0(out_folder, '/out.h5seurat'), dest = "h5ad", overwrite = T)
        if(!is.null(embedding_label)){
                write.csv(srt[[embedding_label]]@cell.embeddings, paste0(out_folder, "/out.embeds.csv"), row.names = T)
        }
        if(save.meta){
                write.csv(srt@meta.data, paste0(out_folder, "/out.metadata.csv"), row.names = T)
        }
}
SaveH5ad <- function(srt, path, name, assay, raw_count_only = F, h5seurat_keep = F, verbose = T){
        for(i in 1:ncol(srt@meta.data)){srt@meta.data[, i] <- as.vector(srt@meta.data[, i])} ## convert factor to vectors
        srt <- DietSeurat(srt, scale.data = F,
                          assays = assay,
                          dimreducs = names(srt@reductions),
                          graphs = names(srt@graphs))
        if(raw_count_only){srt <- SetAssayData(srt, slot = 'data', new.data = GetAssayData(srt, slot = 'counts'))}
        if(verbose){message('Raw matrix:')
                print(GetAssayData(srt, slot = 'counts')[1:20, 1:10])
                message('Data matrix:')
                print(GetAssayData(srt, slot = 'data')[1:20, 1:10])
                message('Scaled Data matrix:')}
        if(sum(dim(GetAssayData(srt, assay = assay, slot = 'scale.data')))==0){message('No scaled data slot')} else{
                print(GetAssayData(srt, assay = assay, slot = 'scale.data')[1:20, 1:10])}
        SaveH5Seurat(object = srt, filename = paste0(path, '/', name, '.h5Seurat'), overwrite = T, verbose = verbose)
        Convert(paste0(path, '/', name, '.h5Seurat'), dest = "h5ad", assay = assay, overwrite = T, verbose = verbose)
        if(!h5seurat_keep){system(paste0('rm ', path, '/', name, '.h5Seurat'))}
}

GetOutlier <- function(vector, iqr_multiplier = 1.5){
        lowerq = quantile(vector)[2]
        upperq = quantile(vector)[4]
        iqr = upperq - lowerq
        threshold.upper = (iqr * iqr_multiplier) + upperq
        threshold.lower = lowerq - (iqr * iqr_multiplier)
        return(c(threshold.lower, threshold.upper))
}

MarkerHeatmap <- function(srt, marker.df, n_cells = 200, top = 10, rescale_all = T,
                          group.cols = mycol_40, disp.min = 0, raster = F){
        topmarker.df <- marker.df |> group_by(cluster) |> top_n(n = top, wt = avg_log2FC)
        cells <- c()
        for(i in levels(srt)){
                if(sum(srt@active.ident == i) >= n_cells){
                        cells <- c(cells, sample(Cells(srt)[srt@active.ident == i], size = n_cells))
                } else {
                        cells <- c(cells, sample(Cells(srt)[srt@active.ident == i], size = sum(srt@active.ident == i)))
                }
        }
        sub.srt <- DietSeurat(srt)
        sub.srt <- sub.srt[U(topmarker.df$gene), cells]
        if(rescale_all){sub.srt <- ScaleData(sub.srt, features = rownames(sub.srt), verbose = F)}
        p <- DoHeatmap(sub.srt, features = topmarker.df$gene,
                       group.colors = group.cols, disp.min = disp.min, size = 5, raster = raster) +
                #scale_fill_distiller(palette = 'GnBu', direction = 1)
                scale_fill_distiller(palette = 'Greys', direction = 1)
        return(p)
}
MarkerVolcano <- function(mk.df, label = T, min.log2FC = 0.25, max_p_adj = 0.01,
                          label_genes = NULL, max.overlaps = 10, line = F, ...){
        df <- mk.df
        df$gene <- rownames(df)
        clean_genes <- RemoveRiboMito(genes = df$gene, human_or_mouse = 'mouse', ribo = T, mito = T)
        df <- df[df$gene %in% clean_genes, ]
        df$log_p <- -log10(df$p_val_adj)
        df$Significance <- 'Not sig.'
        df$Significance[df$log_p >= -log10(max_p_adj) & abs(df$avg_log2FC) >= min.log2FC] <- paste0('Adjusted p < ', max_p_adj)
        df$Significance <- factor(df$Significance, levels = c(paste0('Adjusted p < ', max_p_adj), 'Not sig.'))
        if(is.null(label_genes)){df_label <- df[(df$log_p >= 3 & abs(df$avg_log2FC) >= 10) | df$log_p >= 20, ]
        }else{df_label <- df[df$gene %in% label_genes, ]}
        p <- ggplot(df) +
                geom_point(aes(x = avg_log2FC, y = log_p, color = Significance), size = 0.8) +
                geom_vline(xintercept = 0) +
                scale_color_manual(values = c(mycol_10[2], 'grey65')) +
                labs(x = 'Average Fold Change (log2)', y = 'Adjusted p-value (-log10)') +
                theme_classic() +
                theme(aspect.ratio = 1, panel.grid.major = element_line(colour = 'grey90'))
        if(label){
                p <- p +
                        ggrepel::geom_text_repel(data = df_label, aes(label = gene, x = avg_log2FC, y = log_p),
                                                 color = 'black', size = 3,
                                                 max.overlaps = max.overlaps, 
                                                 min.segment.length = ifelse(line, yes = 0, no = 0.3), ...)
        }
        return(p)
}

MappingHeatmap <- function(srt, que_var, ref_var, percentage = T, log10_scale = F, ref.disp.min = 0,
                           que_order = NULL, ref_order = NULL, cluster = F, center = T){
        df <- srt[[c(ref_var, que_var)]]
        mtx <- as.matrix(Table(df[,1], df[,2]))
        mtx <- mtx+1
        ## normalize by the total number of each query group
        ## i.e. n% of query x cells are mapped to y cell type in the reference
        if(percentage){
                for(i in 1:ncol(mtx)){mtx[,i] <- mtx[,i]/sum(mtx[,i])}
                log10_scale <- F
        }
        if(log10_scale){
                mtx <- log10(mtx+1)
        }
        mtx <- mtx[rowMaxs(mtx) >= ref.disp.min, ]
        mtx2 <- mtx
        mtx <- reshape2::melt(mtx)
        colnames(mtx) <- c('Reference', 'Query', 'Value')
        if(is.null(que_order)){ que_order <- str_sort(U(mtx$Query), numeric = T) }
        if(is.null(ref_order)){ ref_order <- str_sort(U(mtx$Reference), numeric = T) }
        if(!all(que_order %in% U(mtx$Query))){stop('Query order not found in seurat')}
        if(!all(ref_order %in% U(mtx$Reference))){stop('Reference order not found in seurat')}
        mtx <- mtx[mtx$Query %in% que_order & mtx$Reference %in% ref_order, ]
        mtx$Query <- factor(mtx$Query, levels = que_order)
        mtx$Query <- droplevels(mtx$Query)
        mtx$Reference <- factor(mtx$Reference, levels = ref_order)
        mtx$Reference <- droplevels(mtx$Reference)
        if(center){
                if(is.null(ref_order)){stop('ref_order needs to be specified for center = TRUE...')}
                mtx2 <- mtx2[ref_order, ]
                que_index <- rep(NA, ncol(mtx2))
                for(i in 1:ncol(mtx2)){que_index[i] <- which.max(mtx2[, i])}
                que_order <- colnames(mtx2)[order(que_index)]
                mtx$Query <- factor(mtx$Query, levels = que_order)
                mtx$Query <- droplevels(mtx$Query)
                mtx$Reference <- factor(mtx$Reference, levels = ref_order)
                mtx$Reference <- droplevels(mtx$Reference)
        } else if(cluster){
                if(is.null(que_order) | is.null(ref_order)){message('Warning: custom ordering disabled...')}
                mtx2 <- reshape2::dcast(mtx, Reference~Query, value.var = 'Value')
                rownames(mtx2) <- mtx2[,1]
                mtx2[,1] <- NULL
                mtx2
                row_order <- hclust(dist(mtx2))
                row_order <- row_order$labels[row_order$order]
                col_order <- hclust(dist(t(mtx2)))
                col_order <- col_order$labels[col_order$order]
                mtx$Query <- factor(mtx$Query, levels = col_order)
                mtx$Query <- droplevels(mtx$Query)
                mtx$Reference <- factor(mtx$Reference, levels = row_order)
                mtx$Reference <- droplevels(mtx$Reference)
        }
        p <- ggplot(data = mtx, aes(x = Reference, y = Query, fill = Value)) +
                geom_tile() +
                scale_fill_distiller(palette = 'RdYlBu') +
                labs(x = 'Reference groups', y = 'Query groups') +
                theme_classic() +
                theme(aspect.ratio = LU(que_order)/LU(ref_order)) +
                scale_y_discrete(limits = rev) +
                RotatedAxis()
        return(p)
}
MappingIdent <- function(srt, que_var, ref_var, cutoff = 0.1){
        df <- srt[[c(ref_var, que_var)]]
        mtx <- as.matrix(table(df[,1], df[,2]))
        mtx <- mtx+1
        match <- rep(NA, ncol(mtx))
        names(match) <- colnames(mtx)
        cutoff <- cutoff
        for(i in 1:ncol(mtx)){
                mtx[,i] <- mtx[,i]/sum(mtx[,i])
                match[names(match)[i]] <- names(which.max(mtx[,i]))
                if(mtx[names(which.max(mtx[,i])), i] < cutoff){match[names(match)[i]] <- 'NotFound'}
        }
        return(match)
}

DownsampleByMeta <- function(srt, meta_var, n, down_to_min_group = F, random = F, resample = F){
        cells <- split(Cells(srt), srt[[meta_var]])
        if(down_to_min_group){n = min(lengths(cells))}
        for(i in 1:L(cells)){
                if(L(cells[[i]]) > n){
                        if(random){
                                cells[[i]] <- sample(cells[[i]], size = n, replace = resample)
                        } else {
                                cells[[i]] <- cells[[i]][seq(1, L(cells[[i]]), round(L(cells[[i]])/n))]
                        }
                }
        }
        return(cells)
}

AddModuleScore2 <- function(srt, features, names, assay = NULL, return_z = T, ...){
        if(class(features) != 'list'){stop('Features are not in a list...')}
        if(L(features) != L(names)){stop('Features does not match names...')}
        orig.meta <- srt@meta.data
        orig.meta <- orig.meta[, ! colnames(orig.meta) %in% names] ## rewrite existing score
        srt@meta.data <- orig.meta
        if(is.null(assay)){assay <- DefaultAssay(srt)}
        srt <- AddModuleScore(srt, features = features, assay = assay, name = 'X___MODULE___', search = F, ...)
        new.meta <- srt@meta.data[, ! colnames(srt@meta.data) %in% colnames(orig.meta)]
        if(return_z){
                if(L(features) == 1){
                        srt@meta.data <- orig.meta
                        srt@meta.data[, names] <- as.vector(scale(new.meta))
                } else {
                        for(i in 1:L(features)){
                                colnames(new.meta)[colnames(new.meta) == paste0('X___MODULE___', i)] <- names[i]
                                new.meta[,i] <- scale(new.meta[,i])
                        }
                        srt@meta.data <- cbind(orig.meta, new.meta)
                }
        } else {
                if(L(features) == 1){
                        srt@meta.data <- orig.meta
                        srt@meta.data[, names] <- new.meta
                } else {
                        for(i in 1:L(features)){
                                colnames(new.meta)[colnames(new.meta) == paste0('X___MODULE___', i)] <- names[i]
                        }
                        srt@meta.data <- cbind(orig.meta, new.meta)
                }
        }
        return(srt)
}

EmbeddingMappingCorr <- function(ref_seurat, ref_assay = NULL, que_seurat, que_assay = NULL, features,
                                 method = 'pearson'){
        if(is.null(ref_assay)){ref_assay == ref_seurat@active.assay}
        if(is.null(que_assay)){que_assay == que_seurat@active.assay}
        ref_exp <- as.matrix(GetAssayData(ref_seurat, assay = ref_assay, slot = 'data')[features, ])
        que_exp <- as.matrix(GetAssayData(que_seurat, assay = que_assay, slot = 'data')[features, ])
        pcc_mat <- matrix(NA, ncol(ref_seurat), ncol(que_seurat), dimnames = list(Cells(ref_seurat), Cells(que_seurat)))
        for(i in 1:ncol(ref_seurat)){
                if(i %% 50 == 0){message(round(i*100/ncol(ref_seurat), digits = 0), '% data mapped ...')}
                for(j in 1:ncol(que_seurat)){
                        if(method == 'pearson') {pcc_mat[i, j] <- cor(ref_exp[, i], que_exp[, j])
                        } else if(method == 'cosine') {pcc_mat[i, j] <- lsa::cosine(ref_exp[, i], que_exp[, j])
                        } else {stop('Method must be pearson or cosine')}
                }
        }
        message('100% data mapped ...')
        return(pcc_mat)
}

EmbeddingMappingPlot <- function(ref_seurat, ref_cols = mycol_20, ref_group_by = 'ident',
                                 ref_pt_size = 0.5, pt_alpha = 1,
                                 que_seurat, que_cols = 'black', que_group_by = 'Query',
                                 que_pt_size = 1,
                                 ref_reduction = 'umap', k = 20, pcc_mat,
                                 jitter = 0) {
        if(nrow(pcc_mat) != ncol(ref_seurat) | ncol(pcc_mat) != ncol(que_seurat)){
                stop('pcc matrix does not match seurat')}

        knn_df <- matrix(NA, ncol(que_seurat), k)
        rownames(knn_df) <- Cells(que_seurat)
        for(i in 1:ncol(pcc_mat)){
                knn_df[i, ] <- head(rownames(pcc_mat)[order(pcc_mat[, i], decreasing = T)], k)
        }

        ref_emb <- as.data.frame(Embeddings(ref_seurat, reduction = ref_reduction))
        colnames(ref_emb) <- c('X', 'Y')
        if(ref_group_by == 'ident'){ref_emb$Group <- ref_seurat@active.ident
        } else {ref_emb$Group <- ref_seurat[[ref_group_by]][, 1]}
        if(is.null(levels(ref_emb$Group))){ref_emb$Group <- factor(ref_emb$Group)
        } else {ref_emb$Group <- droplevels(ref_emb$Group)}
        ref_emb$Size <- ref_pt_size

        que_emb <- data.frame(X = rep(NA, ncol(que_seurat)),
                              Y = rep(NA, ncol(que_seurat)))
        rownames(que_emb) <- Cells(que_seurat)
        for(i in 1:nrow(que_emb)){
                que_emb$X[i] <- median(ref_emb[knn_df[i, ], 'X'])
                que_emb$Y[i] <- median(ref_emb[knn_df[i, ], 'Y'])
        }
        if(que_group_by %in% colnames(que_seurat@meta.data)){que_emb$Group <- que_seurat[[que_group_by]][, 1]
        } else {que_emb$Group <- factor(que_group_by)}
        que_emb$Group <- droplevels(que_emb$Group)
        que_emb$Size <- que_pt_size

        ref_emb$Group <- factor(paste('Ref.', ref_emb$Group), levels = paste('Ref.', levels(ref_emb$Group)))
        que_emb$Group <- factor(paste('Test', que_emb$Group), levels = paste('Test', levels(que_emb$Group)))
        combine_emb <- rbind(ref_emb, que_emb)
        combine_emb$Group <- factor(combine_emb$Group, levels = c(levels(ref_emb$Group), levels(que_emb$Group)))
        cols <- c(ref_cols[1:L(levels(ref_emb$Group))], que_cols[1:L(levels(que_emb$Group))])
        p <- ggplot() +
                geom_point(data = combine_emb[combine_emb$Group %in% levels(ref_emb$Group), ],
                           aes(x = X, y = Y, color = Group, size = Size), alpha = pt_alpha) +
                geom_jitter(data = combine_emb[combine_emb$Group %in% levels(que_emb$Group), ],
                            aes(x = X, y = Y, color = Group, size = Size), alpha = pt_alpha,
                            width = jitter, height = jitter) +
                scale_color_manual(values = cols) +
                scale_size(breaks = c(min(ref_pt_size, que_pt_size), max(ref_pt_size, que_pt_size)), range = c(1, 2)) +
                theme_classic() +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
        return(p)
}

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

ExplainedVariance <- function(seurat, assay = 'RNA', reduction = 'pca', round = 2){
        mat <- GetAssayData(seurat, assay = assay, slot = 'scale.data')
        pca <- seurat[[reduction]]
        total_variance <- sum(matrixStats::rowVars(mat))
        eigValues <- (pca@stdev)^2
        varExplained <- eigValues / total_variance
        varExplainedPct <- round(varExplained*100, digits = round)
        return(varExplainedPct)
}

MarkerJitterPlot <- function(mk.df, cols = 'GnBu', col.direction = 0, jitter.width = 0.3, pt.size = 0.5){
        mk.df$TMP____logp <- -log10(mk.df$p_val_adj)
        mk.df$TMP____logp[mk.df$TMP____logp > 300] <- 300
        ggplot(mk.df, aes(x = cluster, y = avg_log2FC)) +
                geom_jitter(width = jitter.width, aes(color = TMP____logp), size = pt.size) +
                scale_color_distiller(palette = cols, direction = col.direction) +
                labs(x = 'Clusters', y = 'Mean Log2FC', color = '-Log10 Adj.P') +
                theme_classic() +
                theme(aspect.ratio = 1)
}

RunCellChat <- function(srt, group.by, CellChatDB = CellChatDB.mouse, LR.type = 'all', species = 'mouse'){
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
        cch <- computeCommunProb(cch, type = "truncatedMean", trim = 0.05)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}

AlignmentScore <- function(seurat, batch = 'sample', reduction = 'umap', k = NULL, downsample = 'max', iter = 100){
        x <- rep(0, ncol(seurat))
        names(x) <- Cells(seurat)
        y <- as.vector(seurat[[batch]][,1])
        if(is.null(k)){k <- round(ncol(seurat)/100)}
        N <- length(unique(y))
        knn <- FNN::get.knn(Embeddings(seurat, reduction), k = k, algorithm = "kd_tree")
        for(i in 1:ncol(seurat)) {
                q_ident <- as.character((y))[i]
                x[i] <- table(y[knn[[1]][i,]])[q_ident]
        }
        align_score <- rep(NA, iter)
        for(j in 1:iter){
                if(downsample == 'max') {
                        cell_dnSample <- DownsampleByMeta(seurat, meta_var = batch, down_to_min_group = T, random = T)
                } else {
                        cell_dnSample <- DownsampleByMeta(seurat, meta_var = batch, n = downsample, random = T)
                }
                x_dnSample <- x[unlist(cell_dnSample)]
                x_dnSample <- x_dnSample[!is.na(x_dnSample)]
                align_score[j] <- 1-(mean(x_dnSample)-k/N)/(k-k/N)
        }
        return(align_score)
}

DropMetaLevels <- function(seurat) {
        meta <- seurat@meta.data
        for(i in 1:ncol(meta)) {
                if(class(meta[, i])[1] == 'factor') {
                        meta[, i] <- droplevels(meta[, i])
                        message('Levels dropped for ', colnames(meta)[i])
                }
        }
        seurat@meta.data <- meta
        return(seurat)
}

SearchCellChat <- function(gene, species = 'human', search_term = 'ligand') {
        if(species == 'human') {
                database <- CellChatDB.human
        } else if(species == 'mouse') {
                database <- CellChatDB.mouse
        } else {stop('species must be human or mouse')}
        if(search_term == 'ligand') {
                terms <- database$interaction$ligand
        } else if(search_term == 'receptor') {
                terms <- database$interaction$receptor
        } else if(search_term == 'pathway') {
                terms <- database$interaction$pathway_name
        } else {stop('search_term must be ligand, receptor or pathway')}
        result <- database$interaction[grep(gene, terms), ]
        return(result)
}
