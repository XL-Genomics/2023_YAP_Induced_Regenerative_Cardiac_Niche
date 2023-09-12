##  General functions and variable initialization for
##  ST RNA-seq

####--------------------------------------------------------------------------------------------------------------------
suppressMessages(library('SPOTlight'))
suppressMessages(library('igraph'))
####--------------------------------------------------------------------------------------------------------------------
DimPlotST <-  function(srt, group.by, pt.sizes, cols = mycol_20, legend = 1,
                       ncol = ncol, section.var = 'Sample', na.col = 'grey'){
        names(cols) <- levels(srt@meta.data[, group.by])
        plist <- list()
        srt$TMP_____SECTION <- srt@meta.data[, section.var]
        for(i in 1:L(levels(srt$TMP_____SECTION))){
                sub.srt <- srt[, srt$TMP_____SECTION == levels(srt$TMP_____SECTION)[i]]
                x_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,1])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,1]))
                y_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,2])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,2]))
                asp <- y_scale/x_scale
                plist[[i]] <- DimPlot(sub.srt,
                                      group.by = group.by,
                                      reduction = 'spatial',
                                      pt.size = pt.sizes[i],
                                      raster = F,
                                      cols = cols[levels(droplevels(sub.srt@meta.data[, group.by]))],
                                      na.value = na.col) +
                        NoLegend() +
                        theme(aspect.ratio = asp,
                              axis.text = element_blank(),
                              axis.ticks = element_blank()) +
                        labs(title = levels(srt$TMP_____SECTION)[i], y = '', x = '')
        }
        for(x in legend){plist[[x]] <- plist[[x]] + RestoreLegend()} ## Add legends
        if(i > 1){p <- wrap_plots(plist, ncol = ncol)
        } else {p <- plist[[1]]}
        return(p)
}

FeaturePlotST <- function(srt, features, minvals, maxvals, pt.sizes, cols = 'Spectral', asp = NULL,
                          ncol = ncol, title = features, group.var = 'Sample'){
        plist <- list()
        srt$TMP <- srt@meta.data[, group.var]
        for(i in 1:L(levels(srt$TMP))){
                sub.srt <- srt[, srt$TMP == levels(srt$TMP)[i]]
                x_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,1])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,1]))
                y_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,2])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,2]))
                if(is.null(asp)){asp <- y_scale/x_scale}
                for(j in 1:L(features)){
                        idx <- (j-1)*L(levels(srt$TMP))+i
                        plist[[idx]] <- FeaturePlot(sub.srt,
                                                    order = F,
                                                    features = features[j],
                                                    reduction = 'spatial',
                                                    pt.size = pt.sizes[i],
                                                    min.cutoff = minvals[j],
                                                    max.cutoff = maxvals[j]) +
                                my.scale_colour_distiller(palette = cols,
                                                          limits = c(minvals[j], maxvals[j]),
                                                          aesthetics = 'color') +
                                theme(aspect.ratio = asp,
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank())
                        if(i == 1){
                                plist[[idx]] <- plist[[idx]] + labs(title = '', y = title[j], x = '')
                        } else {
                                plist[[idx]] <- plist[[idx]] + labs(title = '', x = '', y = '')
                        }
                        if(j == 1){
                                plist[[idx]] <- plist[[idx]] + labs(title = levels(srt$TMP)[i])
                        }
                }
        }
        p <- plist[[1]]
        for(x in 2:L(plist)){p <- p + plist[[x]]}

        p <- p + plot_layout(heights = 1, widths =  1,
                             ncol = ncol, nrow = L(features), byrow = T) &
                NoLegend() &
                theme(title = element_text(size = 10),
                      axis.title.y = element_text(face = 'bold', size = 10))
        return(p)
}

FeaturePlotST_Dark <- function(srt, features, minvals, maxvals, pt.sizes, asp = NULL,
                          ncol = ncol, title = features, group.var = 'Sample'){
        plist <- list()
        srt$TMP <- srt@meta.data[, group.var]
        for(i in 1:L(levels(srt$TMP))){
                sub.srt <- srt[, srt$TMP == levels(srt$TMP)[i]]
                x_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,1])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,1]))
                y_scale <- (max(sub.srt@reductions$spatial@cell.embeddings[,2])-
                                    min(sub.srt@reductions$spatial@cell.embeddings[,2]))
                if(is.null(asp)){asp <- y_scale/x_scale}
                for(j in 1:L(features)){
                        idx <- (j-1)*L(levels(srt$TMP))+i
                        plist[[idx]] <- FeaturePlot(sub.srt,
                                                    order = F,
                                                    features = features[j],
                                                    reduction = 'spatial',
                                                    pt.size = pt.sizes[i],
                                                    min.cutoff = minvals[j],
                                                    max.cutoff = maxvals[j]) +
                                scale_colour_gradient(low = 'grey20',
                                                      high = 'green',
                                                      limits = c(minvals[j], maxvals[j]),
                                                      aesthetics = 'color') +
                                theme(aspect.ratio = asp,
                                      text = element_text(color = 'white'),
                                      plot.background = element_rect(fill = "black"),
                                      panel.background = element_rect(fill = "black"),
                                      legend.background = element_rect(fill = "black"),
                                      legend.box.background = element_rect(fill = "black", size = 0),
                                      legend.key = element_rect(fill = "black", size = 0),
                                      strip.background = element_rect(fill = "grey50", colour = NA),
                                      plot.title = element_text(colour = "white"),
                                      plot.subtitle = element_text(colour = "white"),
                                      legend.title = element_text(colour = "white"),
                                      legend.text = element_text(colour = "white"),
                                      strip.text = element_text(colour = "white"),
                                      axis.line.x = element_line(colour = "white"),
                                      axis.line.y = element_line(colour = "white"),
                                      panel.grid = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_blank(),
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank(),
                                )
                        if(i == 1){
                                plist[[idx]] <- plist[[idx]] + labs(title = '', y = title[j], x = '')
                        } else {
                                plist[[idx]] <- plist[[idx]] + labs(title = '', x = '', y = '')
                        }
                        if(j == 1){
                                plist[[idx]] <- plist[[idx]] + labs(title = levels(srt$TMP)[i])
                        }
                }
        }
        p <- plist[[1]]
        for(x in 2:L(plist)){p <- p + plist[[x]]}

        p <- p + plot_layout(heights = 1, widths =  1,
                             ncol = ncol, nrow = L(features), byrow = T) &
                NoLegend() &
                theme(title = element_text(size = 10, color = 'white'),
                      axis.title.y = element_text(face = 'bold', size = 10),
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"))
        return(p)
}

GetColocalProb <- function(Seurat, meta_features){
        ## Make a Spots x Module_scores matrix
        n_feature <- length(meta_features) ## how many features (module_scores) are being tested for colocalization
        mtx <- matrix(NA, nrow = ncol(Seurat), ncol = n_feature)
        rownames(mtx) <- Cells(Seurat)
        colnames(mtx) <- meta_features

        for(i in 1:n_feature){
                scaled_scores <- scale(Seurat@meta.data[, meta_features[i]], center = T, scale = T) ## now mean = 0, sd = 1
                mtx[, i] <- 1 - pnorm(q = scaled_scores,
                                      mean = 0,
                                      sd = 1,
                                      lower.tail = T,
                                      log.p = F) ## replace scores with CDF p values for each module_score column
        }
        # coloc_pval <- -log10(rowProds(mtx[, ])) ## old method for combining p values
        coloc_pval <- -log10((rowMaxs(mtx[, ])^n_feature)) ## new method for combining p values
        return(coloc_pval)
}

PlotColocalProbDist <- function(srt, meta_features_list, split.by.factor, ncol = NULL, cols = mycol_10){
        mtx <- GetColocalProb(srt, meta_features = U(as.vector(unlist(meta_features_list))), return_prob_mat = T)
        df <- as.data.frame(mtx)
        rownames(df) <- Cells(srt)
        df$split <- srt@meta.data[, split.by.factor]
        df3 <- data.frame()
        split_levels <- levels(df$split)
        for(i in 1:L(split_levels)){
                df2 <- df[df$split == levels(df$split)[i], ]
                df4 <- data.frame(
                        Spots = rep((1:nrow(df2))/nrow(df2), L(meta_features_list)),
                        Split = rep(df2$split, L(meta_features_list)),
                        Group = rep(names(meta_features_list), each = nrow(df2)),
                        Prob = NA)
                for(j in 1:L(meta_features_list)){
                        df4$Prob[df4$Group == names(meta_features_list)[j]] <-
                                1-sort(rowProds(as.matrix(df2[, meta_features_list[[j]]])))}
                df3 <- rbind(df3, df4)
        }
        if(is.null(ncol)){ncol <- L(levels(df$split))}
        p <- ggplot(df3) +
                geom_line(aes(x = Spots, y = Prob, group = Group, colour = Group)) +
                scale_colour_manual(values = cols) +
                theme_classic() +
                theme(aspect.ratio = 1) +
                labs(y = 'Colocolization p value', x = 'Fraction of spots') +
                facet_wrap(~Split, scales = 'free', ncol = ncol)
        return(p)
}

PlotColocalPvalBar <- function(srt, group_factor = 'Sample', p_vals, p_vals_new_names = p_vals,
                               p_threshold = 0.05, cols = mycol_14, title = '', ylim = NULL){
        df <- srt@meta.data[, c(group_factor, p_vals)]
        df[,2:(L(p_vals)+1)] <- df[,2:(L(p_vals)+1)] > -log10(p_threshold)
        df2 <- tibble()
        for(i in 1:L(p_vals)){
                mtx <- Table(df[, group_factor], df[, p_vals[i]])
                total <- rowSums(mtx)
                if(!'TRUE' %in% colnames(mtx)){count <- 0
                }else {count <- mtx[, 'TRUE']}
                df2 <- bind_rows(df2,
                                 tibble(Group = factor(levels(srt@meta.data[, group_factor]),
                                                       levels = levels(srt@meta.data[, group_factor])),
                                        Fraction = count/total,
                                        X = p_vals_new_names[i]))
        }
        p <- ggplot(df2) +
                geom_bar(aes(x = Group, y = Fraction, fill = X), stat = 'identity', position = position_dodge()) +
                theme_classic() +
                labs(title = title, y = paste0('Fraction of p < ', p_threshold, ' spots'), x = '', fill = '') +
                RotatedAxis() +
                scale_fill_manual(values = cols)
        if(! is.null(ylim)){
                p <- p + scale_y_continuous(limits = ylim)
        }
        return(p)
}




# Finds the neighbors/surrounding cells
FindSpotNeighbors <- function(object.name) {
        samples <- names(object.name@images)
        n <- length(object.name@images)
        sample.n <- table(object.name$Sample)
        sample.cells <- split(colnames(object.name), object.name$Sample)

        results_list <- c()
        threshold_list <- c()

        message("***** Calculating Threshold *****")

        for(i in 1:n) {
                message("Now iterating through ", samples[i], "...")
                x <- object.name$Coord_x_slide[object.name$Sample == samples[i]]
                y <- object.name$Coord_y_slide[object.name$Sample == samples[i]]
                m <- cbind(x, y)
                matrix <- as.matrix(dist(m))

                distance <- c()

                for(j in 1:sample.n[i]) {
                        if(j %% 500 == 0) {
                                print(j)
                        }
                        spotName <- sample.cells[[samples[i]]][j]
                        mat <- matrix[, j]
                        sDist <- sort(mat)[-1]
                        topN <- head(sDist, 7)
                        distance[[j]] <- topN
                }

                names(distance) <- sample.cells[[samples[i]]]

                distance <- unlist(distance)

                distance_sort <- sort(unique(distance))
                distance_lag <- lag(distance_sort)
                distance_diff <- distance_sort - distance_lag
                distance_merged <- cbind(distance_sort, distance_lag, distance_diff)

                bound <- 2

                distance_bound <- ifelse(distance_diff > bound, distance_diff, NA)
                distance_bound <- distance_bound[!is.na(distance_bound)]

                distance_match <- match(distance_bound[1], distance_merged[, 3])
                threshold_list[[i]] <- mean(c(distance_sort[1], distance_sort[distance_match]))

        }

        message("***** Calculating Neighbors *****")

        for(i in 1:n) {
                message("Now iterating through ", samples[i], "...")
                x <- object.name$Coord_x_slide[object.name$Sample == samples[i]]
                y <- object.name$Coord_y_slide[object.name$Sample == samples[i]]
                m <- cbind(x, y)
                matrix <- as.matrix(dist(m))

                results <- c()

                for(j in 1:sample.n[i]) {
                        if(j %% 500 == 0) {
                                print(j)
                        }
                        spotName <- sample.cells[[samples[i]]][j]
                        mat <- matrix[, j]
                        sDist <- sort(mat)[-1]
                        top6 <- head(sDist, 6)
                        topN <- ifelse(top6 > threshold_list[[i]], NA, top6)
                        topN <- topN[!is.na(topN)]
                        topN <- names(topN)
                        results[[j]] <- topN
                }

                names(results) <- sample.cells[[samples[i]]]
                results_list[[i]] <- results

        }

        results_list <- unlist(results_list, recursive = FALSE)
        return(results_list)
}
