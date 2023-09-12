####--------------------------------------------------------------------------------------------------------------------
####  Mouse Heart scRNA-seq Compendium
####  2021-06-11 by Xiao LI (Texas Heart Institute, US)
####
####  Helper functions for data processing
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load Libraries  ####
####--------------------------------------------------------------------------------------------------------------------
suppressMessages(library('hdf5r'))
suppressMessages(library('data.table'))
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate global variables  ####
####--------------------------------------------------------------------------------------------------------------------
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
markers_lvl1_min <- markers_lvl1[seq(1, L(markers_lvl1), 2)]
markers_lvl2_early <- c('Epcam',  'Krt8',     # Mesoderm progen
                        'Nkx2-5', 'Isl1',     # Mesoderm
                        'Foxa2', 'Pax9',      # Endoderm
                        'Sox2',  'Wnt6')      # Ectoderm
markers_lvl2_immune <- c('Adgre1', 'Fcgr1',   # Mf, Mono
                         'Xcr1',   'Cd209a',  # Dc
                         'Cd79a',  'Ms4a1',   # B
                         'Cd3e',   'Nkg7',    # T+NK
                         'S100a9', 'S100a8')  # Granulocyte
markers_lvl2_bm <- c('Ccl5',   # ILC
                     'Vpreb1', # B
                     'Pf4',  # megakaryocytes
                     'Car1',   # erythrocytes
                     'Prss34',    # basophils
                     'Prg2')  # Geosinophils
markers_lvl2_tc <- c('Cd4',
                     'Trdc', #
                     'Cd8a', #
                     'Foxp3',  # Treg
                     'Ifngr1', # Treg, Th1
                     'Slamf6',   # Tfh
                     'Pdcd1',    # Cd8_exhausted
                     'Ccr7', # Cd8_mem (ccr7+Pdcd1+)
                     'Ctla4',
                     'Xcl1',
                     'Gzmb',
                     'Tox')  # Tfh
markers_lvl2_cm <- c('Nr2f1', 'Cav1', # Atrium
                     'Isl1', 'Tcn', # Outflow tract
                     'Myl2', 'Mpped2', # Ventrical
                     'Shox2', 'Pitx2'
)
## Colors
study_color <- clr_desaturate(shift = 0,
                              c(
                                      colorRampPalette(c("#2077b4", 'white'))(7)[1:5],
                                      colorRampPalette(c("#d62628", 'white'))(7)[1:5],
                                      colorRampPalette(c("#FF8C00", 'white'))(7)[1:5],
                                      colorRampPalette(c("#2ba32b", 'white'))(7)[1:5],
                                      colorRampPalette(c("#9466bd", 'white'))(7)[1:5],
                                      colorRampPalette(c("#e377c1", 'white'))(7)[1:5],
                                      colorRampPalette(c("#18bdcf", 'white'))(7)[1:5],
                                      colorRampPalette(c("#bcbd21", 'white'))(7)[1:5],
                                      colorRampPalette(c("#8b564c", 'white'))(7)[c(1,5)]
                              )
)
age_group_color <- ScaleColor(mycol_RYB, 11)[c(11:8, 4:1)]
condition_color <- c(paste0('grey', round(seq(40, 80, length.out = 9))), ScaleColor(mycol_BuGr, 11), mycol_40[5:8])
genotype_color <- mycol_20[c(3,4,6,14)]
condition_group_color <- c('grey70', rev(ScaleColor(mycol_BuGr, 6))[2:5], mycol_40[5:8])

options(scipen = 999)
options(future.globals.maxSize= 2000*1024^2)
set.seed(505)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate global functions  ####
####--------------------------------------------------------------------------------------------------------------------
InitiateProject <- function(Machine = 'Rivendell', Ver, Part, Catagory = 'mouse', Project_dir, Data_drive){
        if (Machine == 'Rivendell'){
                Data_dir <<- paste0('/Volumes/', Data_drive, '/project/', Project_dir, '/')
                File_dir <<- paste0('~/Documents/Bioinformatics/project/', Project_dir, '/')
        } else if (Machine == 'Gondor') {
                Data_dir <<- paste0('/Volumes/', Data_drive, '/project/', Project_dir, '/')
                File_dir <<- paste0('~/Documents/Bioinformatics/project/', Project_dir, '/')
        } else if (Machine == 'Moria') {
                Data_dir <<- paste0('/moria/', Project_dir, '/')
                File_dir <<- Data_dir
        } else if (Machine == 'bobbyd') {
                Data_dir <<- paste0('~/work/', Project_dir, '/')
                File_dir <<- Data_dir
        } else {stop('Machine not defined')}
        ####------------------------------------------------------------------------------------------------------------
        Rdata_dir <<- paste0(Data_dir, 'rdata/', Catagory, '_v', Ver)
        Docu_dir <<- paste0(File_dir, 'doc/', Catagory, '_v', Ver, '/')
        Plot_dir <<- paste0(File_dir, 'plot/', Catagory, '_v', Ver, '/', Part, '/')
        Meta_dir <<- paste0(File_dir, 'meta/', Catagory, '_v', Ver, '/', Part, '/')
        dir.create(file.path(Rdata_dir), showWarnings = F, recursive = T)
        setwd(Rdata_dir)
        dir.create(file.path('./individual'), showWarnings = F, recursive = T)
        dir.create(file.path('./integrated'), showWarnings = F, recursive = T)
        dir.create(file.path('./analysis'), showWarnings = F, recursive = T)
        dir.create(file.path('./external'), showWarnings = F, recursive = T)
        dir.create(file.path('./tmp'), showWarnings = F, recursive = T)
        dir.create(file.path(Plot_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Meta_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Docu_dir), showWarnings = F, recursive = T)
}
MakeScmultiSrt <- function(gex_filtered_h5_file, atac_frag_file){
        # load both modalities
        inputdata.10x <- Read10X_h5(gex_filtered_h5_file)
        # extract ATAC data
        atac_counts <- inputdata.10x$Peaks
        # only use peaks in standard chromosomes
        grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
        grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
        atac_counts <- atac_counts[as.vector(grange.use), ]
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevelsStyle(annotations) <- 'UCSC'
        genome(annotations) <- "mm10"
        chrom_assay <- CreateChromatinAssay(
                counts = atac_counts,
                sep = c(":", "-"),
                genome = 'mm10',
                fragments = atac_frag_file,
                min.cells = 1,
                min.features = 1,
                annotation = annotations
        )
        srt <- CreateSeuratObject(
                counts = chrom_assay,
                assay = "ATAC"
        )
        # extract RNA data
        rna_counts <- inputdata.10x$`Gene Expression`
        # Create Seurat object
        srt2 <- CreateSeuratObject(counts = rna_counts)
        return(list(srt, srt2))
}
MakeScmultiGexSrt <- function(gex_filtered_h5_file, atac_frag_file){
        # load both modalities
        inputdata.10x <- Read10X_h5(gex_filtered_h5_file)
        # extract ATAC data
        rna_counts <- inputdata.10x$
        # only use peaks in standard chromosomes
        grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
        grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
        atac_counts <- atac_counts[as.vector(grange.use), ]
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevelsStyle(annotations) <- 'UCSC'
        genome(annotations) <- "mm10"
        chrom_assay <- CreateChromatinAssay(
                counts = atac_counts,
                sep = c(":", "-"),
                genome = 'mm10',
                fragments = atac_frag_file,
                min.cells = 1,
                min.features = 1,
                annotation = annotations
        )
        srt <- CreateSeuratObject(
                counts = chrom_assay,
                assay = "ATAC"
        )
        return(srt)
}
MakeSrt <- function(mode, matrix_dir, do_decontX = F, sample, study, method, platform, protocol, processed, tissue,
                    enrichment, preparation, sex, age, genotype, genotype_full, condition, strain, replicate, group) {
        if (mode == '10x') {
                matrix <- paste0(matrix_dir, '/outs/filtered_feature_bc_matrix/')
                srt <- CreateSeuratObject(counts = Read10X(data.dir = matrix),
                                          min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'cellbender') {
                matrix <- paste0(matrix_dir, '/cellbender_filtered.h5')
                srt <- CreateSeuratObject(counts = ReadCB_h5(matrix), min.cells = 1, min.features = 1, project = study)
        }
        ## Use decontX for ambient removal:
        if (do_decontX) {
                dcx.mtx <- decontX(GetAssayData(srt, slot = 'counts', assay = 'RNA'), seed = 505, verbose = F)
                srt <- CreateSeuratObject(counts = round(dcx.mtx$decontXcounts),
                                          min.cells = 1, min.features = 1, project = study)
        }
        srt$sample <- sample
        srt$orig.name <- Cells(srt)
        srt$study <- study
        srt$method <- method
        srt$platform <- platform
        srt$protocol <- protocol
        srt$processed <- processed
        srt$tissue <- tissue
        srt$enrichment <- enrichment
        srt$preparation <- preparation
        srt$sex <- sex
        srt$age <- age
        srt$genotype <- genotype
        srt$genotype_full <- genotype_full
        srt$condition <- condition
        srt$strain <- strain
        srt$replicate <- replicate
        srt$group <- group
        srt <- RenameCells(srt, new.names = paste(srt$study,
                                                  srt$sample,
                                                  srt$orig.name,
                                                  sep = ':'), for.merge = F)
        srt <- PercentageFeatureSet(srt, pattern = '^mt-', col.name = 'pct_mito', assay = 'RNA')
        srt$pct_mito[is.nan(srt$pct_mito)] <- 0
        srt$pct_mito[is.na(srt$pct_mito)] <- 0
        return(srt)
}
MakeDataset <- function(study, study_id, sample_name, mode, matrix_dir, starting_sample = 1, do_decontX = F){
        srt.list <- list()
        for(i in 1:L(sample_name)) {
                sample_id = paste0(study_id, '_S', str_pad(starting_sample - 1 + i, 3, pad = '0'))
                print(paste('Processing sample:', sample_id))
                sample_meta_sub.df <- sample_meta.df[sample_meta.df$study == study &
                                                             sample_meta.df$sample_id == sample_id &
                                                             sample_meta.df$name_on_disk == sample_name[i], ]
                # print('Sample metadata found: ')
                # print(format(sample_meta_sub.df))
                srt.list[[i]] <- MakeSrt(mode = mode,
                                         matrix_dir = matrix_dir[i],
                                         do_decontX = do_decontX,
                                         study = study,
                                         sample = sample_id,
                                         method = sample_meta_sub.df$method,
                                         platform = sample_meta_sub.df$platform,
                                         protocol = sample_meta_sub.df$protocol,
                                         processed = sample_meta_sub.df$processed,
                                         tissue = sample_meta_sub.df$tissue,
                                         enrichment = sample_meta_sub.df$enrichment,
                                         preparation = sample_meta_sub.df$preparation,
                                         sex = sample_meta_sub.df$sex,
                                         age = sample_meta_sub.df$age,
                                         genotype = sample_meta_sub.df$genotype_s,
                                         genotype_full = sample_meta_sub.df$genotype_l,
                                         condition = sample_meta_sub.df$condition,
                                         strain = sample_meta_sub.df$strain,
                                         replicate = sample_meta_sub.df$replicate,
                                         group = paste(sample_meta_sub.df$study,
                                                       sample_meta_sub.df$age,
                                                       sample_meta_sub.df$genotype_s,
                                                       sample_meta_sub.df$tissue,
                                                       sample_meta_sub.df$enrichment,
                                                       sample_meta_sub.df$condition,
                                                       sample_meta_sub.df$replicate,
                                                       sep = '__')
                )
                print('Seurat generated...')
                # print(srt.list[[i]])
                # cat('\n_____________________________________________________\n')
        }
        if(L(srt.list) > 1) {
                merge.srt <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
        } else {
                merge.srt <- srt.list[[1]]
        }
        return(merge.srt)
}


ProcessSrt_std <- function(srt_obj, var.toal = 0.75, assay = 'RNA', do.umap = T, npcs = 50, ...) {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^mt-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        srt.out <- srt.out %>%
                FindVariableFeatures() %>%
                ScaleData(vars.to.regress = c(paste0( c('nFeature_', 'pct_mito_'), assay), 'S.Score', 'G2M.Score'),
                          ...) %>%
                RunPCA(verbose = F, seed.use = 505, npcs = npcs)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)}
        return(srt.out)
}

ProcessSrt_sct <- function(srt_obj, var.toal = 0.75, do.umap = T, assay = 'RNA') {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^mt-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        srt.out <- SCTransform(srt.out,
                               assay = assay,
                               method = "glmGamPoi",
                               seed.use = 505,
                               return.only.var.genes = F,
                               vars.to.regress = c(paste0( c('nFeature_', 'pct_mito_'), assay), 'S.Score', 'G2M.Score')
                               ) %>%
                RunPCA(verbose = F, seed.use = 505)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)}
        return(srt.out)
}

ProcessSrt_hmn <- function(srt_obj, haromize.by = 'sample', harmony_max_iter = 30, var.toal = 0.75, assay = 'RNA',
                           do.umap = T, harmony_early_stop = T, ...) {
        if(harmony_early_stop){
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   max.iter.harmony = harmony_max_iter,
                                   plot_convergence = T,
                                   assay.use = assay, ...)
        } else {
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                                   plot_convergence = T,
                                   max.iter.harmony = harmony_max_iter,
                                   assay.use = assay, ...)
        }
        if(do.umap){
                dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'harmony')
                srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505, reduction = 'harmony',
                                   reduction.name = 'hmn_umap', reduction.key = 'hmnumap_', ...)
        }
        return(srt.out)
}

ProcessSrt_clust <- function(srt_obj, resolution = 0.2, reduction = 'harmony', var.toal = 0.75, ...){
        dimN <- FindDimNumber(srt_obj = srt_obj, var.toal = var.toal, reduction = reduction)
        srt.out <- srt_obj %>%
                FindNeighbors(dims = 1:dimN, reduction = reduction, force.recalc = T) %>%
                FindClusters(resolution = resolution, ...)
}

IntersectSrts <- function(ref.srt, que.srt, ref_data.name, que_data.name) {
        ## Identify shared cells
        shared_cells <- intersect(Cells(ref.srt), Cells(que.srt))
        title <- paste0(ref_data.name, ' has: ', ncol(ref.srt), '... ',
                        que_data.name, ' has: ', ncol(que.srt), '... ',
                        'Shared cells: ',        L(shared_cells))
        ref.srt <- ref.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p1 <- DimPlot(ref.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(ref_data.name, ' has: ', ncol(ref.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        que.srt <- que.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p2 <- DimPlot(que.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(que_data.name, ' has: ', ncol(que.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        print(title)
        ## Create new intersect srt from query data
        que.srt <- que.srt[, shared_cells]
        print('New seurat object generated:')
        print(que.srt)

        ## Match meta data
        ref.meta <- ref.srt@meta.data
        que.meta <- ref.meta[Cells(que.srt),
                             c('pub', 'sample', 'orig.name', 'method', 'platform', 'protocol', 'processed',
                               'tissue', 'enrichment', 'preparation', 'sex', 'age', 'genotype_s', 'genotype_l',
                               'condition', 'strain', 'replicate', 'group', 'batch')]
        print(paste0(nrow(que.meta), ' cells found in reference metadata'))
        que.srt@meta.data <- cbind(que.srt@meta.data, que.meta)
        que.srt$processed <- 'CellBender'
        print(head(que.srt@meta.data, n = 3))
        return(list(p1, p2, que.srt))
}

## A Cellbender bug workaround:
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
        if (!requireNamespace('hdf5r', quietly = TRUE)) {
                stop("Please install hdf5r to read HDF5 files")
        }
        if (!file.exists(filename)) {
                stop("File not found")
        }
        infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
        genomes <- names(x = infile)
        output <- list()
        if (hdf5r::existsGroup(infile, 'matrix')) {
                # cellranger version 3
                message('CellRanger version 3+ format H5')
                if (use.names) {
                        feature_slot <- 'features/name'
                } else {
                        feature_slot <- 'features/id'
                }
        } else {
                message('CellRanger version 2 format H5')
                if (use.names) {
                        feature_slot <- 'gene_names'
                } else {
                        feature_slot <- 'genes'
                }
        }
        for (genome in genomes) {
                counts <- infile[[paste0(genome, '/data')]]
                indices <- infile[[paste0(genome, '/indices')]]
                indptr <- infile[[paste0(genome, '/indptr')]]
                shp <- infile[[paste0(genome, '/shape')]]
                features <- infile[[paste0(genome, '/', feature_slot)]][]
                barcodes <- infile[[paste0(genome, '/barcodes')]]
                sparse.mat <- sparseMatrix(
                        i = indices[] + 1,
                        p = indptr[],
                        x = as.numeric(x = counts[]),
                        dims = shp[],
                        giveCsparse = FALSE
                )
                if (unique.features) {
                        features <- make.unique(names = features)
                }
                rownames(x = sparse.mat) <- features
                colnames(x = sparse.mat) <- barcodes[]
                sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
                # Split v3 multimodal
                if (infile$exists(name = paste0(genome, '/features'))) {
                        types <- infile[[paste0(genome, '/features/feature_type')]][]
                        types.unique <- unique(x = types)
                        if (length(x = types.unique) > 1) {
                                message("Genome ",
                                        genome, "
                                        has multiple modalities, returning a list of matrices for this genome")
                                sparse.mat <- sapply(
                                        X = types.unique,
                                        FUN = function(x) {
                                                return(sparse.mat[which(x = types == x), ])
                                        },
                                        simplify = FALSE,
                                        USE.NAMES = TRUE
                                )
                        }
                }
                output[[genome]] <- sparse.mat
        }
        infile$close_all()
        if (length(x = output) == 1) {
                return(output[[genome]])
        } else{
                return(output)
        }
}

QuickCheck <- function(srt_obj, markers, ...) {
        p1 <- VlnPlot(srt_obj,
                       features = c('nFeature_RNA', 'nCount_RNA', 'pct_mito'),
                       group.by = 'batch', ncol = 3, pt.size = -1) &
                theme(aspect.ratio = 0.5)
        p1 <- wrap_plots(list((
                p1[[1]] + theme(axis.text.x = element_blank())),
                (p1[[2]] + theme(axis.text.x = element_blank())),
                p1[[3]]), ncol = 1)
        p2.1 <- DimPlot(srt_obj, group.by = "batch", reduction = 'hmn_umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p2.2 <- DimPlot(srt_obj, group.by = "batch", reduction = 'umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p3 <- FeaturePlot2(srt_obj,
                           reduction = 'hmn_umap',
                           features = intersect(markers, rownames(srt_obj)),
                           ncol = ceiling(L(intersect(markers, rownames(srt_obj))) / 4))
        return(list(p2.1, p2.2, p1, p3))
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

Subcluster <- function(plot_num, srt, celltype, dimN, resolution, features = NULL, fast_marker = F,
                       regress = c('nFeature_CBN', 'pct_mito_CBN')){
        ## Seurat processing
        message('\n', 'Processing Seurat...')
        srt <- srt |>
                FindVariableFeatures(verbose = F) |>
                ScaleData(verbose = F, vars.to.regress = regress) |>
                RunPCA(seed.use = 505, verbose = F) |>
                RunHarmony(group.by.vars = 'sample', max.iter.harmony = 50, assay.use = 'CBN', verbose = F) |>
                RunUMAP(reduction = 'harmony', dims = 1:dimN, verbose = F)
        srt <- DietSeurat(srt, dimreducs = c('umap', 'harmony'))
        gc()
        p <- wrap_plots(list(DimPlot2(srt, group.by = 'sample', raster = F, cols = mycol_10),
                             DimPlot2(srt, group.by = 'genotype', raster = F, cols = mycol_10)
                             ), ncol = 2)
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.meta'), 30, 20)
        print(p)
        dev.off()
        if(is.null(features)){ features <- markers_lvl1}
        p <- FeaturePlot2(srt, features = U(features), ncol = ceiling(LU(features)^0.5), reduction = 'umap')
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.feature'),
                ceiling(LU(features)^0.5)*2,
                floor(LU(features)^0.5)*2)
        print(p)
        dev.off()
        ## Subcluster evaluation
        message('\n', 'Clustering with ', dimN, ' components, at ', resolution, ' resolution...')
        srt <- srt |>
                FindNeighbors(reduction = 'harmony', dims = 1:dimN) |>
                FindClusters(resolution = resolution)
        for(k in resolution){
                Idents(srt) <- paste0('CBN_snn_res.', k)
                levels(srt) <- str_sort(levels(srt), numeric = T)
                if(fast_marker){
                        marker <- FindAllMarkers(srt, assay = 'CBN', only.pos = T, return.thresh = 0.0001, max.cells.per.ident = 5e3)
                } else {
                        marker <- FindAllMarkers(srt, assay = 'CBN', only.pos = T, return.thresh = 0.0001)
                }
                srt@misc$marker[[paste0('CBN_snn_res.', k)]] <- marker
                cols <- mycol_20
                if(LU(srt@active.ident)>20){cols <- mycol_40}
                p <- DimPlot2(srt, label = T, repel = T, raster = F, cols = cols) + labs(title = paste0('Louvain_', k))
                PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', k, 'res.umap'), 10, 10)
                print(p)
                dev.off()
                p <- MarkerHeatmap(srt, marker.df = marker, top = 20)
                PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', k, 'res.heatmap'),
                        20,
                        L(levels(srt))*20/8)
                print(p)
                dev.off()
                gc()
        }
        return(srt)
}

GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal, threshold = 0.20){
        srt_obj <- ProcessSrt_std(srt_obj, var.toal = dimN.var.toal, assay = 'CBN')
        ## Scrublet (run via reticulate)
        mtx <- srt_obj@assays$CBN@counts
        mtx <- t(mtx)
        scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
        rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                          n_prin_comps = 30L,
                                          min_counts = 2, min_cells = 3)
        rst[[2]] <- scrub_model$call_doublets(threshold = threshold) ## adjusted based on histogram
        sc_doublets <- Cells(srt_obj)[rst[[2]]]
        sc_singlets <- Cells(srt_obj)[!rst[[2]]]
        srt_obj$Scrublet_doublet <- 'Singlet'
        srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
        Scrublet <- rst[[1]]
        names(Scrublet) <- Cells(srt_obj)

        p2 <- DimPlotSplit(srt_obj, split_by = 'Scrublet_doublet', split_order = c('Singlet', 'Doublet'),
                           cols.highlight = mycol_14[c(2, 1)], ncol = 2)
        p2[[1]] <- p2[[1]] + labs(title = paste0('Srub Singlet: ', L(sc_singlets), ' Cells'))
        p2[[2]] <- p2[[2]] + labs(title = paste0('Srub Doublet: ', L(sc_doublets), ' Cells'))
        p <- wrap_plots(
                p2[[1]],
                p2[[2]],
                ncol = 2)
        return(list(
                sc_doublets,
                p,
                Scrublet,
                scrub_model
        ))
}
####--------------------------------------------------------------------------------------------------------------------
