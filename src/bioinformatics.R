##  General functions and variable initialization for
##  general statistics and bioinformatics

####--------------------------------------------------------------------------------------------------------------------
suppressMessages(library('RColorBrewer'))
suppressMessages(library('colorRamps'))
suppressMessages(library('viridisLite'))
suppressMessages(library('plyr'))
suppressMessages(library('prismatic'))
suppressMessages(library('tidyverse'))
suppressMessages(library('pheatmap'))
suppressMessages(library('patchwork'))
suppressMessages(library('DescTools')) ## For Winsorize
suppressMessages(library('Matrix'))
suppressMessages(library('eulerr'))
suppressMessages(library('readxl'))
suppressMessages(library('scales'))
suppressMessages(library('grid'))
suppressMessages(library('ggthemes'))
suppressMessages(library('reshape2'))
suppressMessages(library('ggbeeswarm'))
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
mycol_BuGr <- clr_darken(shift = 0, c('#0B486B', '#3B8686', '#79BD9A', '#A8DBA8', '#CFF09E')) # 1 way gradient
mycol_GG <- clr_darken(shift = 0, colorRampPalette(c("grey30", "green"))(100)) ## Dark grey -- bright green
mycol_RB <- clr_darken(shift = 0, colorRampPalette(c("#2358A3", "white", "#EE2F37"))(100)) ## Blue--White--Red
mycol_RG <- clr_darken(shift = 0, colorRampPalette(c("grey80", 'red3'))(100)) ## Grey--Red
mycol_BG <- clr_darken(shift = 0, colorRampPalette(c("grey80", 'dodgerblue3'))(100)) ## Grey--Blue
mycol_Vir <- clr_desaturate(shift = 0.2, viridis(100))
mycol_Plsm <- clr_desaturate(shift = 0.2,  plasma(100))

## Divergent gradient
mycol_Gfsh <- clr_desaturate(shift = 0.1, colorRampPalette(c('#69D2E7', '#A7DBD8', '#E0E4CC', '#F38630', '#FA6900'))(100))
mycol_Spec <- clr_darken(shift = 0, colorRampPalette(c(brewer.pal(11,"Spectral")))(100)[-c(48:53)])
mycol_RnBw <- clr_lighten(shift = 0.2, clr_desaturate(shift = 0.2, matlab.like(n = 100)))
mycol_RYB <- clr_darken(shift = 0, colorRampPalette(c(rev(brewer.pal(11,"RdYlBu"))))(20))

## Discrete:
mycol_7 <- clr_desaturate(shift = 0.2, c("#CF5556", "darkorange2", "#F9BF2D", "olivedrab3", "#0D8D36", "#47BAC1", "dodgerblue3"))
mycol_14 <- clr_desaturate(shift = 0, c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2')[c(1:3, 5:7)], 'grey20'))
mycol_60 <- clr_desaturate(shift = 0.1,
                           c(
                                   colorRampPalette(c("#2077b4", 'white'))(8)[1:6],
                                   colorRampPalette(c("#d62628", 'white'))(8)[1:6],
                                   colorRampPalette(c("#FF8C00", 'white'))(8)[1:6],
                                   colorRampPalette(c("#56a554", 'white'))(8)[1:6],
                                   colorRampPalette(c("#9466bd", 'white'))(8)[1:6],
                                   colorRampPalette(c("#e377c1", 'white'))(8)[1:6],
                                   colorRampPalette(c("#18bdcf", 'white'))(8)[1:6],
                                   colorRampPalette(c("#bcbd21", 'white'))(8)[1:6],
                                   colorRampPalette(c("#8b564c", 'white'))(8)[1:6],
                                   colorRampPalette(c("#A6D854FF", 'white'))(8)[1:6]
                           )
)
mycol_50 <- clr_desaturate(shift = 0.1,
                           c(
                                   colorRampPalette(c("#2077b4", 'white'))(7)[1:5],
                                   colorRampPalette(c("#d62628", 'white'))(7)[1:5],
                                   colorRampPalette(c("#FF8C00", 'white'))(7)[1:5],
                                   colorRampPalette(c("#56a554", 'white'))(7)[1:5],
                                   colorRampPalette(c("#9466bd", 'white'))(7)[1:5],
                                   colorRampPalette(c("#e377c1", 'white'))(7)[1:5],
                                   colorRampPalette(c("#18bdcf", 'white'))(7)[1:5],
                                   colorRampPalette(c("#bcbd21", 'white'))(7)[1:5],
                                   colorRampPalette(c("#8b564c", 'white'))(7)[1:5],
                                   colorRampPalette(c("#A6D854FF", 'white'))(7)[1:5]
                           )
)
mycol_40 <- clr_desaturate(shift = 0.1,
                           c(
                                   colorRampPalette(c("#2077b4", 'white'))(6)[1:4],
                                   colorRampPalette(c("#d62628", 'white'))(6)[1:4],
                                   colorRampPalette(c("#FF8C00", 'white'))(6)[1:4],
                                   colorRampPalette(c("#56a554", 'white'))(6)[1:4],
                                   colorRampPalette(c("#9466bd", 'white'))(6)[1:4],
                                   colorRampPalette(c("#e377c1", 'white'))(6)[1:4],
                                   colorRampPalette(c("#18bdcf", 'white'))(6)[1:4],
                                   colorRampPalette(c("#bcbd21", 'white'))(6)[1:4],
                                   colorRampPalette(c("#8b564c", 'white'))(6)[1:4],
                                   colorRampPalette(c("#A6D854FF", 'white'))(6)[1:4]
                           )
)
mycol_30 <- clr_desaturate(shift = 0.1,
                           c(
                                   colorRampPalette(c("#2077b4", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#d62628", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#FF8C00", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#56a554", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#9466bd", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#e377c1", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#18bdcf", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#bcbd21", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#8b564c", 'white'))(7)[c(1,3,5)],
                                   colorRampPalette(c("#A6D854FF", 'white'))(7)[c(1,3,5)]
                           )
)
mycol_20 <- mycol_40[seq(1, 40, 2)]
mycol_10 <- mycol_40[seq(2, 40, 4)]

theme_Publication <- function(base_size=14, base_family="Helvetica") {
        (theme_foundation(base_size=base_size, base_family=base_family) +
                 theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
                        text = element_text(),
                        panel.background = element_rect(colour = NA),
                        plot.background = element_rect(colour = NA),
                        panel.border = element_rect(colour = NA),
                        axis.title = element_text(face = "bold",size = rel(1)),
                        axis.title.y = element_text(angle=90,vjust =2),
                        axis.title.x = element_text(vjust = -0.2),
                        axis.text = element_text(),
                        axis.line = element_line(colour="black"),
                        axis.ticks = element_line(),
                        panel.grid.major = element_line(colour="#f0f0f0"),
                        panel.grid.minor = element_blank(),
                        legend.key = element_rect(colour = NA),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.key.size= unit(0.2, "cm"),
                        legend.margin = unit(0, "cm"),
                        legend.title = element_text(face="italic"),
                        plot.margin=unit(c(10,5,5,5),"mm"),
                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                        strip.text = element_text(face="bold")
                ))
}

scale_fill_Publication <- function(...){
        discrete_scale("fill","Publication",
                       manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c",
                                             "#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
        discrete_scale("colour","Publication",
                       manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c",
                                             "#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

####--------------------------------------------------------------------------------------------------------------------
Range01 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
}

O <- function(x, y, uniq = T){
        x <- as.vector(x)
        y <- as.vector(y)
        if(uniq){
                x <- unique(x)
                y <- unique(y)
        }
        out <- c(length(x), length(y), length(intersect(x, y)))
        names(out) <- c('x', 'y', 'common')
        return(out)
}

WinsorizeVector <- DescTools::Winsorize

L <- function(x) {
        length(x)
}

U <- function(x) {
        unique(x)
}

LU <- function(x) {
        length(unique(x))
}

H <- function(x, i = 5, j = 5) {
        x[1:i, 1:j]
}

S <- function(x) {
        cat(x, sep = '\n')
}
Table <- function(...) {table(..., useNA = 'ifany')}

ScaleColor <- function(color_vector = mycol_Spec, length) {
        new_color_vector <-
                color_vector[floor(Range01(1:length) * (length(color_vector) - 1)) + 1]
        return(new_color_vector)
}

ShowPalette <- function(color_vector) {
        col.mtx <- t(matrix(rep(1:length(
                color_vector
        ), 2), ncol = 2))
        pheatmap(
                col.mtx,
                color = color_vector,
                cluster_cols = F,
                legend = F,
                border_color = NA
        )
}

PlotPDF <- function(title, width, height, ...)
{
        pdf(
                paste0(Plot_dir, title, ".pdf"),
                width = width,
                height = height,
                useDingbats = F
        )
}

PlotPNG <- function(title, width, height, res = 300, ...)
{
        png(
                paste0(Plot_dir, title, ".png"),
                bg = "transparent",
                width = width * 300,
                height = height * 300,
                res = res
        )
}

GGPDF <- function(plot, title, width, height, ...)
{
        ggsave(
                filename = paste0(title, '.pdf'),
                plot = plot,
                device = 'pdf',
                path = Plot_dir,
                width = width,
                height = height,
                units = "in",
                dpi = 300,
                limitsize = F,
                ...
        )
}

GGPNG <- function(plot, title, width, height, ...)
{
        ggsave(
                filename = paste0(title, 'png'),
                plot = plot,
                device = 'png',
                path = Plot_dir,
                width = width,
                height = height,
                units = "in",
                dpi = 300,
                limitsize = F,
                ...
        )
}


WriteCSV <- function(x, title, ...)
{
        write_csv(x, file = paste0(Meta_dir, title, '.csv'), ...)
}

ReadCSV <- function(title, ...)
{
        read.csv(file = paste0(Meta_dir, title, '.csv'), ...)
}

# ConvertGeneSpecies <- function(genes, from, to) {
#         require("biomaRt")
#         human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#         mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#         if (from == 'mouse' & to == 'human') {
#                 genes.list = getLDS(
#                         attributes = c("mgi_symbol"),
#                         filters = "mgi_symbol",
#                         values = genes,
#                         mart = mouse,
#                         attributesL = c("hgnc_symbol"),
#                         martL = human,
#                         uniqueRows = T
#                 )
#         } else if (from == 'human' & to == 'mouse') {
#                 genes.list = getLDS(
#                         attributes = c("hgnc_symbol"),
#                         filters = "hgnc_symbol",
#                         values = genes ,
#                         mart = human,
#                         attributesL = c("mgi_symbol"),
#                         martL = mouse,
#                         uniqueRows = T
#                 )
#         } else {
#                 stop('Must be from human to mouse or from mouse to human')
#         }
#         output <- unique(genes.list[, 2])
#         return(output)
# }

## Source data from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
mouse_human_genes = read.csv("~/Documents/Bioinformatics/r/db/HOM_MouseHumanSequence.rpt", sep="\t")
ConvertGeneSpecies <- function(genes, from, to){
        output = c()
        if (from == 'mouse' & to == 'human') {
                for(gene in genes){
                        class_key =
                                (mouse_human_genes %>%
                                         dplyr::filter(Symbol == gene &
                                                        Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
                        if(!identical(class_key, integer(0)) ){
                                human_genes =
                                        (mouse_human_genes %>%
                                                 dplyr::filter(DB.Class.Key == class_key &
                                                                Common.Organism.Name=="human"))[,"Symbol"]
                                for(human_gene in human_genes){
                                        output = append(output, human_gene)
                                }
                        }
                }
        }
        else if (from == 'human' & to == 'mouse') {
                for(gene in genes){
                        class_key =
                                (mouse_human_genes %>%
                                         dplyr::filter(Symbol == gene &
                                                        Common.Organism.Name=="human"))[['DB.Class.Key']]
                        if(!identical(class_key, integer(0)) ){
                                mouse_genes =
                                        (mouse_human_genes %>%
                                                 dplyr::filter(DB.Class.Key == class_key &
                                                                Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
                                for(mouse_gene in mouse_genes){
                                        output = append(output, mouse_gene)
                                }
                        }
                }
        }
        return (output)
}

ConvertGeneID <- function(genes, species, from = 'id', to = 'symbol') {
        if (species == 'human') {
                library('EnsDb.Hsapiens.v86')
                if(from == 'id' & to == 'symbol'){
                        output <- ensembldb::select(EnsDb.Hsapiens.v86,
                                                    keys= genes,
                                                    keytype = "GENEID",
                                                    columns = c("SYMBOL","GENEID"))
                } else if (from == 'symbol' & to == 'id'){
                        output <- ensembldb::select(EnsDb.Hsapiens.v86,
                                                    keys= genes,
                                                    keytype = "SYMBOL",
                                                    columns = c("SYMBOL","GENEID"))
                } else {
                        stop('Wrong "from" or "to" parameters!')
                }
        } else if (species == 'mouse') {
                library('EnsDb.Mmusculus.v79')
                if(from == 'id' & to == 'symbol'){
                        output <- ensembldb::select(EnsDb.Mmusculus.v79,
                                                    keys= genes,
                                                    keytype = "GENEID",
                                                    columns = c("SYMBOL","GENEID"))
                } else if (from == 'symbol' & to == 'id'){
                        output <- ensembldb::select(EnsDb.Mmusculus.v79,
                                                    keys= genes,
                                                    keytype = "SYMBOL",
                                                    columns = c("SYMBOL","GENEID"))
                } else {
                        stop('Wrong "from" or "to" parameters!')
                }
        } else {
                stop('Species must be human or mouse')
        }
        return(output)
}

HyperGeometric <- function(overlap, x, y, total){
        out <- phyper(overlap, y, total-y, x, lower.tail= F)
        return(out)
}
Plot2WayVenn <- function(all_genes, x, y, x_name = 'x', y_name = 'y'){
        df <- data.frame(gene = U(all_genes),
                         X = F,
                         Y = F)
        for(i in 1:LU(all_genes)){
                if(df$gene[i] %in% x){df$X[i] <- T}
                if(df$gene[i] %in% y){df$Y[i] <- T}
        }
        pvl <- HyperGeometric(x = LU(x),
                              y = LU(y),
                              overlap = O(x, y)[3],
                              total = LU(all_genes))
        colnames(df)[2:3] <- c(x_name, y_name)
        plot(euler(df[, 2:3], shape = "ellipse"), quantities = T, main = paste0('p < ', pvl), col = mycol_10)
}

my.scale_colour_distiller <- function(..., type = "seq", palette = 1, direction = -1,
                                      values = NULL, space = "Lab", na.value = "grey50",
                                      guide = "colourbar", aesthetics = "colour", n_pal = 11) {
        # warn about using a qualitative brewer palette to generate the gradient
        type <- match.arg(type, c("seq", "div", "qual"))
        if (type == "qual") {
                warning("Using a discrete colour palette in a continuous scale.\n  Consider using type = \"seq\" or type = \"div\" instead", call. = FALSE)
        }
        continuous_scale(aesthetics, "distiller",
                         gradient_n_pal(brewer_pal(type, palette, direction)(n_pal), values, space), na.value = na.value, guide = guide, ...)
        # NB: 6 colours per palette gives nice gradients; more results in more saturated colours which do not look as good
}

GetProb <- function(x){
        z <- scale(x)
        p <- pnorm(z, lower.tail = T)
        return(p)
}

GetKeggGene <- function(kegg_id){
        x <- x <- KEGGREST::keggGet(kegg_id)[[1]]$GENE
        x <- x[seq(0, L(x), 2)]
        x <- gsub("\\;.*","", x)
        return(x)
}
