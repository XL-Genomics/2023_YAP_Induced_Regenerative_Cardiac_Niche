db_dir <- '~/Documents/Bioinformatics/r/db/'


####--------------------------------------------------------------------------------------------------------------------
####  General genesets  ####
####--------------------------------------------------------------------------------------------------------------------
mhc_genes <- list(MHC1 = as.vector(read.csv(paste0(db_dir, 'mouse_mhc_gene.csv'), header = T)[, 1]),
                  MHC2 = as.vector(read.csv(paste0(db_dir, 'mouse_mhc_gene.csv'), header = T)[1:7, 2]))

chrX_genes <- readRDS(paste0(db_dir, 'mouse_chrX_genes.rds'))
chrY_genes <- readRDS(paste0(db_dir, 'mouse_chrY_genes.rds'))


####--------------------------------------------------------------------------------------------------------------------
####  Markersets from publication  ####
####--------------------------------------------------------------------------------------------------------------------
Yap_target <- list(
        Yap_target_cm_dev = as.vector(
                read.table(paste0(db_dir, 'yap_target_cm_atac_2019_DevCell.csv'), header = T)[, 1]),
        Yap_target_cm = as.vector(
                read.table(paste0(db_dir, 'yap_target_cm_scenic_Yuka_Unpub.txt'), header = F)[, 1]),
        Yap_target_cf = as.vector(
                read.table(paste0(db_dir, 'yap_target_cf_scenic_2019_GenesDev.txt'), header = T)[, 1]),
        Yap_target_yap5sa = as.vector(
                readRDS(paste0(db_dir, 'mouse_yap_target_yap5sa.rds'))),
        Yap_target_yap5sa_union = as.vector(
                readRDS(paste0(db_dir, 'mouse_yap_target_yap5sa_union.rds'))),
        Yap_target_cf_cutrun = readRDS(paste0(paste0(db_dir,
                                              'cf_yap_target.latscko_yap1cutrun_based_on_genedev_paper.rds')))[[1]]
)
Yap_target$Yap_target_cm[Yap_target$Yap_target_cm == 'Ptrf'] <- 'Cavin1' ## update symbol

#### Bmp4 induced genes XXXX ####
Bmp4_target <- list(
        Bmp4_target = read.csv(paste0(db_dir, 'Bmp4_induced_gene_list.csv'))$Ensembl_Gene_Symbol
)

#### Neonatal MI induced regenerative CM markers, 2020 DevCell EOlson ####
Neo_MI_CM_marker <- list(
        CM1 = c('Gm30382', 'Fhl2', 'Fgf13', 'Mhrt', 'Ank2'),
        CM2 = c('Casc5', 'Lockd', 'Top2a', 'Gm14091'),
        CM3 = c('Mid1', 'Ddc', '1700042O10Rik'),
        CM4 = c('Mb', 'Sod2', 'Myl3', 'Mdh2', 'Atp5b'),
        CM5 = c('Enah', 'Ankrd1', 'Cd44', 'Gm13601', 'Xirp2')
)

#### Atrial CF LatsCKO CF markers ####
LatsCKO_aCF_marker <- list(
        aCF1 = c('Rora', 'Mast4', 'Gsn', 'Cfh', 'Dcn', 'Abca8a', 'Mgp', 'Fmo2', 'Apoe', 'Zbtb20'),
        aCF2 = c('Igf1r', 'Fmn1', 'Large1', 'Tead1', 'Phldb2', 'Rock2', 'Scel', 'Ext1', 'Ptchd4', 'Clca3a1')
)

#### Ventricular YAP5SA CM C3ar1+Mac, C3+CF markers ####
Yap5sa_CM_Mac_CF_marker <- list(
        CM1 = c("Atp5j","Cox6b1","Atp5g1","Chchd10","Atp5c1","Cox7b","Ndufb10","Mpc2","Sdhb","Uqcrb","Ndufv2",
                "Etfa","Acadm","Ckmt2","Mgst3","Cox7a1","Mrpl42","Eci1","Phyh","Fabp3-ps1"),
        CM2 = c("Flnc","Sorbs2","Itga7","Ahnak","Tgm2","Prnp","Tns1","Ddb1","Eef2","Rtn4","Serpinh1","Ccn2",
                "Jam2","Rras2","Cavin1","Lmcd1","Rock2","Acta2","Parm1"),
        C3ar1_Mac = c('C3ar1', 'Cx3cr1', 'C5ar1', 'Rbpj', 'F13a1', 'Mrc1'),
        C3_FB = c('C3', 'Fstl1', 'Serpina3n', 'Cxcl14', 'Serping1')
)

#### Cardaic fibroblast subtype markers ####
CF_marker <- list(
        Resident_CF = c('Scara5', 'Gfpt2', 'Sema3c', 'Uap1', 'Adamts5'),
        Matrifibrocytes = c('Chad', 'Chrdl1', 'Ecrg4', 'Cilp', 'Cilp2', 'Comp')
)

#### Purkinje cell signatures ####
Purkinje_marker <- readRDS(paste0(db_dir, 'mouse_purkinje_genes.rds'))

#### Mono or bi-nucleated CM signatures ####
Mono_binucleated_adult_cm_marker <- readRDS(paste0(db_dir, 'mono_bi_nucleus_adult.list.rds'))
Mono_binucleated_p7_cm_marker <- readRDS(paste0(db_dir, 'mono_bi_nucleus_p7.list.rds'))


#### Post MI Zone signatures -- Calcagno et al Nat Card Res 2022 ####
# data <- as.data.frame(read_excel('external/44161_2022_160_MOESM4_ESM.xlsx')) ## Calcagno et al
# data$cluster...7 <- revalue(as.character(data$cluster...7),
#                             replace = c('1' = 'RZ', '2' = 'BZ1', '3' = 'BZ2', '4' = 'IZ'))
# PostMI_zone_marker <- split(data$gene, data$cluster...7)
# PostMI_zone_marker <- PostMI_zone_marker[ c('RZ', 'BZ1', 'BZ2', 'IZ') ]
# saveRDS(PostMI_zone_marker, paste0(db_dir, 'adult_mouse_post_mi_zone_marker.list.rds'))
PostMI_zone_marker <- readRDS(paste0(db_dir, 'adult_mouse_post_mi_zone_marker.list.rds'))
