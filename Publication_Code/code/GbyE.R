library(tidyverse); library(readxl) 

# Testing GbE

#Loading support packages
source("https://raw.githubusercontent.com/liu-xinrui/GbyE/main/GbyE.R")
#source("https://raw.githubusercontent.com/liu-xinrui/data/main/gapit_functions.txt")
source("~/Uni/Doctorate/gapit_functions_080425.txt")

# Individuals for experiment
sample_meta_drought <- read_xlsx("C:/Users/swirl/OneDrive/Documents/Uni/Doctorate/Samples/All_samples.xlsx") %>% filter(!is.na(NSWID) & `Genotyping Purpose`=="Drought_experiment") %>% dplyr::select(NSWID, Drought_Exp_ID)

# Data import
# GAPIT says: Missing data is not allowed for numerical genotype data

gt_datafile="~/Uni/Doctorate/Samples/Genotyping/Report-DMela25-10753/Report-DMela25-10753/raw/Report-DMela25-10753_SNP.csv"
dart_onerow <- read.csv(file = gt_datafile,  header = FALSE)

gt_indivs_col_keep <- which(unname(unlist(dart_onerow[7,-c(1:16)])) %in% sample_meta_drought$NSWID)
dart_onerow_rminf <- dart_onerow[,-c(1:16)]

## Numeric gt data
gt_only_num <- dart_onerow_rminf[-c(1:7),]
gt_only_num[gt_only_num == "-"] <- NA; gt_only_num[gt_only_num==1] <- 3; gt_only_num[gt_only_num==2] <- 1; gt_only_num[gt_only_num==3] <- 2
gt_only_num[] <- lapply(gt_only_num, function(x) as.numeric(x))

indivs <- unname(unlist(dart_onerow_rminf[7,]))

## SNP information
SNP_name <- dart_onerow[-c(1:7),1]; SNP_name <- gsub("-", "_", SNP_name)
SNP_CHR <- sub("^(MqA_CHR[0-9]{2}).*$", "\\1", dart_onerow[-c(1:7), 1])
SNP_pos <- as.numeric(sub(".*:(\\d+).*", "\\1",  dart_onerow[-c(1:7), 1])) + 150
SNP_data <- data.frame(SNP_name, as.numeric(sub("^MqA_CHR(\\d+)$", "\\1", SNP_CHR)), as.numeric(SNP_pos))

## Phenotype data
# 1) Having 30 v 60, Low v High
# [1] "The first is assumed to be additive effect and the second is assumed to be interactive effect"

Pheno_raw_1=read.csv("~/Uni/Doctorate/Ch2 Stressor/Publication/Data/MR_score_reinoc_comp_meta.csv")
Pheno_raw_2=read.csv("~/Uni/Doctorate/Ch2 Stressor/Publication/Data/MR_score_GT_pred.csv")

Pheno_raw_Drought_len <- left_join(Pheno_raw_2, Pheno_raw_1) %>% 
  dplyr::select(NSWID, Difference, Drought_len) %>% 
  filter(!is.na(NSWID), !is.na(Difference))

Pheno_raw_Drought_sev <- left_join(Pheno_raw_2, Pheno_raw_1) %>% 
  dplyr::select(NSWID, Difference, Drought_sev) %>% 
  filter(!is.na(NSWID), !is.na(Difference))

Pheno_raw_wider_Drought_len <- pivot_wider(Pheno_raw_Drought_len, id_cols = NSWID, values_from = Difference, names_from = Drought_len)
Pheno_raw_wider_Drought_sev <- pivot_wider(Pheno_raw_Drought_sev, id_cols = NSWID, values_from = Difference, names_from = Drought_sev)

colnames(Pheno_raw_wider_Drought_len)[1] <- "Taxa"; colnames(Pheno_raw_wider_Drought_sev)[1] <- "Taxa" 

# 2) Droughted v WR

Pheno_raw_1=read.csv("~/Uni/Doctorate/Ch2 Stressor/Publication/Data/MR_score_reinoc_comp_meta.csv")
Pheno_raw_2=read.csv("~/Uni/Doctorate/Ch2 Stressor/Publication/Data/MR_score_GT_pred.csv")

Pheno_raw_DRvWR <- left_join(Pheno_raw_2, Pheno_raw_1) %>% 
  dplyr::select(NSWID, COI_Scores_Inoc1, COI_Scores_Inoc2) %>% 
  filter(!is.na(NSWID), !is.na(COI_Scores_Inoc1))
colnames(Pheno_raw_DRvWR)[1] <- "Taxa"

int_taxa = intersect(Pheno_raw_wider_Drought_len$Taxa, indivs)
Pheno_raw_wider_Drought_len_only_intersect <- Pheno_raw_wider_Drought_len[Pheno_raw_wider_Drought_len$Taxa %in% int_taxa,]
##### Filtering and finalising datasets
#Remove genotypic data without phenotypes to maintain data consistency
indivs_int_GT_cols <- match(Pheno_raw_wider_Drought_len_only_intersect$Taxa, indivs)
colnames(gt_only_num) <- indivs; rownames(gt_only_num) <- SNP_name
gt_only_intersect_ord <- gt_only_num[, indivs_int_GT_cols]
gt_only_intersect_ord_t <- t(gt_only_intersect_ord)

## Filter missingness 
dim(gt_only_intersect_ord_t)
gt_only_intersect_ord_t_filtloc <- gt_only_intersect_ord_t[,(colSums(is.na(gt_only_intersect_ord_t))/nrow(gt_only_intersect_ord_t))<0.1] # Only keep loci with less than 90% missingness
gt_only_intersect_ord_t_filtind <- gt_only_intersect_ord_t_filtloc[(rowSums(is.na(gt_only_intersect_ord_t_filtloc))/ncol(gt_only_intersect_ord_t_filtloc))<0.1,] # Only keep indvs with less than 90% missingness
gt_only_intersect_ord_t <- gt_only_intersect_ord_t_filtind
dim(gt_only_intersect_ord_t)

## Impute remaining missingness (4699 out of 508354) - find 
# Find most common allele across each snp
comm_loci <- apply(gt_only_intersect_ord_t, 2, function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
})
# For each missing loci, impute with the most common loci
  # Find index of missingness and replace
gt_only_intersect_ord_t_imputed <- gt_only_intersect_ord_t

na_idx <- is.na(gt_only_intersect_ord_t_imputed)
gt_only_intersect_ord_t_imputed[na_idx] <- rep(comm_loci, each = nrow(gt_only_intersect_ord_t))[na_idx]
gt_only_intersect_ord_t <- gt_only_intersect_ord_t_imputed

##### Filter pheno and SNPs to matching
Pheno_raw_wider_Drought_len_only_intersect <- Pheno_raw_wider_Drought_len[Pheno_raw_wider_Drought_len$Taxa %in% rownames(gt_only_intersect_ord_t),]
Pheno_raw_wider_Drought_sev_only_intersect <- Pheno_raw_wider_Drought_sev[Pheno_raw_wider_Drought_sev$Taxa %in% rownames(gt_only_intersect_ord_t),]

SNP_data_intersect <- SNP_data[SNP_data$SNP_name %in% colnames(gt_only_intersect_ord_t),]

### Append extra rows onto genotype

row.names(gt_only_intersect_ord_t) <- Pheno_raw_wider_Drought_len_only_intersect$Taxa
colnames(gt_only_intersect_ord_t) <- SNP_data_intersect$SNP_name
#gt_only_intersect_ord_t_extrw <- rbind(colnames(gt_only_intersect_ord_t), gt_only_intersect_ord_t)
gt_only_intersect_ord_t_extrw <- rownames_to_column(as.data.frame(gt_only_intersect_ord_t), var = "taxa")

##### Drought severity as env

#Remove phenotypic data without genotypes to maintain data consistency
setwd("~/Uni/Doctorate/Ch2 Stressor/GBE")
#Run GbyE

## 1)
myGbyE_Drought_len=GbyE(GD=gt_only_intersect_ord_t_extrw,
                        GM=SNP_data_intersect,
                        Y=as.data.frame(Pheno_raw_wider_Drought_len_only_intersect),
                        PCA.total=3,
                        gwas=T,
                        gs=T,
                        plot=T,
                        gwas.model="MLM",
                        method="gapit")

myGbyE_Drought_sev=GbyE(GD=gt_only_intersect_ord_t_extrw,
                        GM=SNP_data_intersect,
                        Y=as.data.frame(Pheno_raw_wider_Drought_sev_only_intersect),
                        PCA.total=3,
                        gwas=T,
                        gs=T,
                        plot=T,
                        file.output = T,
                        gwas.model="MLM",
                        method="gapit")

## 2)
myGbyE_DrVWR=GbyE(GD=gt_only_intersect_ord_t_extrw,
                        GM=SNP_data_intersect,
                        Y=as.data.frame(Pheno_raw_DRvWR),
                        PCA.total=3,
                        gwas=T,
                        gs=T,
                        plot=T,
                        file.output = T,
                        gwas.model="MLM",
                        method="gapit")
tmp_GBE <- myGbyE_DrVWR$GWAS$GbyE  # Interactive effect
tmp_Orig <- myGbyE_DrVWR$GWAS$adde # Additive effect

GBE_DrVWR_Summ <- data.frame(cbind(SNP=tmp_GBE$SNP, Chromosome=tmp_GBE$Chromosome, SNP_POS=tmp_GBE$Position., Int_eff = tmp_GBE$effect, Add_eff = tmp_Orig$effect, IntMinAdd_eff = tmp_GBE$effect-tmp_Orig$effect))
write.csv(GBE_DrVWR_Summ, file="GBE_DrVWR_Summ.csv", row.names=T)
save.image(file="myGbyE_Drought_sev_len_res.Rd")


## Troubleshooting
# - gt in numeric (number gt), rows and rownames as taxa, cols and colnames as SNP_names
# - Error in xtfrm.data.frame(x) : cannot xtfrm data frames; conver Y as dataframe 
# - Error in GM.L2[, 2] + max(GM.L2[, 2]) ; Ensure chromosome and positions are numeric only
# - No '-' in SNPNames
