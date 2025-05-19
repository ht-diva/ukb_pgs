
# I am using this LDpred2 tutorial:
# https://choishingwan.github.io/PRS-Tutorial/ldpred/

#library(remotes)
#remotes::install_github("https://github.com/privefl/bigsnpr.git")

# 0. Prepare workspace
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#----------------#
# 3. Load and transform the summary statistic file

path_gwas_chr22 <- "/scratch/dariush.ghasemi/projects/nf-pipeline-regenie/gwas_sample_output/gwas_sample/results/gwas/target.gwas.regenie.gz"
path_geno_chr22 <- "/processing_data/shared_datasets/ukbiobank/GWAS_pipeline/step1_pruned_vars/by_chromosome/chr22.filtered.bed"
path_rds_chr22  <- "/processing_data/shared_datasets/ukbiobank/GWAS_pipeline/step1_pruned_vars/by_chromosome/chr22.filtered.rds"
path_geno_test <- "/scratch/dariush.ghasemi/projects/ukb_pgs/data/test_chr22.bed"
path_rds_test  <- "/scratch/dariush.ghasemi/projects/ukb_pgs/data/test_chr22.rds"

# Read in the summary statistic file
gwas_chr22 <- bigreadr::fread2(path_gwas_chr22)

# remove columns
sumstats <- gwas_chr22 |>
  #dplyr::filter(GENPOS > 17096864, GENPOS < 25426545) |>
  dplyr::select(- c(TEST, EXTRA, CHISQ)) |>
  dplyr::mutate(
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
    )
  
# LDpred 2 require the header to follow the exact naming
names(sumstats) <- c(
  "chr",
  "pos",
  "rsid",
  "a0",
  "a1",
  "a1_freq",
  "INFO",
  "n_eff",
  "beta",
  "beta_se",
  "p",
  "MAF"
  )

# Transform the OR into log(OR)
#sumstats$beta <- log(sumstats$OR)

# Filter out hapmap SNPs
#sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

#----------------#
# 3. Calculate the LD matrix

# Get maximum amount of cores
NCORES <- bigstatsr::nb_cores()

# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL

# We want to know the ordering of samples in the bed file 
fam.order <- NULL

# preprocess the bed file (only need to do once for each data set)
bigsnpr::snp_readBed(path_geno_chr22)
#bigsnpr::snp_readBed(path_geno_test)


# now attach the genotype object
obj.bigSNP <- bigsnpr::snp_attach(path_rds_chr22)

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")

# perform SNP matching
info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE) #, return_flip_and_rev = TRUE, join_by_pos = FALSE
#  mutate(freq = ifelse(`_REV_`, 1 - freq, freq), `_REV_` = NULL, `_FLIP_`= NULL)

# to subset the genotype using Plink
#write.table(info_snp$rsid, "snps_chr22.list", quote = F, row.names = F)

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes

# Rename the data structures
CHR <- map$chr
POS <- map$pos

# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- bigsnpr::snp_asGeneticPos(CHR, POS, dir = ".")

#----------------#
# calculate LD
for (chr in 22:22) {
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- bigsnpr::snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

# file size in GB
#file.size(corr$sbk) / 1024^3

# We assume the fam order is the same across different chromosomes
#fam.order <- data.table::as.data.table(obj.bigSNP$fam)

# Rename fam order
# data.table::setnames(
#   fam.order,
#   c("family.ID", "sample.ID"),
#   c("FID", "IID")
#   )

#----------------#
# 4. Perform LD score regression
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

ldsc <- bigsnpr::snp_ldsc(
  ld, 
  length(ld), 
  chi2 = (df_beta$beta / df_beta$beta_se)^2,
  sample_size = df_beta$n_eff, 
  blocks = NULL
  )

ldsc <- bigsnpr::snp_ldsc2(corr0, df_beta)
h2_est <- ldsc[["h2"]]


#----------------#
# 5. Calculate the null R2 (quantitative trait)

# Reformat the phenotype file such that y is of the same order as the sample ordering in the genotype file
#y <- pheno[fam.order, on = c("FID", "IID")]

# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)

# null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
#   paste0("Height~Sex+", .) %>%
#   as.formula %>%
#   lm(., data = y) %>%
#   summary
 
# null.r2 <- null.model$r.squared

# infinitesimal model
beta_inf <- bigsnpr::snp_ldpred2_inf(corr, df_beta, h2 = h2_est)


#----------------#
# 7. Obtain model PRS

# calculate PRS for all samples
ind.test <- 1:nrow(genotype)

pred_inf <- bigstatsr::big_prodVec(
  genotype,
  beta_inf,
  #ind.row = ind.test,
  ind.col = info_snp$`_NUM_ID_`
  )


