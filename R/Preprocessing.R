library(reshape2)
library(kableExtra)
library(plotly)
library(vsn)
library(tibble)
library(pheatmap)
library(SummarizedExperiment)
library(dplyr)
library(wateRmelon)
library(minfi)

source('~/MOGDx2.0/R/preprocess_functions.R')

setwd('~/MOGDx2.0/') 

dataset <- 'PPMI'
project <- 'All'
trait <- 'CONCOHORT_DEFINITION'
TimeStep = 'V08'

ALL <- c('Hyposmia' , '' , 'Sporadic' , 'RBD' , 'RBD/Hyposmia' , 'Genetic')
nonGenetic <- c('Hyposmia' , '' , 'Sporadic' , 'RBD' , 'RBD/Hyposmia')
Genetic <- c('' , 'Genetic')
  
# -------------------------------------------------------------------------
# META File Generation ----------------------------------------------------
# -------------------------------------------------------------------------

# The meta data is coded into the head folder of the dataset

meta <- read.csv(paste0('./data/',dataset,'/datMeta.csv'))
meta <- meta %>% filter(EVENT_ID == TimeStep) %>% filter(!((CONCOHORT_DEFINITION %in% c('Prodromal' , "Parkinson's Disease")) & (Subgroup == '')))
meta <- meta %>% filter(PATNO %in% common_idx$PATNO)
print(dim(meta))
write.csv(meta , file = paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv'))

# -------------------------------------------------------------------------
# mRNA pre-processing -----------------------------------------------------
# -------------------------------------------------------------------------

# Pull in Count Matrices --------------------------------------------------
data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/mRNA/mRNA_',TimeStep,'.csv') , row.names = 1 , check.names = FALSE)

count_mtx <- data
count_mtx <- count_mtx[, !(duplicated(colnames(count_mtx)))]

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

# Get intersection of count and meta
common_idx <- intersect(colnames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[ , common_idx]
datMeta <- datMeta[common_idx , ]

# Perform differential expression analysis
diff_expr_res <- diff_expr(count_mtx , datMeta , trait , 500 , 'mRNA')

# Save differential expression results
datExpr <- diff_expr_res$datExpr
datMeta <- diff_expr_res$datMeta
dds <- diff_expr_res$dds
top_genes <- diff_expr_res$top_genes
save(datExpr, datMeta, dds, top_genes, file=paste0('./data/',dataset,'/',TimeStep,'/mRNA/mRNA_processed.RData'))
rm(datExpr , datMeta , dds , top_genes , data , count_mtx)

# -------------------------------------------------------------------------
# miRNA preprocessing -----------------------------------------------------
# -------------------------------------------------------------------------
 
#Pull in count data
data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/miRNA/miRNA_',TimeStep,'.csv') , row.names = 1 , check.names = FALSE)

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

# Get Intersection of ID's
count_mtx <- floor(data)
count_mtx <- count_mtx[ , !(duplicated(colnames(count_mtx)))] 
common_idx <- intersect(colnames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[,common_idx]
datMeta <- datMeta[common_idx , ]

# Perform differential expression analysis
trait <- 'CONCOHORT_DEFINITION'
diff_expr_res <- diff_expr(count_mtx , datMeta , trait , 200 , 'mRNA')

# Save differential expression results
datExpr <- diff_expr_res$datExpr
datMeta <- diff_expr_res$datMeta
dds <- diff_expr_res$dds
top_genes <- diff_expr_res$top_genes
save(datExpr, datMeta, dds, top_genes, file=paste0(paste0('./data/',dataset,'/',TimeStep,'/miRNA/miRNA_processed.RData')))
rm(datExpr , datMeta , dds , top_genes , data , count_mtx)

# -------------------------------------------------------------------------
# DNAm preprocessing ------------------------------------------------------
# -------------------------------------------------------------------------

# Load CpG Counts ---------------------------------------------------------
data <- readRDS(paste0('./data/',dataset,'/',TimeStep,'/DNAm/Methyl_',TimeStep,'.danet.rds'))
count_mtx <- getBeta(data)

to_keep = complete.cases(count_mtx) #removed 0 cpg sites
length(to_keep) - sum(to_keep)

count_mtx <- t(count_mtx[to_keep,])
count_mtx <- count_mtx[data$Basename , ]
rownames(count_mtx) <- as.character(data$PATNO)

# Compute the variance across CpG sites
cpg_variances <- colVars(count_mtx)

# Sort the variances in descending order and get the indices
sorted_indices <- order(cpg_variances, decreasing = TRUE)

# Select the top 100000 most variable CpG sites
num_top_cpg <- 300000
top_cpg_indices <- sorted_indices[1:num_top_cpg]

count_mtx <- count_mtx[ , top_cpg_indices]

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

# Get Intersection of ID's
count_mtx <- count_mtx[!(duplicated(rownames(count_mtx))) , ] 
common_idx <- intersect(rownames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[common_idx , ]
datMeta <- datMeta[common_idx , ]

# Run glmnet to get CpG's associated with phenotypes of interest ------------
trait <- "CONCOHORT_DEFINITION"
phenotypes <- datMeta[,c('PATNO' , trait)]
colnames(phenotypes)

traits <- c(trait)

traitResults <- lapply(traits, function(trait) {
  cvTrait(count_mtx, phenotypes, traits, nFolds = 10)
})

cpg_sites <- c()
for (res in traitResults) { 
  trait_coefs <- coef(res$model , s = "lambda.min")
  cpg_sites_tmp <- c()
  for (coefs in trait_coefs) {
    class_coefs <- rownames(coefs)[which(coefs != 0)]
    if (length(class_coefs) > 1) {
      class_coefs <- class_coefs[2:length(class_coefs)]
      cpg_sites_tmp <- unique(c(cpg_sites_tmp ,class_coefs ))
    }
  }
  cpg_sites[[res$trait]] <- cpg_sites_tmp
}


# Save CpG sites and expression data
datExpr <- count_mtx

save(cpg_sites , datExpr , datMeta , file = paste0('./data/',dataset,'/',TimeStep,'/DNAm/DNAm_processed.RData'))
rm(datExpr , datMeta  , cpg_sites , data , count_mtx)

# -------------------------------------------------------------------------
# CSF Biomarkers ----------------------------------------------------------
# -------------------------------------------------------------------------
# Read in Count Matrix
data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/CSF/CSF_',TimeStep,'.csv') , row.names = 1 , check.names = FALSE)

count_mtx <- data
count_mtx <- count_mtx[,colSums(is.na(count_mtx))<0.5*nrow(count_mtx)] #remove columns with more than 50% NA

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

#Intersect to common IDs
count_mtx <- count_mtx[!(duplicated(rownames(count_mtx))) , ] 
common_idx <- intersect(rownames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[common_idx , ]
datMeta <- datMeta[common_idx , ]

# Replace NA's with 0's 
count_mtx[is.na(count_mtx)] <- 0 

# Perform log transform on count matrix to give normal distribution resemblance
#count_mtx_log <- log(count_mtx)

# Run glmnet R2 regression to get CNV's of interest
trait <- "CONCOHORT_DEFINITION"
phenotypes <- datMeta[,c('PATNO' , trait)]
colnames(phenotypes)
traits <- c(trait)

traitResults <- lapply(traits, function(trait) {
  cvTrait(as.matrix(count_mtx), phenotypes, traits, nFolds = 10)
})

csf_sites <- c()
for (res in traitResults) { 
  trait_coefs <- coef(res$model , s = "lambda.min")
  csf_sites_tmp <- c()
  for (coefs in trait_coefs) {
    class_coefs <- rownames(coefs)[which(coefs != 0)]
    if (length(class_coefs) > 1) {
      class_coefs <- class_coefs[2:length(class_coefs)]
      csf_sites_tmp <- unique(c(csf_sites_tmp ,class_coefs ))
    }
  }
  csf_sites[[res$trait]] <- csf_sites_tmp
}

# Code to handle the error here
csf_sites <- c()
csf_sites[[trait]] <- colnames(count_mtx)

# Code to execute whether there's an error or not
# Save CNV's and expression
datExpr <- count_mtx
save(csf_sites , datExpr , datMeta , file = paste0('./data/',dataset,'/',TimeStep,'/CSF/CSF_processed.RData'))
rm(datExpr , datMeta , csf_sites , data , count_mtx)

# ------------------------------------------------------------------------------------------
# MDS-UPDRS - Motor Symptoms ---------------------------------------------------------------
# ------------------------------------------------------------------------------------------

# Read in Count Matrix
data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/MDS-UPDRS/MDS-UPDRS_',TimeStep,'.csv') , row.names = 1 )

count_mtx <- data
count_mtx <- count_mtx[,colSums(is.na(count_mtx))<0.5*nrow(count_mtx)] #remove columns with more than 50% NA

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

#Intersect to common IDs
count_mtx <- count_mtx[!(duplicated(rownames(count_mtx))) , ] 
common_idx <- intersect(rownames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[common_idx , ]
datMeta <- datMeta[common_idx , ]

# Replace NA's with 0's 
count_mtx[is.na(count_mtx)] <- 0 

# Run glmnet R2 regression to get CNV's of interest
trait <- "CONCOHORT_DEFINITION"
phenotypes <- datMeta[,c('PATNO' , trait)]
colnames(phenotypes)
traits <- c(trait)

traitResults <- lapply(traits, function(trait) {
  cvTrait(as.matrix(count_mtx), phenotypes, traits, nFolds = 10)
})

p_sites <- c()
for (res in traitResults) { 
  trait_coefs <- coef(res$model , s = "lambda.min")
  p_sites_tmp <- c()
  for (coefs in trait_coefs) {
    class_coefs <- rownames(coefs)[which(coefs != 0)]
    if (length(class_coefs) > 1) {
      class_coefs <- class_coefs[2:length(class_coefs)]
      p_sites_tmp <- unique(c(p_sites_tmp ,class_coefs ))
    }
  }
  p_sites[[res$trait]] <- p_sites_tmp
}

# Save predictive features and expression
datExpr <- count_mtx
save(p_sites , datExpr , datMeta , file = paste0('./data/',dataset,'/',TimeStep,'/MDS-UPDRS/MDS-UPDRS_processed.RData'))
rm(datExpr , datMeta , p_sites , data , count_mtx)


# -------------------------------------------------------------------------
# Clinical Descriptor Features --------------------------------------------
# -------------------------------------------------------------------------
data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/Clinical/Clinical_',TimeStep,'.csv') , row.names = 1 )

count_mtx <- data
count_mtx <- count_mtx[,colSums(is.na(count_mtx))<0.5*nrow(count_mtx)] #remove columns with more than 50% NA

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

#Intersect to common IDs
count_mtx <- count_mtx[!(duplicated(rownames(count_mtx))) , ] 
common_idx <- intersect(rownames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[common_idx , ]
datMeta <- datMeta[common_idx , ]

# Save predictive features and expression
clin_sites <- c()
clin_sites[[trait]] <- colnames(count_mtx)
datExpr <- count_mtx
save(clin_sites , datExpr , datMeta , file = paste0('./data/',dataset,'/',TimeStep,'/Clinical/Clinical_processed.RData'))
rm(datExpr , datMeta , clin_sites , data , count_mtx)

# -------------------------------------------------------------------------
# SNP Features ------------------------------------------------------------
# -------------------------------------------------------------------------

data <- read.csv(paste0('./data/',dataset,'/',TimeStep,'/SNP/SNP_',TimeStep,'.csv') , row.names = 1 )

count_mtx <- data
count_mtx <- count_mtx[,colSums(is.na(count_mtx))<0.5*nrow(count_mtx)] #remove columns with more than 50% NA

# Pull in Meta File
datMeta <- as.data.frame(read.csv(paste0('./data/',dataset,'/',TimeStep,'/',TimeStep,'_datMeta.csv')))
datMeta <- datMeta %>% filter(EVENT_ID == TimeStep , CONCOHORT_DEFINITION %in% c('Healthy Control' , "Parkinson's Disease" , 'Prodromal'),Subgroup %in% Genetic)
datMeta$CONCOHORT_DEFINITION <- as.factor(datMeta$CONCOHORT_DEFINITION)
levels(datMeta$CONCOHORT_DEFINITION) <- c('HC' , 'PD' , 'PL')
datMeta <-  datMeta[!(duplicated(datMeta$PATNO)) , ]
rownames(datMeta) <- as.character(datMeta$PATNO)

#Intersect to common IDs
count_mtx <- count_mtx[!(duplicated(rownames(count_mtx))) , ] 
common_idx <- intersect(rownames(count_mtx) , rownames(datMeta))
count_mtx <- count_mtx[common_idx , ]
datMeta <- datMeta[common_idx , ]

# Save predictive features and expression
snp_sites <- c()
snp_sites[[trait]] <- colnames(count_mtx)
datExpr <- count_mtx
save(snp_sites , datExpr , datMeta , file = paste0('./data/',dataset,'/',TimeStep,'/SNP/SNP_processed.RData'))
rm(datExpr , datMeta , snp_sites , data , count_mtx)

