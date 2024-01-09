library(clusterProfiler)
library(DESeq2 )
library(enrichplot)
library(biomaRt)
library(pathview)

project <- 'PPMI'
trait <- 'CONCOHORT_DEFINITION'
TimeStep = 'V08'
modality <- 'mRNA'
  
setwd('./Year2/MOGDx/')
load(paste0('./data/',project , '/raw/Genetic_Subset/' , TimeStep,'/',modality,'_processed.RData'))

subtypes <- levels(as.factor(datMeta[[trait]]))
gene_list <- list()
gene_list[['Condition']] <- datMeta[colnames(datExpr) , 'CONCOHORT_DEFINITION']
top_genes = c()
for (subtype1 in subtypes[1:length(subtypes)-1]) {
  subtypes = subtypes[subtypes != subtype1]
  for (subtype2 in subtypes)  {
    if (subtype1 != subtype2) {
      res <- results(dds , contrast = c(trait , subtype1 , subtype2))
    }
    res <- res[!is.na(res$padj) , ]
    res <- res[((res$padj < 0.05) & (abs(res$log2FoldChange) > 0.1)) , ]
    print(dim(res))
    top_genes <-  res$log2FoldChange
    names(top_genes) <- substr(rownames(res) , 1,15)
    gene_list[[paste0(subtype1,'vs.',subtype2)]] <- top_genes

  }a
}

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "ensembl_gene_id", values = rownames(datExpr), mart = ensembl)
gene_info <-gene_info[ (gene_info$external_gene_name != '') , ]
rownames(gene_info) <- gene_info$ensembl_gene_id

datExpr_filt <- datExpr[gene_info$ensembl_gene_id , ]
rownames(datExpr_filt) <- gene_info$external_gene_name

write.csv(datExpr_filt , file = paste0('./data/PPMI/raw/Genetic_Subset/Feature_Importance/',modality,'/',TimeStep,'_',modality,'_datExpr.csv'))

for (i in 2:4) {
  common_genes <- intersect(gene_info$ensembl_gene_id , names(gene_list[[i]]))
  gene_list[[i]] = gene_list[[i]][common_genes]
  names(gene_list[[i]]) <- gene_info[common_genes , "external_gene_name"]
}

saveRDS(gene_list, file = paste0('./data/PPMI/raw/Genetic_Subset/Feature_Importance/',modality,'/',TimeStep,'_',modality,'_DEGs.rds'))
