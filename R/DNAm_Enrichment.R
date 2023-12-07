library(mCSEA)

dataset <- 'TCGA'
project <- 'KIPAN'
trait <- 'subtype'
modality <- 'DNAm'

setwd('~/MOGDx2.0/')
load(paste0('/data/',dataset , '/raw/', project ,'/',modality,'_processed.RData'))
## -----------------------------------------------------------------------------
betaTest <- t(datExpr)
subtypes <- levels(as.factor(datMeta[[trait]]))
for (subtype in subtypes) {
  print(subtype)
  phenoTest <- as.matrix(datMeta[[trait]])
  colnames(phenoTest) <- trait
  rownames(phenoTest) <- colnames(betaTest)
  phenoTest[ which(phenoTest[ , trait] != subtype) , ] <- 'Other'
  #phenoTest[,trait] <- as.factor(phenoTest[,trait])
  phenoTest <- as.data.frame(phenoTest)
  
  myRank <- rankProbes(betaTest, phenoTest, refGroup = 'Other' , caseGroup = subtype)
  
  myRank[!(names(myRank) %in% cpg_sites[[trait]])] <- 0.000001
  
  myResults <- mCSEATest(myRank, betaTest, phenoTest, 
                         regionsTypes = "genes", platform = "EPIC")
  
  save(myResults , file = paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/',modality,'/',subtype,'_mCSEA.RData'))
}

myResults_all <- list()
for (subtype in subtypes) {
  load(paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/',modality,'/',subtype,'_mCSEA.RData'))
  myResults_all[[subtype]] <- myResults
} 

top_genes <- c()
all_genes <- c()
out_table <- c()
for (result in myResults_all) {
  
  all_genes <- c(all_genes , rownames(result$genes))
  top_genes <- c(top_genes , rownames(result$genes)[which(result$genes$padj < 0.05)])
  result_table <- cbind(result$genes[order(result$genes$padj)[which(result$genes$padj < 0.05)] , c(1,2,5,6)] , 'Subtype' = rep(unique(result$pheno[result$pheno != 'Other']) , sum(result$genes$padj < 0.05)))
  out_table <- rbind(out_table , result_table )
  
}

write.table(top_genes[!duplicated(top_genes)] , file =  paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/DNAm_enrichment.txt') , 
            , quote = FALSE , row.names = FALSE , col.names = FALSE)
write.table(all_genes[!duplicated(all_genes)] , file =  paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/DNAm_background_enrichment.txt') , 
            , quote = FALSE , row.names = FALSE , col.names = FALSE)

#write.csv(out_table[order(out_table$NES , decreasing = TRUE) ,] , file = paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/DNAm/DNAm_enrichment.csv'))

#write.table(data.frame(out_table[order(out_table$NES , decreasing = TRUE) , 'NES'] ) , 
#            file = paste0('/data/',dataset , '/raw/',project,'/Feature_Importance/',modality,'/',modality,'_DEGs.rnk'),
#            sep = '\t' ,row.names = rownames(out_table) , quote = FALSE , col.names = FALSE)          
