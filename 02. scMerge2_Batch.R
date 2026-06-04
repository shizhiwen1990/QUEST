suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
})

sce <- as.SingleCellExperiment(seu_obj)
sce$batch <- sce$patient_id
log_expr_mat <- log1p(GetAssayData(seu_obj, assay = "RNA", slot = "counts"))
log_expr_mat <- as.matrix(log_expr_mat)
seg_custom <- intersect(c("MALAT1","B2M","RPL28","RPL10","RPS27","RPS12","TMSB4X","RPLP1" ), rownames(log_expr_mat))

scMerge2_res <- scMerge2(
  exprsMat = log_expr_mat,        
  batch = sce$batch,              
  ctl = seg_custom,               
  k_pseudoBulk = 50,              
  verbose = TRUE
)

expr_corrected <- scMerge2_res$newY
expr_corrected <- as.matrix(expr_corrected)
expr_corrected <- t(expr_corrected)

write.table(data.frame(ID = rownames(expr_corrected), expr_corrected), 
            file = "scMerge2.txt", 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

