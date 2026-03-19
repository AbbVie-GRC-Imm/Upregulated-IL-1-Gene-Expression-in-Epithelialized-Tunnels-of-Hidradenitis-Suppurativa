# Load the Spatial Transcriptomics dataset
HS_RCTD <- readRDS("RCTD_all_53_combined.rds")

# Get unique sample names and add "X" prefix to each as per naming convention
imagenames <- unique(Hs_all$sample)
imagenames <- paste0("X", imagenames)

# Initialize empty data frames to store correlation results for each cytokine
cor_TNF <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL1B <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL1A <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL17A <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL17F <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL23A <- data.frame(matrix(NA, nrow = 53, ncol = 26))
cor_IL12B <- data.frame(matrix(NA, nrow = 53, ncol = 26))  # Note: cor_IL12B is included but not used further

# Assign row and column names for better interpretability
rownames(cor_TNF) = imagenames
rownames(cor_IL1B) = imagenames
rownames(cor_IL1A) = imagenames
rownames(cor_IL17A) = imagenames
rownames(cor_IL17F) = imagenames
rownames(cor_IL23A) = imagenames

# Set column names based on cell type columns (columns 7 to 32 in the metadata)
colnames(cor_TNF) = colnames(HS_RCTD@meta.data)[7:32]
colnames(cor_IL1B) = colnames(HS_RCTD@meta.data)[7:32]
colnames(cor_IL1A) = colnames(HS_RCTD@meta.data)[7:32]
colnames(cor_IL17A) = colnames(HS_RCTD@meta.data)[7:32]
colnames(cor_IL17F) = colnames(HS_RCTD@meta.data)[7:32]
colnames(cor_IL23A) = colnames(HS_RCTD@meta.data)[7:32]

# Loop through each sample to compute correlations
for(i in 1:53) {
  # Subset data for current sample (removes X prefix)
  data1 <- subset(x = HS_RCTD, subset = sample == gsub("X", "", imagenames[i]))
  
  # Extract relevant metadata and add cytokine expression as new columns
  meta <- data1@meta.data
  meta <- meta[,7:32]  # cell type abundance columns
  meta <- cbind(meta, FetchData(data1, vars = "TNF"))
  meta <- cbind(meta, FetchData(data1, vars = "IL1B"))
  meta <- cbind(meta, FetchData(data1, vars = "IL1A"))
  meta <- cbind(meta, FetchData(data1, vars = "IL17A"))
  meta <- cbind(meta, FetchData(data1, vars = "IL17F"))
  meta <- cbind(meta, FetchData(data1, vars = "IL23A"))
  
  # For each cell type (column), compute Pearson's correlation with cytokine expression and store
  for(j in 1:26){
    cor_TNF[i,j]   = cor(meta[,j], meta$TNF)
    cor_IL1B[i,j]  = cor(meta[,j], meta$IL1B)
    cor_IL1A[i,j]  = cor(meta[,j], meta$IL1A)
    cor_IL17A[i,j] = cor(meta[,j], meta$IL17A)
    cor_IL17F[i,j] = cor(meta[,j], meta$IL17F)
    cor_IL23A[i,j] = cor(meta[,j], meta$IL23A)
  }
}

# Save correlation matrices to CSV files for downstream analyses
write.table(cor_IL17A, "Correlation_between_celltypes_and_IL17A_53.csv", quote=F, sep=",")
write.table(cor_IL17F, "Correlation_between_celltypes_and_IL17F_53.csv", quote=F, sep=",")
write.table(cor_IL1A,  "Correlation_between_celltypes_and_IL1A_53.csv",  quote=F, sep=",")
write.table(cor_IL1B,  "Correlation_between_celltypes_and_IL1B_53.csv",  quote=F, sep=",")
write.table(cor_IL23A, "Correlation_between_celltypes_and_IL23A_53.csv", quote=F, sep=",")
write.table(cor_TNF,"Correlation_between_celltypes_and_TNF_53.csv",   quote=F, sep=",")
