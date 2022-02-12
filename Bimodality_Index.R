# This code can be used to identify bimodally expressed genes.
# The bimodality index is calculated using six methods:
# 1) K-means (KM)
# 2) Mixture modeling (MM)
# 3) Mixture modeling with k-means (MK_0)
# 4) Controlled mixture modeling with k-means (CMK_0) 
# 5) Mixture modeling with k-means excluding clusters <10% (MK_10)
# 6) Controlled mixture modeling with k-means excluding clusters <10% (CMK_10)
# The example given here shows how to download breast cancer 
# RNA-seq data from GDC, log2 transform the gene expression, 
# and calculate the bimodality index for each gene.
# If the user desires to use their own data, they may skip to the 
# section entitled 'Data Prep'. Cancer and control gene expression 
# in RPKM should be uploaded in 2 separate files and formatted such that 
# rows = genes, columns = samples, and first column = gene names.



#################### Data Download ####################

library(TCGAbiolinks)

# Download breast cancer gene expression data
query.mrna_cancer <- GDCquery(project = "TCGA-BRCA", 
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - FPKM",
                              sample.type = c("Primary Tumor"))
GDCdownload(query.mrna_cancer)
mrna_cancer <- GDCprepare(query = query.mrna_cancer,
                          save = TRUE, summarizedExperiment = FALSE,
                          save.filename = "TCGA-BRCA_gene.rda") 


# Download breast control tissue gene expression
query.mrna_control <- GDCquery(project = "TCGA-BRCA", 
                               experimental.strategy = "RNA-Seq",
                               data.category = "Transcriptome Profiling", 
                               data.type = "Gene Expression Quantification",
                               workflow.type = "HTSeq - FPKM",
                               sample.type = c("Solid Tissue Normal"))
GDCdownload(query.mrna_control)
mrna_control <- GDCprepare(query = query.mrna_control,
                           save = TRUE, summarizedExperiment = FALSE,
                           save.filename = "Control_TCGA-BRCA_gene.rda")


#################### Data Prep ####################


################ If loading data ###########################
## Cancer and control expression data should be uploaded 
## in 2 separate files and formatted such that rows = genes, 
## columns = samples, and first column = gene names.
## Use next two lines to read in files.
############################################################


# mrna_cancer = read.csv('myfile_cancer.csv') # cancer expression
# mrna_control = read.csv('myfile_control.csv') # control expression

# Log transform data
integer_cancer = apply(mrna_cancer[ , -1], 2, as.integer)
log_transformed_cancer = log2(integer_cancer + 1)
gene_name = as.list(as.matrix((mrna_cancer[ , 1])))

integer_control = apply(mrna_control[ , -1], 2, as.integer)
log_transformed_control = log2(integer_control + 1)

x = as.matrix(log_transformed_cancer)
y = as.matrix(log_transformed_control)

# Remove genes with very low expression
filtered = apply(x, 1, max)
cancer = x[which(filtered > 1) ,]
control = y[which(filtered > 1) ,]
saved_genes = gene_name[which(filtered > 1)]


#################### Bimodality Index ####################

library(mclust)

# Create output matrix
output_bimodal = matrix(nrow = nrow(cancer), ncol = 6)
rownames(output_bimodal) = saved_genes
colnames(output_bimodal) = c("KM", "MM", "MK_0", "CMK_0", "MK_10", "CMK_10")
set.seed(1)

for (n in 1:nrow(cancer)) {
  print(n)
  
  ## KM and MM
  
  km = kmeans(cancer[n , ], 2) # K-means for cancer
  pi1 = km$size[1] / length(km$cluster)
  pi2 = km$size[2] / length(km$cluster)
  u1 = km$centers[1]
  u2 = km$centers[2]
  sig1 = var(cancer[which(km$cluster == 1)])
  sig2 = var(cancer[which(km$cluster == 2)])
  KM = (sqrt(pi1 * pi2) * abs(u1 - u2)
        / sqrt(pi2 * sig1 + pi1 * sig2))
  output_bimodal[n, 1] = KM
  
  BIC_cancer <- mclustBIC(cancer[n , ], # Mixture model for cancer
                          G = 1:2)
  mod_cancer <- Mclust(cancer[n , ],
                       x = BIC_cancer,
                       G = 1:2)
  
  if (mod_cancer$G == 1) {
    MM = 0
    MK_0 = 0
    CMK_0 = 0
    MK_10 = 0
    CMK_10 = 0
    output_bimodal[n, 2] = MM # MM = 0 if unimodal
    output_bimodal[n, 3] = MK_0 # MK = 0 if unimodal
    output_bimodal[n, 4] = CMK_0 # CMK = 0 if unimodal
    output_bimodal[n, 5] = MK_10 # MK = 0 if unimodal
    output_bimodal[n, 6] = CMK_10 # CMK = 0 if unimodal
    
  } else {
    pi_MM_1 = length(which(mod_cancer$classification == 1)) / mod_cancer$n
    pi_MM_2 = length(which(mod_cancer$classification == 2)) / mod_cancer$n
    u_MM_1 = mod_cancer$parameters$mean[1]
    u_MM_2 = mod_cancer$parameters$mean[2]
    sig_MM_1 = var(cancer[which(mod_cancer$classification == 1)])
    sig_MM_2 = var(cancer[which(mod_cancer$classification == 2)])
    MM = (
      sqrt(pi_MM_1 * pi_MM_2) * abs(u_MM_1 - u_MM_2)
      / sqrt(pi_MM_2 * sig_MM_1 + pi_MM_1 * sig_MM_2)
    )
    output_bimodal[n, 2] = MM
    MK_0 = KM
    output_bimodal[n, 3] = MK_0
    
    ## MK_10
    
    if (pi1 < 0.10) {
      MK_10 = 0
      output_bimodal[n, 5] = MK_10 # MK_10 = 0 if cluster < 10%
      
    } else {
      if (pi2 < 0.10) {
        MK_10 = 0
        output_bimodal[n, 5] = MK_10 # MK_10 = 0 if cluster < 10%
        
      } else {
        MK_10 = KM
        output_bimodal[n, 5] = MK_10
      }
    }
    
    ## CMK_0 and CMK_10
    
    if (max(control[n ,]) < 1) {
      CMK_0 = MK_0
      CMK_10 = MK_10
      output_bimodal[n, 4] = CMK_0 # CMK = MK if control is very low expression
      output_bimodal[n, 6] = CMK_10 # CMK = MK if control is very low expression
      
    } else {
      BIC_control <- mclustBIC(control[n ,], # Mixture model for control
                               G = 1:2)
      mod_control <- Mclust(control[n ,],
                            x = BIC_control,
                            G = 1:2)
      
      if (mod_control$G == 1) {
        CMK_0 = MK_0
        CMK_10 = MK_10
        output_bimodal[n, 4] = CMK_0 # CMK = MK if control is unimodal
        output_bimodal[n, 6] = CMK_10 # CMK = MK if control is unimodal
        
      } else {
        km.control = kmeans(control[n ,], 2) # k-means for control
        pi1.control = km.control$size[1] / length(km.control$cluster)
        pi2.control = km.control$size[2] / length(km.control$cluster)
        u1.control = km.control$centers[1]
        u2.control = km.control$centers[2]
        sig1.control = var(control[which(km.control$cluster == 1)])
        sig2.control = var(control[which(km.control$cluster == 2)])
        
        u.control = c(u1.control, u2.control)
        pi.control = c(pi1.control, pi2.control)
        sig.control = c(sig1.control, sig2.control)
        u = c(u1, u2)
        pi = c(pi1, pi2)
        sig = c(sig1, sig2)
        
        smaller = which(u == min(u))
        larger = which(u == max(u))
        smaller.c = which(u.control == min(u.control))
        larger.c = which(u.control == max(u.control))
        
        CMK_pi = (MK_0 - (1 / (
          abs(pi.control[smaller.c] - pi[smaller])
          + abs(pi.control[larger.c] - pi[larger])
        )))
        CMK_u = (MK_0 - (1 / (
          abs(u.control[smaller.c] - u[smaller])
          + abs(u.control[larger.c] - u[larger])
        )))
        CMK_0 = max(c(CMK_pi, CMK_u))
        output_bimodal[n, 4] = CMK_0
        
        if (pi1.control < 0.10) {
          CMK_10 = MK_10
          output_bimodal[n, 6] = CMK_10 # CMK = MK if control cluster <10%
          
        } else {
          if (pi2.control < 0.10) {
            CMK_10 = MK_10
            output_bimodal[n, 6] = CMK_10 # CMK = MK if control cluster <10%
            
          } else {
            if (pi1 < 0.10) {
              CMK_10 = 0
              output_bimodal[n, 6] = CMK_10 # CMK_10 = 0 if cancer cluster < 10%
              
            } else {
              if (pi2 < 0.10) {
                CMK_10 = 0
                output_bimodal[n, 6] = CMK_10 # CMK_10 = 0 if cancer cluster < 10%
                
              } else {
                CMK_10 = CMK_0
                output_bimodal[n, 6] = CMK_10
              }
            }
          }
        }
      }
    }
  }
}

# Print a table of bimodality indices
output_table = matrix(nrow = 5, ncol = 6)
rownames(output_table) = c("BI > 0", "BI > 1", "BI > 1.4", "BI > 1.5", "BI > 2")
colnames(output_table) = c("KM", "MM", "MK_0", "CMK_0", "MK_10", "CMK_10")
cutoffs = c(0, 1, 1.4, 1.5, 2)

for (i in 1:5){
  for (j in 1:6){
    output_table[i,j] = length(which(output_bimodal[, j] > cutoffs[i]))
  }
}

output_table

# Save matrix of bimodality indices for each gene using each method
write.csv(output_bimodal, file = "breast_cancer_genes_bimodality.csv")


# Plot histograms of top 3 genes using each method 
ordered_KM = order(output_bimodal[, 1], na.last = TRUE, decreasing = TRUE)
ordered_MM = order(output_bimodal[, 2], na.last = TRUE, decreasing = TRUE)
ordered_MK = order(output_bimodal[, 5], na.last = TRUE, decreasing = TRUE)
ordered_CMK = order(output_bimodal[, 6], na.last = TRUE, decreasing = TRUE)

par(mfrow = c(3,4)) # Each column is different method
for (i in 1:3){
  hist(c(cancer[ordered_KM[i],]), # KM
       prob = FALSE, xlim = c(-2,5), 
       yaxt="n", xaxt="n",
       main = "", xlab = "log2(RPKM)", ylab = "", 
       breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_KM[i],]), 
               from = -2, to = 5, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
  
  hist(c(cancer[ordered_MM[i],]), # MM
       prob = FALSE, xlim = c(-2,5), 
       yaxt="n", xaxt="n",
       main = "", xlab = "log2(RPKM)", ylab = "", 
       breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_MM[i],]), 
               from = -2, to = 5, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
  
  hist(c(cancer[ordered_MK[i],]), # MK
       prob = FALSE, xlim = c(-2,15), 
       yaxt="n", xaxt="n",
       main = "", xlab = "log2(RPKM)", ylab = "", breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_MK[i],]), 
               from = -2, to = 15, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
  
  hist(c(cancer[ordered_CMK[i],]), # CMK
       prob = FALSE, xlim = c(-2,15), 
       yaxt="n", xaxt="n",
       main = "", xlab = "log2(RPKM)", ylab = "", 
       breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_CMK[i],]), 
               from = -2, to = 15, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
}

