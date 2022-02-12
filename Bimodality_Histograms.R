# This code was used to create histograms for top bimodally expressed breast cancer genes.
# 
# Data includes:
# 1) mrna_cancer.csv: log2 transformed cancer expression
# 2) mrna_control.csv: log2 transformed control expression
# 3) bimodal_index.csv: maxtrix of bimodality indices using MM and CM


library(ggpubr)
library(cowplot)

# Load data
load(file = "histograms.RData")

# Plot histograms for top 3 genes using each method 
ordered_MM = order(BI[, 1], na.last = TRUE, decreasing = TRUE)
ordered_CM = order(BI[, 2], na.last = TRUE, decreasing = TRUE)
gene_names = rownames(BI)

# Top genes (MM)
par(mfrow = c(1,3))
for (i in 1:3){
  hist(c(cancer[ordered_MM[i],]), # MM
       prob = FALSE, xlim = c(-2,15), 
       yaxt="n", xaxt="n",
       main = gene_names[i], xlab = "log2(RPKM)", ylab = "", 
       breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_MM[i],]), 
               from = -2, to = 15, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
}

# Top genes (CM)
par(mfrow = c(1,3)) # Each column is different method
for (i in 1:3){
  hist(c(cancer[ordered_CM[i],]), # CM
       prob = FALSE, xlim = c(-2,15), 
       yaxt="n", xaxt="n",
       main = gene_names[i], xlab = "log2(RPKM)", ylab = "", 
       breaks = 20)
  axis(side = 4)
  par(new=TRUE)
  plot(density(c(cancer[ordered_CM[i],]), 
               from = -2, to = 15, bw = .5),
       lwd = 1, col = "red", 
       main = "", xlab = "", ylab = "")
}


# Plot distributions for cancer and control for GSTM1
z = 18642
Group = c(rep("Cancer", 1102), rep("Control", 113))
gene = c(cancer[z ,], control[z ,])
data = data.frame(cbind(Group, gene))
data$gene = as.numeric(as.character(data$gene))

ggdensity(data, x = "gene", y = "..density..",
          rug = TRUE, fill = "Group", xlab = "log2(RPKM)",
          palette = c("red3", "deepskyblue")) + 
  xlim(-4, 10) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))

