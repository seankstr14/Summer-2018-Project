library(ggplot2)
library(reshape2)
library(gridExtra)

# same as generate_vector_from_search_terms with another name for ease of remembering
generate_vector_from_search_terms <- function(v_search_terms)  {
  return (paste("^", v_search_terms, "$", sep = "", collapse = "|"))
}

# read in data and mapping file
dat.CodeSet = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/Normalized_CodeString_data.txt", header = T, sep = "\t")
dat.map = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/map.txt", header = T, sep = "\t")
at.map$gender <- "male"
dat.map$gender[c(25:33, 37:39, 43:48)] = "female"
dat.genes = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/gene_table.txt", header = T, sep = "\t")

# ID each set of data points (e.g. GF SI) 
v_TxGrp = paste(dat.map$Microbiota, dat.map$Diet, dat.map$Tissue, sep = "__")
v_TxGrp.unique = unique(v_TxGrp)
v_TxGrp.unique = v_TxGrp.unique[grep("BSH_low", v_TxGrp.unique,invert = TRUE)]

#Housekeeping Genes
v_indices.housekeepingGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Housekeeping", dat.genes$Classification)]), colnames(dat.CodeSet))

#Endogenous Genes
v_indices.endogenousGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Endogenous", dat.genes$Classification)]), colnames(dat.CodeSet))

for (i_groupA in 1:(length(v_TxGrp.unique)-1))  {
  for (i_groupB in (i_groupA+1):length(v_TxGrp.unique))  {
    treatment_1 = v_TxGrp.unique[i_groupA]
    treatment_2 = v_TxGrp.unique[i_groupB]
    print(treatment_1)
    print(treatment_2)
    
    treatment_1 = "GF__Turmeric__SI"
    treatment_2 = "BSH_combined__Turmeric__SI"
    
#Treatment 1 Sample IDs
v_indices.matches1 = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_1, v_TxGrp)]), dat.CodeSet$Sample.ID)

#Treatment 2 Sample IDs
v_indices.matches2 = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_2, v_TxGrp)]), dat.CodeSet$Sample.ID)

# compare the two treatments
v_pval1 = rep(1, length(v_indices.endogenousGenes))
log_pval = rep(1, length(v_pval1))
for (i in 1:length(v_indices.endogenousGenes))  {
  v_pval1[i] = t.test(dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_1, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i]],
                      dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_2, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i]])$p.value
  log_pval[i] = -log10(v_pval1[i])
}
#table(v_pval1 <= 0.05)  ## all FALSE
#colnames(dat.CodeSet)[v_indices.genes[v_pval1 <= 0.05]]

#generate endogenous & Housekeeping gene mean + Standard Error + P-Value for treatment 1 & 2
matrix1 = dat.CodeSet[v_indices.matches1, -1]
geneMean1 = rep(1, length(v_indices.endogenousGenes))
names(geneMean1) = colnames(matrix1[,-c(69:75)])
matrix2 = dat.CodeSet[v_indices.matches2, -1]
geneMean2 = rep(1, length(v_indices.endogenousGenes))
geneMean2_housekeeping = rep(a, length(v_indices.housekeepingGenes))
for (i in 1:length(v_indices.endogenousGenes)) {
  geneMean1[i] = mean(matrix1[ ,i])
  geneMean2[i] = mean(matrix2[ ,i])
}
names(log_pval) = names(geneMean1)

#Generate fold change
foldChange = rep(1, length(v_indices.endogenousGenes))
absfoldChange = foldChange
for(i in 1:length(meanData)){
  foldChange[i] = geneMean1[i]/geneMean2[i]
  absfoldChange[i] = abs(foldChange[i])
}
names(foldChange) = names(geneMean1)
names(absfoldChange) = names(foldChange)
both = FALSE
for(i in 1:length(v_indices.endogenousGenes)) {
  if(log_pval[i] > -log10(0.05) & (foldChange[i] > 1.5 | foldChange[i] < 2/3)) {
    both = TRUE
  } else if (both == TRUE) {
    break
  }
}
checkPval = FALSE
for(i in 1:length(v_indices.endogenousGenes)){
  if(log_pval[i] > -log10(0.05) & both == FALSE) {
    checkPval = TRUE
  } else if (checkPval == TRUE) {
    break
  }
}
checkFoldChange = FALSE
for(i in 1:length(v_indices.endogenousGenes)){
  if((foldChange[i] > 1.5 | foldChange[i] < 2/3)  & both == FALSE ) {
    checkFoldChange = TRUE
  } else if (checkFoldChange == TRUE) {
    break
  }
}

#volcanoPlot Scratch
pvalThreshold = -log10(0.05)
volcanoPlotData = as.data.frame(t(rbind(foldChange, log_pval)))
if(both == TRUE) {
ggplot(volcanoPlotData, aes(x = foldChange, y = log_pval), colour = variable) + 
  geom_point() +
  geom_point(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), aes(col = "fold change past threshold")) +
  geom_text(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), 
            aes(label = row.names(subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3)), col = "fold change past threshold")) +
  geom_point(data = subset(volcanoPlotData, log_pval > pvalThreshold), aes(col = "pval < 0.05")) +
  geom_text(data = subset(volcanoPlotData, log_pval > pvalThreshold), 
            aes(label = row.names(subset(volcanoPlotData, log_pval > pvalThreshold)), col = "pval < 0.05")) +
  geom_point(data = subset(volcanoPlotData, log_pval > pvalThreshold & (foldChange > 1.5 | foldChange < 2/3)), aes(col = "both")) +
  geom_text(data = subset(volcanoPlotData, log_pval > pvalThreshold & (foldChange > 1.5 | foldChange < 2/3)), 
            aes(label = row.names(subset(volcanoPlotData, log_pval > pvalThreshold & (foldChange > 1.5 | foldChange < 2/3))), 
            col = "both")) + xlab("fold change") + ylab("-log(pval)") +
    geom_text(data = volcanoPlotData[30,], aes(label = "Glp2r", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[68,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[68,], aes(label = "Vip", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[12,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[12,], aes(label = "DCX", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[2,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[2,], aes(label = "Asbt", col = "Important Genes")) +
  geom_vline(xintercept = c(2/3, 1.5)) + geom_hline(yintercept = (pvalThreshold)) +
  ggtitle(paste(treatment_1, " compared to ", treatment_2, sep = "")) + theme(legend.position="bottom") + 
  scale_colour_manual(name = "Legend", values = c("blue", "green", "purple", "red"))
ggsave(filename = paste("volcano plot ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
       device = "pdf", width = 11, height = 8.5, units = "in")
} else if(checkPval == TRUE & checkFoldChange == TRUE) {
  ggplot(volcanoPlotData, aes(x = foldChange, y = log_pval), colour = variable) + 
    geom_point() +
    geom_point(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), aes(col = "fold change past threshold")) +
    geom_text(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), 
              aes(label = row.names(subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3)), col = "fold change past threshold")) +
    geom_point(data = subset(volcanoPlotData, log_pval > pvalThreshold), aes(col = "pval < 0.05")) +
    geom_text(data = subset(volcanoPlotData, log_pval > pvalThreshold), 
              aes(label = row.names(subset(volcanoPlotData, log_pval > pvalThreshold)), col = "pval < 0.05")) + 
    xlab("fold change") + ylab("-log(pval)") + 
    geom_point(data = volcanoPlotData[30,], aes(col = "Important Genes")) + 
    geom_text(data = volcanoPlotData[30,], aes(label = "Glp2r", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[68,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[68,], aes(label = "Vip", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[12,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[12,], aes(label = "DCX", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[2,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[2,], aes(label = "Asbt", col = "Important Genes")) +
    geom_vline(xintercept = c(2/3, 1.5)) + geom_hline(yintercept = (pvalThreshold)) +
    ggtitle(paste(treatment_1, " compared to ", treatment_2, sep = "")) + theme(legend.position="bottom") +
    scale_colour_manual(name = "Legend", values = c("green", "purple", "red"))
  ggsave(filename = paste("volcano plot ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
         device = "pdf", width = 11, height = 8.5, units = "in")
} else if(checkPval == TRUE & checkFoldChange == FALSE) {
  ggplot(volcanoPlotData, aes(x = foldChange, y = log_pval), colour = variable) + 
    geom_point() +
       geom_point(data = subset(volcanoPlotData, log_pval > pvalThreshold), aes(col = "pval")) +
    geom_text(data = subset(volcanoPlotData, log_pval > pvalThreshold), 
              aes(label = row.names(subset(volcanoPlotData, log_pval > pvalThreshold)), col = "pval")) + 
    xlab("fold change") + ylab("-log(pval)") +
    geom_point(data = volcanoPlotData[30,], aes(col = "Important Genes")) + 
    geom_text(data = volcanoPlotData[30,], aes(label = "Glp2r", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[68,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[68,], aes(label = "Vip", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[12,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[12,], aes(label = "DCX", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[2,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[2,], aes(label = "Asbt", col = "Important Genes")) +
    geom_vline(xintercept = c(2/3, 1.5)) + geom_hline(yintercept = (pvalThreshold)) +
    ggtitle(paste(treatment_1, " compared to ", treatment_2, sep = "")) + theme(legend.position="bottom") + 
    scale_colour_manual(name = "genes past threshold", values = c("red", "purple"))
  ggsave(filename = paste("volcano plot ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
         device = "pdf", width = 11, height = 8.5, units = "in")
} else if (checkFoldChange == TRUE & checkPval == FALSE){
  ggplot(volcanoPlotData, aes(x = foldChange, y = log_pval), colour = variable) + 
    geom_point() +
    geom_point(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), aes(col = "fold change")) +
    geom_text(data = subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3), 
              aes(label = row.names(subset(volcanoPlotData, foldChange > 1.5 | foldChange < 2/3)), col = "fold change")) +
    xlab("fold change") + ylab("-log(pval)") +
    geom_point(data = volcanoPlotData[30,], aes(col = "Important Genes")) + 
    geom_text(data = volcanoPlotData[30,], aes(label = "Glp2r", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[68,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[68,], aes(label = "Vip", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[12,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[12,], aes(label = "DCX", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[2,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[2,], aes(label = "Asbt", col = "Important Genes")) +
    geom_vline(xintercept = c(2/3, 1.5)) + geom_hline(yintercept = (pvalThreshold)) +
    ggtitle(paste(treatment_1, " compared to ", treatment_2, sep = "")) + theme(legend.position="bottom") + 
    scale_colour_manual(name = "genes past threshold", values = c("green", "purple"))
  ggsave(filename = paste("volcano plot ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
         device = "pdf", width = 11, height = 8.5, units = "in")
} else {
  ggplot(volcanoPlotData, aes(x = foldChange, y = log_pval), colour = variable) + 
    geom_point() + xlab("fold change") + ylab("-log(pval)") +
    geom_point(data = volcanoPlotData[30,], aes(col = "Important Genes")) + 
    geom_text(data = volcanoPlotData[30,], aes(label = "Glp2r", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[68,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[68,], aes(label = "Vip", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[12,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[12,], aes(label = "DCX", col = "Important Genes")) +
    geom_point(data = volcanoPlotData[2,], aes(col = "Important Genes")) +
    geom_text(data = volcanoPlotData[2,], aes(label = "Asbt", col = "Important Genes")) +
    geom_vline(xintercept = c(2/3, 1.5)) + geom_hline(yintercept = (pvalThreshold)) +
    ggtitle(paste(treatment_1, " compared to ", treatment_2, sep = "")) + theme(legend.position="bottom") + 
    scale_colour_manual(name = "genes past threshold", values = c("green", "purple"))
  ggsave(filename = paste("volcano plot ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
         device = "pdf", width = 11, height = 8.5, units = "in")
}


#Absolute value fold change histogram plot scratch
pdf(paste("histogram ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""))
hist(foldChange, main= paste(treatment_1, "vs", treatment_2, sep = " "))
dev.off()
  }
}



