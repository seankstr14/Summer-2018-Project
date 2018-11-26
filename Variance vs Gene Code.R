install.packages('C:/Users/5p336/Documents/R/win-library/Rttf2pt1_1.3.7.zip', repos = NULL)
library(ggplot2)
library(extrafont)


# same as generate_vector_from_search_terms with another name for ease of remembering
generate_vector_from_search_terms <- function(v_search_terms)  {
  return (paste("^", v_search_terms, "$", sep = "", collapse = "|"))
}

# read in data and mapping file
dat.CodeSet = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/Normalized_CodeString_data.txt", header = T, sep = "\t")
dat.map = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/map.txt", header = T, sep = "\t")
#dat.map$Microbiota = as.character(dat.map$Microbiota)
#dat.map$Microbiota[c(4:9, 16:21, 28:33, 43:48)] = "colonized"
#dat.map$Microbiota = factor(dat.map$Microbiota)
dat.genes = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/gene_table.txt", header = T, sep = "\t")

# ID each set of data points (e.g. GF SI) 
v_TxGrp = paste(dat.map$Microbiota, dat.map$Diet, dat.map$Tissue, sep = "__")
v_TxGrp_WO_TIssue = paste(dat.map$Microbiota, dat.map$Diet, sep = "__")
v_TxGrp.unique = unique(v_TxGrp_WO_TIssue)
#v_TxGrp.unique = v_TxGrp.unique[grep("BSH_low", v_TxGrp.unique,invert = TRUE)]

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

treatment_1 =  "BSH_high__Turmeric"
treatment_2 = "GF__Turmeric"

#Treatment 1 Sample IDs
v_indices.matches1_SI = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_1, "__SI", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)
v_indices.matches1_PC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_1, "__PC", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)
v_indices.matches1_DC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_1, "__DC", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)

#Treatment 2 Sample IDs
v_indices.matches2_SI = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_2, "__SI", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)
v_indices.matches2_PC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_2, "__PC", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)
v_indices.matches2_DC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(paste(treatment_2, "__DC", sep = ""), v_TxGrp)]), dat.CodeSet$Sample.ID)

#generate endogenous & Housekeeping gene mean + Standard Error + P-Value for treatment 1
matrix1SI = dat.CodeSet[v_indices.matches1_SI, -1]
matrix1PC = dat.CodeSet[v_indices.matches1_PC, -1]
matrix1DC = dat.CodeSet[v_indices.matches1_DC, -1]
se1SI = rep(1, length(v_indices.endogenousGenes))
logse1SI = se1SI
se1PC = rep(1, length(v_indices.endogenousGenes))
logse1PC = se1PC
se1DC = rep(1, length(v_indices.endogenousGenes))
logse1DC = se1DC
names(se1SI) = colnames(matrix1SI[,-c(69:75)])
names(logse1SI) = colnames(matrix1SI[,-c(69:75)])
names(se1PC) = colnames(matrix1PC[,-c(69:75)])
names(logse1PC) = colnames(matrix1PC[,-c(69:75)])
names(se1DC) = colnames(matrix1DC[,-c(69:75)])
names(logse1DC) = colnames(matrix1DC[,-c(69:75)])
for (i in 1:length(v_indices.endogenousGenes)) {
  se1SI[i] = sd(matrix1SI[ ,i]/sqrt(length(v_indices.matches1_SI)))
  logse1SI[i] = log2(se1SI[i])
  se1PC[i] = sd(matrix1PC[ ,i]/sqrt(length(v_indices.matches1_PC)))
  logse1PC[i] = log2(se1PC[i])
  se1DC[i] = sd(matrix1DC[ ,i]/sqrt(length(v_indices.matches1_DC)))
  logse1DC[i] = log2(se1DC[i])
}


#generate endogenous & housekeeping gene mean + Standard Error + P-Value for treatment 2
matrix2SI = dat.CodeSet[v_indices.matches2_SI, -1]
matrix2PC = dat.CodeSet[v_indices.matches2_PC, -1]
matrix2DC = dat.CodeSet[v_indices.matches2_DC, -1]
se2SI = rep(1, length(v_indices.endogenousGenes))
logse2SI = se2SI
se2PC = rep(1, length(v_indices.endogenousGenes))
logse2PC = se2PC
se2DC = rep(1, length(v_indices.endogenousGenes))
logse2DC = se2DC
for (i in 1:length(v_indices.endogenousGenes)) {
  se2SI[i] = sd(matrix2SI[ ,i]/sqrt(length(v_indices.matches2_SI)))
  logse2SI[i] = log2(se2SI[i])
  se2PC[i] = sd(matrix2PC[ ,i]/sqrt(length(v_indices.matches2_PC)))
  logse2PC[i] = log2(se2PC[i])
  se2DC[i] = sd(matrix2DC[ ,i]/sqrt(length(v_indices.matches2_DC)))
  logse2DC[i] = log2(se2DC[i])
}

#Plot standard error vs gene
seDataSI = as.data.frame(t(rbind(se1SI, se2SI)))
seDataPC = as.data.frame(t(rbind(se1PC, se2PC)))
seDataDC = as.data.frame(t(rbind(se1DC, se2DC)))
se_Data = as.data.frame(cbind(seDataSI, seDataPC, seDataDC))
seData = as.data.frame(se_Data[order(rowSums(se_Data), decreasing=TRUE),])
ggplot(seData, aes(x = reorder(row.names(seData), -seData[,3])), colour = variable, shape = variable) + 
  geom_point(aes(y = se1SI, col = treatment_1, shape = "SI")) + 
  geom_point(aes(y = se2SI, col = treatment_2, shape = "SI")) +
  geom_point(aes(y = se1PC, col = treatment_1, shape = "PC")) + 
  geom_point(aes(y = se2PC, col = treatment_2, shape = "PC")) +
  geom_point(aes(y = se1DC, col = treatment_1, shape = "DC")) + 
  geom_point(aes(y = se2DC, col = treatment_2, shape = "DC")) +
  theme(axis.text.x = element_text(angle = 60), legend.position="bottom") + 
  ggtitle(paste(treatment_1, "vs", treatment_2, sep = " ")) + xlab("gene") + ylab("standard error variance") +
  scale_colour_manual(name = "Treatment", values = c("blue", "red"))
ggsave(filename = paste("SE ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
       device = "pdf", width = 11, height = 8.5, units = "in")

#plot log2 standard error vs gene
log2seDataSI = as.data.frame(t(rbind(logse1SI, logse2SI)))
log2seDataPC = as.data.frame(t(rbind(logse1PC, logse2PC)))
log2seDataDC = as.data.frame(t(rbind(logse1DC, logse2DC)))
log2_seData = as.data.frame(cbind(log2seDataSI, log2seDataPC, log2seDataDC))
log2seData = as.data.frame(log2_seData[order(rowSums(se_Data), decreasing=TRUE),])
ggplot(log2seData, aes(x = reorder(row.names(log2seData), -log2seData[,3])), colour = variable, shape = variable) + 
  geom_point(aes(y = logse1SI, col = treatment_1, shape = "SI")) + 
  geom_point(aes(y = logse2SI, col = treatment_2, shape = "SI")) +
  geom_point(aes(y = logse1PC, col = treatment_1, shape = "PC")) + 
  geom_point(aes(y = logse2PC, col = treatment_2, shape = "PC")) +
  geom_point(aes(y = logse1DC, col = treatment_1, shape = "DC")) + 
  geom_point(aes(y = logse2DC, col = treatment_2, shape = "DC")) +
  theme(axis.text.x = element_text(angle = 60), legend.position="bottom") + 
  ggtitle(paste(treatment_1, " vs ", treatment_2, sep = "")) + xlab("gene") + ylab("Log2 (standard error variance)") +
  scale_colour_manual(name = "Treatment", values = c("blue", "red"))
ggsave(filename = paste("log2 SE", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), 
       device = "pdf", width = 11, height = 8.5, units = "in")
  }
}
