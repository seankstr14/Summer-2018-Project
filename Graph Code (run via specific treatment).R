# same as generate_vector_from_search_terms with another name for ease of remembering
generate_vector_from_search_terms <- function(v_search_terms)  {
  return (paste("^", v_search_terms, "$", sep = "", collapse = "|"))
}

# read in data and mapping file
dat.CodeSet = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/Normalized_CodeString_data.txt", header = T, sep = "\t")
dat.map = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/map.txt", header = T, sep = "\t")
dat.genes = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/gene_table.txt", header = T, sep = "\t")

# ID each set of data points (e.g. GF SI) 
v_TxGrp = paste(dat.map$Microbiota, dat.map$Diet, dat.map$Tissue, sep = "__")
v_TxGrp.unique = unique(v_TxGrp)
v_TxGrp.unique = v_TxGrp.unique[grep("BSH_low", v_TxGrp.unique,invert = TRUE)]

#Housekeeping Genes
v_indices.housekeepingGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Housekeeping", dat.genes$Classification)]), colnames(dat.CodeSet))

#Endogenous Genes
v_indices.endogenousGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Endogenous", dat.genes$Classification)]), colnames(dat.CodeSet))


    treatment_1 =  "BSH_combined__Turmeric__PC"
    treatment_2 = "BSH_combined__Non_turmeric__PC"
    
    #Treatment 1 Sample IDs
    v_indices.matches1 = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_1, v_TxGrp)]), dat.CodeSet$Sample.ID)
    
    #Treatment 2 Sample IDs
    v_indices.matches2 = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_2, v_TxGrp)]), dat.CodeSet$Sample.ID)

    # compare the two treatments
    v_pval1 = rep(1, length(v_indices.endogenousGenes))
    pval_housekeeping = rep(1, length(v_indices.housekeepingGenes))
    for (i in 1:length(v_indices.endogenousGenes))  {
      v_pval1[i] = t.test(dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_1, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i]],
                          dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_2, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i]])$p.value
    }
    for (i in 1:length(v_indices.housekeepingGenes))  {
      pval_housekeeping[i] = t.test(dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_1, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i+68]],
                                    dat.CodeSet[grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep(treatment_2, v_TxGrp)]), dat.CodeSet$Sample.ID), v_indices.genes[i+68]])$p.value
    }
    table(v_pval <= 0.05)  ## all FALSE
    colnames(dat.CodeSet)[v_indices.genes[v_pval <= 0.05]]
    
    #generate endogenous & Housekeeping gene mean + Standard Error + P-Value for treatment 1
    matrix1 = dat.CodeSet[v_indices.matches1, -1]
    geneMean1 = rep(1, length(v_indices.endogenousGenes))
    geneMean1_housekeeping = rep(1, length(v_indices.housekeepingGenes))
    names(geneMean1) = colnames(matrix1[,-c(69:75)])
    names(geneMean1_housekeeping) = colnames(matrix1[,-c(1:68)])
    se1 = rep(1, length(v_indices.endogenousGenes))
    se1_housekeeping = rep(1, length(v_indices.housekeepingGenes))
    names(se1) = colnames(matrix1[,-c(69:75)])
    names(se1_housekeeping) = colnames(matrix1[,-c(1:68)])
    for (i in 1:length(v_indices.endogenousGenes)) {
      geneMean1[i] = mean(matrix1[ ,i])
      se1[i] = sd(matrix1[ ,i]/sqrt(length(v_indices.matches1)))
    }
    for (i in 1:length(v_indices.housekeepingGenes)) {
      geneMean1_housekeeping[i] = mean(matrix1[ ,i+68])
      se1_housekeeping[i] = sd(matrix1[ ,i+68]/sqrt(length(v_indices.matches1)))
    }
    
    #generate endogenous & housekeeping gene mean + Standard Error + P-Value for treatment 2
    matrix2 = dat.CodeSet[v_indices.matches2, -1]
    geneMean2 = rep(1, length(v_indices.endogenousGenes))
    geneMean2_housekeeping = rep(a, length(v_indices.housekeepingGenes))
    se2 = rep(1, length(v_indices.endogenousGenes))
    se2_housekeeping = rep(1, length(v_indices.housekeepingGenes))
    for (i in 1:length(v_indices.endogenousGenes)) {
      geneMean2[i] = mean(matrix2[ ,i])
      se2[i] = sd(matrix2[ ,i]/sqrt(length(v_indices.matches2)))
    }
    for (i in 1:length(v_indices.housekeepingGenes)) {
      geneMean2_housekeeping[i] = mean(matrix2[ ,i+68])
      se2_housekeeping[i] = sd(matrix2[ ,i+68]/sqrt(length(v_indices.matches2)))
    }
    
    #order data
    meanData1 = rbind(geneMean1, geneMean2)
    meanData = as.data.frame(meanData1[,order(colSums(meanData1), decreasing=TRUE)])
    meanData1_housekeeping = rbind(geneMean1_housekeeping, geneMean2_housekeeping)
    meanData_housekeeping = as.data.frame(meanData1_housekeeping[,order(colSums(meanData1_housekeeping), decreasing=TRUE)])
    seData1 = rbind(se1, se2)
    seData = as.data.frame(seData1[,order(colSums(meanData1), decreasing=TRUE)])
    seData1_housekeeping = rbind(se1_housekeeping, se2_housekeeping)
    seData_housekeeping = as.data.frame(seData1_housekeeping[,order(colSums(meanData1_housekeeping), decreasing=TRUE)])
    v_pval = v_pval1[order(colSums(meanData1), decreasing=TRUE)]
    v_pval_housekeeping = pval_housekeeping[order(colSums(meanData1_housekeeping), decreasing=TRUE)]
    
    #generate plot of endogenous genes w/ < 100 counts
    c = 0
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) <= 100) {
        c = c + 1
      }
    }
    graph_meanData = as.data.frame(matrix(nrow = 2, ncol =  c))
    graph_seData = matrix(nrow = 2, ncol = c)
    d = 1
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) <= 100) {
        graph_meanData[,d] = meanData[,i]
        nameofGene = colnames(meanData[i])
        names(graph_meanData)[d] = nameofGene
        graph_seData[,d] = seData[,i]
        d = d + 1
      }
    }
    pdf(paste("endogenous(100) ", treatment_1," vs ", treatment_2, ".pdf", sep = ""), width = 11, height = 8.5)
    bp = barplot(as.matrix(graph_meanData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Genes", ylab = "Counts",
                 ylim= c(0, max(graph_meanData) + 20), col=c("blue","green"), beside=TRUE, names.arg = names(graph_meanData), las=2)
    arrows(x0 = bp, y0 = as.matrix(graph_meanData) - graph_seData, y1= as.matrix(graph_meanData) + graph_seData, code=3, angle=90, length = 0.03)
    legend("topright", legend = c(treatment_1, treatment_2), fill = c("blue", "green"), cex = 0.75)
    dev.off()
    
    #generate plot of endogenous genes w/ 100 < x <= 1000 counts
    c = 0
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) >100 & mean(meanData[[i]]) <= 1000) {
        c = c + 1
      }
    }
    graph_meanData = as.data.frame(matrix(nrow = 2, ncol =  c))
    graph_seData = matrix(nrow = 2, ncol = c)
    d = 1
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) >100 & mean(meanData[[i]]) <= 1000) {
        graph_meanData[,d] = meanData[,i]
        nameofGene = colnames(meanData[i])
        names(graph_meanData)[d] = nameofGene
        graph_seData[,d] = seData[,i]
        d = d + 1
      }
    } 
    pdf(paste("endogenous(100-1000) ", treatment_1, " vs " , treatment_2, ".pdf", sep =""), width = 11, height = 8.5)
    bp = barplot(as.matrix(graph_meanData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Genes", ylab = "Counts",
                 ylim= c(0, max(graph_meanData) + 200), col=c("blue","green"), beside=TRUE, names.arg = names(graph_meanData), las = 2)
    arrows(x0 = bp, y0 = as.matrix(graph_meanData) - graph_seData, y1= as.matrix(graph_meanData) + graph_seData, code=3, angle=90, length = 0.03)
    legend("topright", legend = c(treatment_1, treatment_2), fill = c("blue", "green"), cex = 0.75)
    dev.off()
    
    #generate plot of endogenous genes > 1000
    c = 0
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) > 1000) {
        c = c + 1
      }
    }
    graph_meanData = as.data.frame(matrix(nrow = 2, ncol =  c))
    graph_seData = matrix(nrow = 2, ncol = c)
    d = 1
    for(i in 1:length(v_indices.endogenousGenes)){
      if(mean(meanData[[i]]) > 1000) {
        graph_meanData[,d] = meanData[,i]
        nameofGene = colnames(meanData[i])
        names(graph_meanData)[d] = nameofGene
        graph_seData[,d] = seData[,i]
        d = d + 1
      }
    }
    pdf(paste("endogenous (1000) ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), width = 11, height = 8.5)
    bp = barplot(as.matrix(graph_meanData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Genes", ylab = "Counts",
                 ylim= c(0, max(graph_meanData) + 1500), col=c("blue","green"), beside=TRUE, names.arg = names(graph_meanData))
    arrows(x0 = bp, y0 = as.matrix(graph_meanData) - graph_seData, y1= as.matrix(graph_meanData) + graph_seData, code=3, angle=90, length = 0.1)
    legend("topright", legend = c(treatment_1, treatment_2), fill = c("blue", "green"), cex = 0.75)
    dev.off()
    
    
    #generate significant difference plot
    a = 0
    for(i in 1:length(v_indices.endogenousGenes)){
      if(v_pval[i] <= 0.05) {
        a = a + 1
      }
    }
    sigData = as.data.frame(matrix(nrow = 2, ncol = a))
    se_sigData = matrix(nrow = 2, ncol = a)
    d = 1
    for(i in 1:length(v_indices.endogenousGenes)){
      if(v_pval[i] <= 0.05){
        sigData[,d] = meanData[,i]
        nameofGene = colnames(meanData[i])
        names(sigData)[d] = nameofGene
        se_sigData[,d] = seData[,i]
        d = d + 1
      }
    }
    pdf(paste("significant ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), width = 11, height = 8.5)
    bp = barplot(as.matrix(sigData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Significant Genes", ylab = "Counts", 
                 ylim= c(0, max(sigData) + 200), col=c("cyan","red"), beside=TRUE, names.arg = names(sigData), las=2)
    arrows(x0 = bp, y0 = as.matrix(sigData) - se_sigData, y1= as.matrix(sigData) + se_sigData, code=3, angle=90, length = 0.1)
    legend("topleft", legend = c(treatment_1, treatment_2), fill = c("cyan", "red"), cex = 0.75)
    dev.off()
    
    #Plot housekeeping genes w/ <= 1500 counts
    c = 0
    for(i in 1:length(v_indices.housekeepingGenes)){
      if(mean(meanData_housekeeping[[i]]) <= 1500) {
        c = c + 1
      }
    }
    graph_meanData = as.data.frame(matrix(nrow = 2, ncol =  c))
    graph_seData = matrix(nrow = 2, ncol = c)
    d = 1
    for(i in 1:length(v_indices.housekeepingGenes)){
      if(mean(meanData_housekeeping[[i]]) <= 1500) {
        graph_meanData[,d] = meanData_housekeeping[,i]
        nameofGene = colnames(meanData_housekeeping[i])
        names(graph_meanData)[d] = nameofGene
        graph_seData[,d] = seData_housekeeping[,i]
        d = d + 1
      }
    }
    pdf(paste("housekeeping (low) ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), width = 11, height = 8.5)
    bp = barplot(as.matrix(graph_meanData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Housekeeping Genes", ylab = "Counts", 
                 ylim= c(0, max(graph_meanData) + 200), col=c("blue","gold"), beside=TRUE, names.arg = names(graph_meanData))
    arrows(x0 = bp, y0 = as.matrix(graph_meanData) - graph_seData, y1= as.matrix(graph_meanData) + graph_seData, code=3, angle=90, length = 0.1)
    legend("topright", legend = c(treatment_1, treatment_2), fill = c("blue", "gold"), cex = 0.75)
    dev.off()
    
    #Plot housekeeping genes w/ > 1500 counts
    c = 0
    for(i in 1:length(v_indices.housekeepingGenes)){
      if(mean(meanData_housekeeping[[i]]) > 1500) {
        c = c + 1
      }
    }
    graph_meanData = as.data.frame(matrix(nrow = 2, ncol =  c))
    graph_seData = matrix(nrow = 2, ncol = c)
    d = 1
    for(i in 1:length(v_indices.housekeepingGenes)){
      if(mean(meanData_housekeeping[[i]]) > 1500) {
        graph_meanData[,d] = meanData_housekeeping[,i]
        nameofGene = colnames(meanData_housekeeping[i])
        names(graph_meanData)[d] = nameofGene
        graph_seData[,d] = seData_housekeeping[,i]
        d = d + 1
      }
    }
    pdf(paste("housekeeping (high) ", treatment_1, " vs ", treatment_2, ".pdf", sep = ""), width = 11, height = 8.5)
    bp = barplot(as.matrix(graph_meanData), main = paste(treatment_1, treatment_2, sep = " vs. "), xlab = "Housekeeping Genes", ylab = "Counts",
                 ylim= c(0, max(graph_meanData) + 10000), col=c("blue","gold"), beside=TRUE, names.arg = names(graph_meanData))
    arrows(x0 = bp, y0 = as.matrix(graph_meanData) - graph_seData, y1= as.matrix(graph_meanData) + graph_seData, code=3, angle=90, length = 0.1)
    legend("topright", legend = c(treatment_1, treatment_2), fill = c("blue", "gold"), cex = 0.75)
    dev.off()
