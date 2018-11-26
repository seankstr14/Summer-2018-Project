
# same as generate_vector_from_search_terms with another name for ease of remembering
generate_vector_from_search_terms <- function(v_search_terms)  {
  return (paste("^", v_search_terms, "$", sep = "", collapse = "|"))
}

# read in data and mapping file
dat.CodeSet = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/Normalized_CodeString_data.txt", header = T, sep = "\t")
dat.map = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/map.txt", header = T, sep = "\t")
#dat.map$gender <- "male"
#dat.map$gender[c(25:33, 37:39, 43:48)] = "female"
dat.genes = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/gene_table.txt", header = T, sep = "\t")

# ID each set of data points (e.g. GF SI) 
v_TxGrp = paste(dat.map$Microbiota, dat.map$Diet, dat.map$Tissue, sep = "__")
v_TxGrp.unique = unique(v_TxGrp)
v_TxGrp.unique = v_TxGrp.unique[grep("BSH_low", v_TxGrp.unique,invert = TRUE)]

#Housekeeping Genes
v_indices.housekeepingGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Housekeeping", dat.genes$Classification)]), colnames(dat.CodeSet))

#Endogenous Genes
v_indices.endogenousGenes = grep(generate_vector_from_search_terms(dat.genes$Gene[grep("Endogenous", dat.genes$Classification)]), colnames(dat.CodeSet))

v_Tissue = dat.map$Tissue
v_indices.SI = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep("SI", v_Tissue)]), dat.CodeSet$Sample.ID)
v_indices.PC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep("PC", v_Tissue)]), dat.CodeSet$Sample.ID)
v_indices.DC = grep(generate_vector_from_search_terms(dat.map$Sample.ID[grep("DC", v_Tissue)]), dat.CodeSet$Sample.ID)

#generate correlation between Glp2r and other genes
for(i in 60:76) {
  if(i == 31) {
    next
  }
  glp2r.SI = rep(1,length(v_indices.SI))
  glp2r.PC = rep(1,length(v_indices.SI))
  glp2r.DC = rep(1,length(v_indices.SI))
  gene.SI = rep(1,length(v_indices.SI))
  gene.PC = rep(1,length(v_indices.SI))
  gene.DC = rep(1,length(v_indices.SI))
  for(j in 1:16){
    glp2r.SI[j] = dat.CodeSet$Glp2r[3*j - 2]
    glp2r.PC[j] = dat.CodeSet$Glp2r[3*j - 1]
    glp2r.DC[j] = dat.CodeSet$Glp2r[3*j]
    gene.SI[j] = dat.CodeSet[3*j - 2, i]
    gene.PC[j] = dat.CodeSet[3*j - 1, i]
    gene.DC[j] = dat.CodeSet[3*j, i]
  }
  
  #Linear Regression Analysis
  x = dat.CodeSet$Glp2r
  y = as.matrix(dat.CodeSet[i])
  ylab = colnames(dat.CodeSet[i])
  colonglp2r = c(glp2r.PC, glp2r.DC) 
  colongene = c(gene.PC, gene.DC)
  gene.lm = lm(y ~ x)
  linearRegression = summary(gene.lm)
  SIgene.lm = lm(gene.SI ~ glp2r.SI)
  SIlinearRegression = summary(SIgene.lm)
  PCgene.lm = lm(gene.PC ~ glp2r.PC)
  PClinearRegression = summary(PCgene.lm)
  DCgene.lm = lm(gene.DC ~ glp2r.DC)
  DClinearRegression = summary(DCgene.lm)
  colon.lm = lm(colongene ~colonglp2r)
  colonlinearRegression = summary(colon.lm)
  f <- linearRegression$fstatistic
  a <- SIlinearRegression$fstatistic
  b <- PClinearRegression$fstatistic
  c <- DClinearRegression$fstatistic
  d <- colonlinearRegression$fstatistic
  
  #Plot all Tissues in one plot 
  ymax = range(c(gene.SI, gene.PC, gene.DC))
  halfrange = (ymax[1] + ymax[2])/2
  dev.new(width = 600, height = 330, unit = "px")
  par(xpd = F, mar = par()$mar + c(0,0,0,7))
  plot(glp2r.SI, gene.SI, xlim = range(c(glp2r.SI, glp2r.PC, glp2r.DC)), ylim = ymax, 
       title(paste("Glp2r vs. ", ylab, sep = "")), col = "red", xlab = "Glp2r", ylab = ylab, pch = 15)
  abline(SIgene.lm, col = "red")
  points(glp2r.PC, gene.PC, col = "gold", pch = 17)
  abline(PCgene.lm, col = "gold")
  points(glp2r.DC, gene.DC, col = "blue", pch = 16)
  abline(DCgene.lm, col = "blue")
  abline(gene.lm)
  abline(colon.lm, col = "green")
  par(xpd=T)
  legend(225, ymax[1],legend = c("All", "SI", "PC", "DC", "Colon"), col = c("black", "red", "gold", "blue", "green"),
         title = "Key",pch = c(18, 15, 17, 16, 18), cex = 0.8, ncol = 2 )
  text(250, ymax[2], paste(" All slope = ",round(coef(linearRegression)[2], digits = 5), "\n All R^2 = ",
                       round(linearRegression$r.squared, digits = 5),",\n adj All R^2 = ", 
                       round(linearRegression$adj.r.squared, digits = 5), "\n All p-val = ", 
                       pf(f[1], f[2], f[3], lower=FALSE), sep = ""), cex = 2/3)
  text(250, (3/2) * halfrange, paste(" SI slope = ",round(coef(SIlinearRegression)[2], digits = 5), "\n SI R^2 = ",
              round(SIlinearRegression$r.squared, digits = 5),",\n adj SI R^2 = ", 
              round(SIlinearRegression$adj.r.squared, digits = 5), "\n SI p-val = ", 
              pf(a[1], a[2], a[3], lower=FALSE), sep = ""), cex = 2/3, col = "red")
  text(250, (7/6) * halfrange, paste(" PC slope = ",round(coef(PClinearRegression)[2], digits = 5), "\n PC R^2 = ",
                       round(PClinearRegression$r.squared, digits = 5),",\n adj PC R^2 = ", 
                       round(PClinearRegression$adj.r.squared, digits = 5), "\n PC p-val = ", 
                       pf(b[1], b[2], b[3], lower=FALSE), sep = ""), cex = 2/3, col = "gold")
  text(250, (6/7) * halfrange, paste(" DC slope = ",round(coef(DClinearRegression)[2], digits = 5), "\n DC R^2 = ",
                      round(DClinearRegression$r.squared, digits = 5),",\n adj DC R^2 = ", 
                      round(DClinearRegression$adj.r.squared, digits = 5), "\n DC p-val = ", 
                      pf(c[1], c[2], c[3], lower=FALSE), sep = ""), cex = 2/3, col = "blue")
  text(250, (1.6/3) * halfrange, paste(" Colon slope = ",round(coef(colonlinearRegression)[2], digits = 5), "\n Colon R^2 = ",
                      round(colonlinearRegression$r.squared, digits = 5),",\n adj Colon R^2 = ", 
                      round(colonlinearRegression$adj.r.squared, digits = 5), "\n Col. p-val = ", 
                      pf(d[1], d[2], d[3], lower=FALSE), sep = ""), cex = 2/3, col = "green")
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.copy(pdf, paste("Glp2r vs. ", ylab, ".pdf", sep = ""))
  dev.off()
}
  

