# by neelendu

# objective: deploy Random Forests to identify ENS/BAM genes most discriminatory for (i) microbiota, (ii) biogeography, (iii) diet, and (iv) each combination of these variables.
#   ENS = enteric nervous system
#   BAM = bile acid metabolism
# microbiota = germ-free, BSH-high, or BSH-combined
# biogeography = small intestine, proximal colon, or distal colon
# diet = turmeric or non-turmeric

# load requisite libraries
#source("~/Dey_lab/scripts/library_general.R")
#library(pheatmap)
library(randomForest)
sem <- function (x)  {
  return (sd(x, na.rm=T) / (sqrt(length(na.omit(x)))))
}

# read in data files
dat.norm = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/Normalized_CodeString_data.txt", header = T, sep = "\t")
dat.norm = dat.norm[-c(10:12, 22:24, 34:36, 40:42),]
#dat.norm = dat.norm[-c(1:3, 10:15, 22:27, 34:42),]
dat.map = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/map.txt", header = T, sep = "\t")
dat.map$Microbiota = as.character(dat.map$Microbiota)
dat.map$Microbiota[c(4:9, 16:21, 28:33, 43:48)] = "colonized"
dat.map$Microbiota = factor(dat.map$Microbiota)
dat.map = dat.map[-c(10:12, 22:24, 34:36, 40:42),]
#dat.map = dat.map[-c(1:3, 10:15, 22:27, 34:42),]
dat.map$Microbiota = factor(dat.map$Microbiota)
dat.genes = read.table(file = "file:///C:/Users/5p336/Documents/Research/First Project/gene_table.txt", header = T, sep = "\t")

# USER CHECK
table(dat.norm$Sample.ID == dat.map$Sample.ID)  # should be all TRUE

# scratch / RF model, testing on prediction of DCA
n.trees = 10000  # testing with 50-100; use 10000 for final product
n.iterations_RF = 100  # testing with 5-10; use 100 for final product

# define the response vector and matrix of predictors -- USER: pick one of each and comment out the rest
 v.responseVector = as.factor(dat.map$Microbiota); s.response = "microbiota"
# v.responseVector = as.factor(dat.map$Tissue); s.response = "tissue"
# v.responseVector = as.factor(dat.map$Diet); s.response = "diet"
m.predictors = dat.norm

# machine learning (Random Forests)
rf.discrimGenes = randomForest(v.responseVector~., m.predictors, ntree=n.trees, importance=TRUE, na.action=na.omit, proximity=T)
#pdf(paste("Response Variable = ", s.response, ", Trees = ", n.trees, ", Runs = 1", ".pdf", sep = ""), width = 11, height = 8.5)
#varImpPlot(x=rf.discrimGenes, type=1, main = paste("Response Variable = ", s.response, ", Trees = ", n.trees, ", Runs = 1", ".pdf", sep = ""))  # this will show you the results from a single run
dev.off()

# let's run this a bunch more times to determine which genes, on average, are most important to the machine learning model 
imp.iterOne = importance(rf.discrimGenes, type = 1, scale=T)
tab.impSco = matrix(NA, nrow=length(imp.iterOne), ncol=n.iterations_RF)
tab.impSco[,1] = imp.iterOne
rownames(tab.impSco) = rownames(imp.iterOne)
for (i in 2:n.iterations_RF)  {
  rf.discrimGenesDuJour = randomForest(v.responseVector~., m.predictors, ntree=n.trees, importance=TRUE, na.action=na.omit, proximity=T)
  if (table(rownames(importance(rf.discrimGenesDuJour, type = 1, scale=T)) == rownames(tab.impSco)) != nrow(tab.impSco))  {warning("ERROR: order of factors not matched properly to table!")}
  tab.impSco[,i] = importance(rf.discrimGenesDuJour, type = 1, scale=T)
}
df.impSco = data.frame(Gene = rownames(tab.impSco), Mean_Importance_Scores = rowMeans(tab.impSco), SEM_Importance_Scores = apply(X=tab.impSco, 1, sem))
df.impSco = df.impSco[order(-df.impSco$Mean_Importance_Scores),]
df.impSco = df.impSco[-11,]
pdf(paste("Response Variable = ", s.response, ", Trees = ", n.trees, ", Runs = ", n.iterations_RF, ".pdf", sep = ""), width = 11, height = 8.5)
bp = barplot(rev(df.impSco$Mean_Importance_Scores[c(1:24)]), col = "red", horiz = TRUE,  names.arg = rev(df.impSco$Gene[c(1:24)]), las = 1,cex.names = 0.75, 
             xlab = "Gene Importance", ylab = "Gene", xlim = c(-1, 40),
             main = paste("Response Variable = ",s.response, ", # Trees = ", n.trees, ", # Runs = ", n.iterations_RF, sep = ""))
arrows(rev(df.impSco$Mean_Importance_Scores[c(1:24)]) - rev(df.impSco$SEM_Importance_Scores[c(1:24)]), bp, 
       rev(df.impSco$Mean_Importance_Scores[c(1:24)]) + rev(df.impSco$SEM_Importance_Scores[c(1:24)]), bp, code=3, angle=90, length = 0.025)
dev.off()

# consider writing code to save data/figures; if doing so, integrate s.response into the file naming