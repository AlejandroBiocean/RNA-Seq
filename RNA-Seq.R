#1. Read count genes data
unSortedCounts <- read.table(file = "C:\\Users/alexb/Desktop/Bioinformática/Algoritmos para el análisis de secuencias/TP3_RNASeq_Counts.txt", sep ="\t", stringsAsFactors = TRUE, header = TRUE, row.names = "Geneid")
countsFile <- unSortedCounts[ , order(names(unSortedCounts))]


#2. Read and order Sample Sheet
unSortedSheet <- read.table(file = "C:\\Users/alexb/Desktop/Bioinformática/Algoritmos para el análisis de secuencias/TP3_RNASeq_SampleSheet.txt", sep ="\t", stringsAsFactors = FALSE, header = TRUE)
targetsFile <- unSortedSheet[order(unSortedSheet$Sample),]
rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample 


#3. Load "DESeq2" and create DESeqDataSet
library("DESeq2")
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsFile, colData = targetsFile, design = ~ Strain + Treatment)


#4. DESeq2 analysis. 
dds <- DESeq(ddsFullCountTable) 
resultsNames(dds)


#5. Create Control & Treatment objects. 
factor_condition <- as.factor(dds$condition)
Control <- factor_condition == "Control"
Treatment <- factor_condition == "Treatment"


#6. Establish control group. 
ddsFullCountTable$Treatment <- relevel(ddsFullCountTable$Treatment, "Control")


#7. dds results. 
res <- results(dds)
head(res)  
summary(res)


#8. Dispersion plot. 
plotDispEsts(dds)


#9. Differentially expressed genes number. 
is.na(res)
sum(is.na(res))
na.omit(res)

genes_differentially_expr <- sum(res$padj < 0.05, na.rm = TRUE) 
genes_differentially_expr


#10. Order differentially expressed genes by p-value adjusted (padj). 
res_orden <- res[order(na.omit(res$padj)), ]
write.table(res_orden, file = "analysis_differentially_genes_results.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

differentially_expression_data <- read.table("analysis_differentially_genes_results.tsv", header = TRUE, sep = "\t")
differentially_expression_data 
 


#11. Generate MAplots with different rates. 
plotMA(res, ylim=c(-2,2), alpha = 0.05)
plotMA(res, ylim=c(-2,2), alpha = 0.001)
plotMA(res, ylim=c(-2,2), alpha = 0.000001) 


