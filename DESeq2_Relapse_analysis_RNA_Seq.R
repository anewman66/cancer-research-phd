#DESeq Testing Script using whole RNA_Seq cohort with relapse information.
#Excludes:

library(BiocInstaller)
?BiocUpgrade
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("DESeq2")
library("DESeq")
library("DESeq2")


################################################################################################################

#First, loading data in is different, need to use DESeqDataSetFromHTSeqCount rather than newCountDataSet...
#Also add design=condition
BNHL.relapse.countdesign = read.csv("H:/Alex/RNA_Seq/Relapse_Non-relapse_CountDesign.txt", sep='\t')
?DESeqDataSetFromHTSeqCount
View(BNHL.relapse.countdesign)
DESeqDataSet.relapse = DESeqDataSetFromHTSeqCount(BNHL.relapse.countdesign, directory="H:/Alex/RNA_Seq/HTSeq", design= ~ condition)
#The above line created a value "first_time"??
DESeqDataSet.relapse = DESeqDataSet.relapse[ rowSums(counts(DESeqDataSet.relapse)) > 1, ] #pre-filt, can exclude
#First, set the factor levels correctly, or the package will treat them alphabetically.
#DESeqDataSet.relapse$condition <- factor(DESeqDataSet.relapse$condition, levels ="Non-Relapse","Relapse")
#factor level above lead to errors so instead of trying it, go straight to relevel below.
View(DESeqDataSet.relapse)

#DESeqDataSet.relapse$condition <- relevel(DESeqDataSet.relapse$condition, ref="NonRelapse")
DESeqDataSet.relapse <- DESeq(DESeqDataSet.relapse)
result.relapse <- results(DESeqDataSet.relapse)
resultvector2 <- c("condition","NonRelapse","Relapse")
result.relapse1 <- results(DESeqDataSet.relapse, contrast=resultvector2)
head(result.relapse)
setwd("H:/Alex/RNA_Seq")
write.csv(result.relapse,"result1_relapse.txt")
#As far as I can tell, the rest should be as above
#Can also order our results by smallest adj.pvalue:
pvalordered.result.relapse <- result.relapse[order(result.relapse$padj),]
pvalordered2.result.relapse <- result.relapse[order(result.relapse$pvalue),]
head(pvalordered2.result.relapse)
#Summarise the results:
summary(result.relapse) #at the moment this gives no up- or down-regulated genes?!
sum(result.relapse$padj < 0.1, na.rm=TRUE)
sum(result.relapse$pvalue < 0.1, na.rm=TRUE)
?plotMA
head(result.relapse)
is.data.frame(DESeqDataSet.relapse)
plotMA(DESeqDataSet.relapse)
plotMA(data.frame(result.relapse))
hist(result.relapse$pvalue, breaks=100, col="skyblue", border="slateblue", main="")
hist(result.relapse$padj, breaks=100, col="skyblue", border="slateblue", main="")

#Using Bonferroni padj
bonferroni.pvals_relapse= p.adjust(result.relapse$pvalue, method='bonferroni')
write.csv(bonferroni.pvals_relapse, 'bonferroni_pvals_relapse.txt')

#Now how do I get my gene list? 
biocLite('biomaRt')
library('biomaRt')
??biomaRt
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)

#Need to introduce 3 features for getBM(), filters, attributes and values. Also need to specify mart, as we have.
filters=listFilters(ensembl)
head(filters)
attributes=listAttributes(ensembl)
head(attributes)
gene_list_relapse=read.csv("H:/Alex/RNA_Seq/result1_relapse_bpadj.txt", sep='\t')
#Below, p05 denotes p=0.05. Choosing all rows with b.padj below 0.05
subtype.gene.listp05 = gene_list1[(gene_list1$b.padj < 0.05),]
ensembl.gene.listp05 = subtype.gene.listp05$X
sig.diff.genes=getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
                     filters = 'ensembl_gene_id',
                     values = ensembl.gene.listp05 ,
                     mart=ensembl)

write.csv(sig.diff.genes, 'significant_genes_subtype_p05.txt')
?write.csv














