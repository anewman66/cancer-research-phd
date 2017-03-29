#DESeq Testing Script using whole RNA_Seq cohort with relapse information.
#Excludes:

library(BiocInstaller)
?BiocUpgrade
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("BiocUpgrade")
biocLite("DESeq2", dependencies=TRUE)
install.packages("data.table",dependencies=TRUE)
library("data.table")
library("DESeq")
library("DESeq2")

setwd("H:/Alex/RNA_Seq")

################################################################################################################

#First, loading data in is different, need to use DESeqDataSetFromHTSeqCount rather than newCountDataSet...
#Also add design=condition
BNHL.subtype.countdesign = read.csv("H:/Alex/RNA_Seq/BL_vs_DLBCL.txt", sep='\t')
?DESeqDataSetFromHTSeqCount
DESeqDataSet.Subtype = DESeqDataSetFromHTSeqCount(BNHL.subtype.countdesign, directory="H:/Alex/RNA_Seq/HTSeq", design= ~ subtype)
#The above line created a value "first_time"??
DESeqDataSet.Subtype = DESeqDataSet.Subtype[ rowSums(counts(DESeqDataSet.Subtype)) > 1, ] #pre-filt, can exclude
#First, set the factor levels correctly, or the package will treat them alphabetically.
#DESeqDataSet.relapse$condition <- factor(DESeqDataSet.relapse$condition, levels ="Non-Relapse","Relapse")
#factor level above lead to errors so instead of trying it, go straight to relevel below.
head(DESeqDataSet.Subtype)
DESeqDataSet.Subtype$subtype <- relevel(DESeqDataSet.Subtype$subtype, ref="DLBCL")
DESeqDataSet.Subtype <- DESeq(DESeqDataSet.Subtype)
result.subtype <- results(DESeqDataSet.Subtype)
head(result.subtype)
write.csv(result.subtype,"result1_subtype.txt")
#As far as I can tell, the rest should be as above

#Here going to try another method of generating the padj values:
bonferroni.pvals= p.adjust(result.subtype$pvalue, method='bonferroni')
write.csv(bonferroni.pvals, 'bonferroni_pvals_subtype.txt')

#Can also order our results by smallest adj.pvalue:
#pvalordered.result.subtype <- result.subtype[order(result.subtype$padj),]
#pvalordered2.result.subtype <- result.subtype[order(result.subtype$pvalue),]
#head(pvalordered2.result.subtype)
#Summarise the results:
#summary(result.subtype) #at the moment this gives no up- or down-regulated genes?!
#sum(result.subtype$padj < 0.1, na.rm=TRUE)
#sum(result.subtype$pvalue < 0.1, na.rm=TRUE)
#?plotMA
#head(result.subtype)
#head(pvalordered.result.subtype)
#is.data.frame(DESeqDataSet.Subtype)
#plotMA(DESeqDataSet.Subtype)
#plotMA(data.frame(result.subtype))
#hist(result.subtype$padj, breaks=100, col="skyblue", border="slateblue", main="")
#hist(result.subtype$pvalue, breaks=100, col="skyblue", border="slateblue", main="")


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
gene_list1=read.csv("H:/Alex/RNA_Seq/result1_subtype_bpadj.txt", sep='\t')
#Below, p05 denotes p=0.05. Choosing all rows with b.padj below 0.05
subtype.gene.listp05 = gene_list1[(gene_list1$b.padj < 0.05),]
ensembl.gene.listp05 = subtype.gene.listp05$X
sig.diff.genes=getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
      filters = 'ensembl_gene_id',
      values = ensembl.gene.listp05 ,
      mart=ensembl)

write.csv(sig.diff.genes, 'significant_genes_subtype_p05.txt')
?write.csv
