#DESeq Testing Script using 10_10065 F1 and F2 to test differential expression between relapse and non-relapse.
library(BiocInstaller)
?BiocUpgrade
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("DESeq2")
library("DESeq")
library("DESeq2")
??DESeq2
# Making a descriptive dataframe
#If you have used htseq-count (see above) to obtain your read counts, you can use the funtion newCountDataSet-
#FromHTSeqCount to obtain a CountDataSe directly from the count files produced by htseq-count. To this end, produce
#a data frame, similar to pasillaDesign constructed above, with the sample names in the first column, the file names
#in the second column, and the condition in the third column (or, if you have more than one factor in your design
#matrix, all these factors in the third and folling columns), and pass this data frame to the function. See the help page
#(\?newCountDataSetFromHTSeqCount") for further details.
?read.csv
??newCountDataSetFromHTSeqCount
countdesign = read.csv("H:/Alex/RNA_Seq/10_10065_Count_Design.txt", header=TRUE, sep= "\t")
cds = newCountDataSetFromHTSeqCount(countdesign, directory= "H:/Alex/RNA_Seq/")
#We have now loaded in the two count files using newCountDataSetFromHTSeqCount 
#Next step is normalisation using "estimateSizeFactors" to get size factors to remove the impact of 
#differing sequencing depth on the samples.
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized=TRUE)) #shows the values of each count divided by samples respective sizeFactor
cds = estimateDispersions(cds, method='blind') #estimating variance. As no repeats at the moment for counts,
#the method='blind' function is needed.
str(fitInfo(cds))
plotDispEsts(cds) #draws the plot comparing dispersion to mean of normalised counts
head(fData(cds)) #to be used for quality control alongside plot... not entirely sure what really to look for
res = nbinomTest(cds,"Diagnostic","Relapse") #as in the test we had one of each condition and had to use 
#method='blind', this will not yield any useful results. 
head(res)
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

#Now to start again with my own data
BNHL.countdesign = read.csv("H:/Alex/RNA_Seq/
