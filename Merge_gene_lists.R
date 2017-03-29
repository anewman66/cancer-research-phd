### Merging lists of genes with counts and incidence in % from two different sources. One in this case is 
## non-relapse genes and one is relapse genes. Gene lists are not identical but share common genes.
# This practice is to allow for comparison between subgroups. 
getwd()
setwd("N:/Alex Newman/Exome Sequencing Analysis/MuTect1_Variant")
m=merge(CN_Samples,Exome_Seq_Samples,by="Name",all=TRUE)
#Needed to set all=TRUE because otherwise R only outputs the rows that were merged.
m2=merge(m,RNA_Seq_Samples,by="Name",all=TRUE)

write.csv(m2, file = "CN_RNA_Exome_Samples.txt")
