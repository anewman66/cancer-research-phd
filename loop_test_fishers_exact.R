library(BiocInstaller)

setwd("N:/Alex Newman/Genomics paper NEWMAN/Post_Path_Review/Paper_Writing")

counts.pBL.aBL = read.table("pBL_aBL_GISTIC_counts_16032017.txt", sep = "\t",header = TRUE, row.names = 1)
#counts.pBL.aBL = read.table("pBL_aBL_inc_counts_gistic2.txt", sep = "\t",header = TRUE, row.names = 1)

#From Vikki
#Two groups:
#Group 1 has 63 samples - 50 do not have the CNA, 13 have the CNA 
#Group 2 has 70 samples - 65 do not have the CNA, 5 have the CNA

MIR17HG_rel_test2 = matrix(c(9,6,57,23), nrow=2)
chisq.test(MIR17HG_rel_test2)
fisher.test(MIR17HG_rel_test2)

TP53_rel_test = matrix(c(5,10,65,15), nrow=2)
chisq.test(TP53_rel_test)
fisher.test(TP53_rel_test)

counts.pDLBCL.aDLBCL = read.table("pDLBCL_vs_aDLBCL_fisher_test_counts.txt", sep = "\t", header=TRUE, row.names = 1)

#Other tactic using apply

alltables = counts.pBL.aBL

pBL_aBL_pvalues = apply(alltables,1, function(x) fisher.test(matrix(x,nr=2))$p.value)
print(pBL_aBL_pvalues)

pBL_aBL_pvalues=data.frame(pBL_aBL_pvalues)
write.table(pBL_aBL_pvalues, file = "pvalues_pBL_aBL_16032017.txt", sep = '\t')

alltables1 = counts.pDLBCL.aDLBCL

pDLBCL_aDLBCL_pvalues = apply(alltables1,1, function(x) fisher.test(matrix(x,nr=2))$p.value)
print(pDLBCL_aDLBCL_pvalues)

#THIS ABSOLUTELY WORKS.
#Edited for paeds_166 files.




apoptosis_for_R_p_vals = apply(apoptosis_for_R_sepd,1, function(x) chisq.test(matrix(x,nr=2))$p.value)
print(apoptosis_for_R_p_vals)

apoptosis_for_R_p_vals=data.frame(apoptosis_for_R_p_vals)
write.table(pBL_aBL_pvalues, file = "pvalues_pBL_aBL_16032017.txt", sep = '\t')