# cancer-research-phd
A collection of scripts used in my analysis of DNA- and RNA-Seq
Some are R scripts, some are for submitting to the unix cluster for larger jobs. I will eventually be working on scripting these together so that they are run in one workflow rather than separate scripts.

Main functions:

Full RNA-Seq QC and alignment pipeline - see the BSU presentation slide I made for reference.
  * FASTQC on fastq files
  * Genome generation with STAR
  * Alignment to reference genome using STAR
  * sam to bam conversion
  * Sort bam files by genomic co-ordinates for indexing
  * Index bam files
  * htseq for transcript quantification
