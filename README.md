# `DMRNAseq`

This workflow is is intended for analysis of direct RNA data. It combines both the differential methylation and differential expression pipelines into one. The structure of this combined workflow is as follows:
1. Alignment with `minimap2`
2. Transcript counting and read to feature assignment using `HTSeq`    
    a. Transcript Counting used for differential expression     
    b. Read to feature assignment is used for differential methylation     
3. Differential Expression Using `DESeq2`
4. Differential Methylation of Transcripts from `HTSeq` Assignments