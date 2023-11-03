# InfiniumEPICMethylation.hg19
Gene annotation for Illumina Infinium EPIC DNA methylation probes.

## Introduction
This work is to create the annotation of CpG sites to genes. 

## Algorithm
The annotations for CpG sites on the Illumina Methylation EPIC arrays were obtained from the "RnBeads" R package (Mueller, 2019). Each CpG site was associated with the closest transcript, determined by the shortest distance from the CpG site to the transcript start site (TSS). Additionally, relevant gene information from the associated transcript was included. The source of the transcript and gene information is the "EnsDb.Hsapiens.v75" R package (Rainer, 2017). Only the transcripts with biotypes "processed_transcript" and "protein_coding" were used.

The computation of distance considers the orientation of the transcripts. Consequently, the TSS corresponds to the 5' end of a transcript on the "+" strand and to the 3' end on the "-" strand.
