library(RnBeads)
library(RnBeads.hg19)
library(EnsDb.Hsapiens.v75)

rnb.anno.epic <- rnb.get.annotation("probesEPIC")

### mergen them into one
rnb.anno.epic.all <- rnb.anno.epic[[1]]

for(i in 2:length(rnb.anno.epic)) {
  rnb.anno.epic.all <- c(rnb.anno.epic.all, rnb.anno.epic[[i]])
}


edb <- EnsDb.Hsapiens.v75
allgenes <- genes(edb)

transcripts <- transcripts(edb)
transcripts <- transcripts[transcripts$tx_biotype %in% c("processed_transcript", "protein_coding")]
transcripts <- transcripts[grepl("^[1-9]$|^1[0-9]$|^2[0-2]$|^X$|^Y$", seqnames(transcripts))]
seqlevels(transcripts) <- as.character(unique(seqnames(transcripts)))
transcripts <- renameSeqlevels(transcripts, paste0("chr", seqlevels(transcripts)))
genome(transcripts) <- "hg19"

############# modifying the transcripts and running the distanceToNearest
transcripts.mod <- transcripts

end(transcripts.mod[which(strand(transcripts.mod) == "+")]) = start(transcripts.mod[which(strand(transcripts.mod) == "+")]) 
start(transcripts.mod[which(strand(transcripts.mod) == "-")]) = end(transcripts.mod[which(strand(transcripts.mod) == "-")]) 

nearest_transcripts <- distanceToNearest(rnb.anno.epic.all, transcripts.mod, ignore.strand = TRUE)
closest_transcript_ids <- names(transcripts.mod)[nearest_transcripts@to]
closest_transcript_dist <- nearest_transcripts@elementMetadata$distance
v.distTxStart <- closest_transcript_dist
names(v.distTxStart) <- closest_transcript_ids

### add annotation to the grange object
rnb.anno.epic.all.new <- rnb.anno.epic.all
mcols(rnb.anno.epic.all.new) <- data.frame(mcols(rnb.anno.epic.all.new),
                                           tx_id=names(v.distTxStart),
                                           tx_chr=seqnames(transcripts[names(v.distTxStart)]),
                                           tx_start=start(transcripts[names(v.distTxStart)]),
                                           tx_end=end(transcripts[names(v.distTxStart)]),
                                           tx_width=width(transcripts[names(v.distTxStart)]),
                                           tx_strand=strand(transcripts[names(v.distTxStart)]),
                                           distTxStart=v.distTxStart,
                                           gene_id=transcripts[names(v.distTxStart)]$gene_id,
                                           gene=allgenes[transcripts[names(v.distTxStart)]$gene_id]$gene_name,
                                           gene_biotype=allgenes[transcripts[names(v.distTxStart)]$gene_id]$gene_biotype,
                                           seq_coord_system=allgenes[transcripts[names(v.distTxStart)]$gene_id]$seq_coord_system,
                                           symbol=allgenes[transcripts[names(v.distTxStart)]$gene_id]$symbol
                                           )

mcols(rnb.anno.epic.all.new)$entrezid <- allgenes[transcripts[names(v.distTxStart)]$gene_id]$entrezid

### using the findOverlaps methods to find the overlapped genes (only "protein coding genes")
overlaps_genes_cpgs <- findOverlaps(rnb.anno.epic.all, allgenes.protein, ignore.strand=TRUE)
df.overlaps_genes_cpgs <- as.data.frame(overlaps_genes_cpgs)
df.overlaps_genes_cpgs$gene_id <- names(allgenes.protein[df.overlaps_genes_cpgs$subjectHits])
df.overlaps_genes_cpgs$gene <- allgenes.protein[df.overlaps_genes_cpgs$subjectHits]$gene_name

df.overlaps_genes_cpgs <- df.overlaps_genes_cpgs[!duplicated(df.overlaps_genes_cpgs[c("queryHits", "gene")]), ]

df.overlaps_genes_cpgs.agg <- aggregate(gene_id ~ queryHits, df.overlaps_genes_cpgs, paste, collapse = ";")
df.overlaps_genes_cpgs.agg.gene <- aggregate(gene ~ queryHits, df.overlaps_genes_cpgs, paste, collapse = ";")


v.ifTxIN <- rep(NA, length(rnb.anno.epic.all))
v.ifTxIN[df.overlaps_genes_cpgs.agg$queryHits] <- df.overlaps_genes_cpgs.agg$gene_id

v.ifTxIN.gene <- rep(NA, length(rnb.anno.epic.all))
v.ifTxIN.gene[df.overlaps_genes_cpgs.agg.gene$queryHits] <- df.overlaps_genes_cpgs.agg.gene$gene

mcols(rnb.anno.epic.all.new)$gene_id_in <- v.ifTxIN
mcols(rnb.anno.epic.all.new)$gene_in <- v.ifTxIN.gene

### the final annotation is stored in a data frame
df.rnb.anno.epic.all <- as.data.frame(rnb.anno.epic.all.new)
