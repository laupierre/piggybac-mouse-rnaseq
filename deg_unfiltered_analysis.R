## 201110_A00558_0095_BHVNWGDMXX 

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub ("star.IIT_RNM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)


annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]


torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
colnames (a) <- gsub ("star.IIT_X_", "", colnames (a))
a <- a[ ,-1]


pheno <- data.frame (matrix (nrow=dim (a)[2], ncol=2))
colnames (pheno) <- c("sample", "genotype")
row.names (pheno) <- pheno$sample <- colnames (a)
pheno$genotype <- gsub ("_.*", "", pheno$sample)
pheno$genotype [pheno$genotype == "PGBD5"] <- "shPGBD5" 
pheno

idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]

stopifnot (colnames (a) == pheno$sample)


dds <- DESeqDataSetFromMatrix(countData = round (a), colData = pheno, design = ~ genotype)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "shPGBD5", "shCONTROL"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Pgbd5", ] 

write.xlsx (res, "piggybac_shrna_mouse_in-vitro.xlsx", rowNames=F)


## heatmap plot
# Th, DDC, NR4A2, SLC6A3 (known target), SNCA, Pitx3, SLC18A2 (known target, see doi: 10.1111/j.1471-4159.2009.06404.x), Tcf7l2
select <- c("ENSMUSG00000000214.12", "ENSMUSG00000020182.17", "ENSMUSG00000026826.14", "ENSMUSG00000021609.7", "ENSMUSG00000025889.14", "ENSMUSG00000025229.16", "ENSMUSG00000025094.9", "ENSMUSG00000024985.21")

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Heatmap plot.pdf")
pheatmap( log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select, ], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off ()








