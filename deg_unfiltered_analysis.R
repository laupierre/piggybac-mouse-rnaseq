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
# remove a pair of samples
pheno <- pheno[!pheno$sample %in% c("PGBD5_sh8_4","shCONTROL_4"), ]
pheno

a <- a[ ,colnames (a) %in% pheno$sample]
idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]

stopifnot (colnames (a) == pheno$sample)


dds <- DESeqDataSetFromMatrix(countData = round (a), colData = pheno, design = ~ genotype)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "shPGBD5", "shCONTROL"))

res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Pgbd5", ] 

write.xlsx (res, "piggybac_shPB_vs_shCTR_mouse_in-vitro.xlsx", rowNames=F)



## heatmap plot
# Neurod1, Neurod2, Dcx, RhoA, RhoB, Rac1
#select <- c("ENSMUSG00000034701.10", "ENSMUSG00000038255.7", "ENSMUSG00000031285.15", "ENSMUSG00000007815.14", "ENSMUSG00000054364.6", "ENSMUSG00000001847.15", "ENSMUSG00000050751.16")

select <- c("Pgbd5", "Pax6", "Tbr2", "Hes1", "Hes5" ,"Nes", "Ccnd2", "Cdk2", "Neurod1", "Neurod2", "Neurod6", "Sox11", "Fezf2", "Ctip2", "Satb2", "Dcx", "Rnd2", "Rac1", "Rac2", "Rac3", "Rhoa", "Rhob", "Cdc42")
select <- annot[annot$gene_name %in% select, ]$Geneid

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Heatmap plot.pdf")
pheatmap( log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select, ], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off ()



## heatmap plot of top 75 genes

select <- res$Geneid[1:75]

pdf ("Heatmap 75 top genes plot.pdf")
pheatmap( log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select, ], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
dev.off ()


## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, label=sample)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		geom_text_repel()  + 
		  coord_fixed () 

ggsave ("PCA plot.pdf")









