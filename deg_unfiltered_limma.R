## 201110_A00558_0095_BHVNWGDMXX 
# /Volumes/king/piggybac/mouse2020

library (limma)
library (edgeR)
library (openxlsx)


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
## remove a pair of samples

# remove all samples from batch 6
pheno <- pheno[grep ("_6", pheno$sample, invert=TRUE), ]
# remove an additional sample: MsCTX_PGBD5_sh8_7_S23
pheno <- pheno[grep ("PGBD5_sh8_7", pheno$sample, invert=TRUE), ]
pheno

a <- a[ ,colnames (a) %in% pheno$sample]
idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]

stopifnot (colnames (a) == pheno$sample)



x <- DGEList(counts=a) 

# filter for cpm > 6 in at least 3 samples
isexpr <- rowSums(cpm(x) > 6) >= 3
x <- x[isexpr, ]
dim (x$counts)
#  11846     9


## paired limma test (the paired factor is treated as a batch factor)

celltype <- factor (pheno$genotype)
batch <- factor (gsub (".*_", "", gsub ("_S.*", "", colnames (x$counts)) ))
batch <- factor (as.numeric (gsub ("\\..*", "", batch)))
batch
# Levels: 4 7 8

design <- model.matrix (~ 0 + celltype + batch) 
colnames (design) <- gsub ("celltype", "", colnames (design))

contr.matrix <- makeContrasts (shpgbvsctrl = shPGBD5-CONTROL, shpgbdvsshctrl = shPGBD5-shCONTROL, levels = colnames(design))

v <- voomWithQualityWeights(x, design=design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, trend=TRUE)

res1 <- topTable(efit,coef=1,sort.by="P", n= Inf)
res2 <- topTable(efit,coef=2,sort.by="P", n= Inf)

res1[row.names (res1) == anno[anno$gene_name == "Pgbd5", ]$gene_id, ]
#                          logFC  AveExpr         t      P.Value  adj.P.Val
#ENSMUSG00000050751.16 -1.181046 6.163766 -4.949085 0.0007322531 0.01160692

res2[row.names (res2) == anno[anno$gene_name == "Pgbd5", ]$gene_id, ]
#                          logFC  AveExpr         t      P.Value  adj.P.Val
#ENSMUSG00000050751.16 -1.580888 6.163766 -7.013265 5.456355e-05 0.02951769

boxplot (res1$logFC)
abline (h=0)
abline (h=1)
abline (h=-1)

boxplot (res2$logFC)
abline (h=0)
abline (h=1)
abline (h=-1)


colnames (res1) <- paste (colnames (res1), "shpgbvsctrl", sep=".")
colnames (res2) <- paste (colnames (res2), "shpgbvsshctrl", sep=".")
resall <- merge (res1, res2, by="row.names")
colnames (resall)[1] <- "gene_id"
res <- merge (resall, anno, by="gene_id")

res[res$gene_name == "Pgbd5", ]


## see individual paired log2 fold changes. These are not normalized values
## Pairwise comparisons of batches 4 and 8 (these are the 2 working experiments) !!!
## Not sure why they are not against the CONTROL originally !!!
norm.exprs <- v$E
norm.dif1 <- data.frame (D04= norm.exprs[ ,"PGBD5_sh8_4"] -  norm.exprs[ ,"shCONTROL_4"])
norm.dif2 <- data.frame (D081= norm.exprs[ ,"PGBD5_sh8_8.1"] -  norm.exprs[ ,"shCONTROL_8"])
norm.dif3 <- data.frame (D082= norm.exprs[ ,"PGBD5_sh8_8.2"] -  norm.exprs[ ,"shCONTROL_8"])

#norm.dif1 <- data.frame (D04= norm.exprs[ ,"PGBD5_sh8_4"] -  norm.exprs[ ,"CONTROL_4"])
#norm.dif2 <- data.frame (D081= norm.exprs[ ,"PGBD5_sh8_8.1"] -  norm.exprs[ ,"CONTROL_8"])
#norm.dif3 <- data.frame (D082= norm.exprs[ ,"PGBD5_sh8_8.2"] -  norm.exprs[ ,"CONTROL_8"])

norm.dif <- cbind (norm.dif1, norm.dif2, norm.dif3)


norm.dif$consistent <- "No"
norm.dif$consistent [apply (norm.dif[ ,1:3], 1, function (x) all (x > 0) )] <- "Up"
norm.dif$consistent [apply (norm.dif[ ,1:3], 1, function (x) all (x < 0) )] <- "Down"
table (norm.dif$consistent)
#Down   No   Up 
#3964 5267 2791 

res <- merge (res, norm.dif, by.x="gene_id", by.y="row.names")
res <- res[ ,c(1:17, 19:22, 18)]
res <- res[order (res$adj.P.Val.shpgbvsctrl), ]

table (res$adj.P.Val.shpgbvsctrl < 0.05 & res$consistent != "No")
#FALSE  TRUE 
# 9140  2882  
table (res$adj.P.Val.shpgbvsshctrl < 0.05 & res$consistent != "No")
#FALSE  TRUE 
#11939    83  

res[res$gene_name == "Pgbd5", ]

write.xlsx (res, "deg_unfiltered_piggybac_mouse_shRNA_limma_new_pipeline.xlsx", rowNames=F)









