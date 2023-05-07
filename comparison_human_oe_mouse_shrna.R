library (openxlsx)

sh <- read.xlsx ("/Volumes/king/piggybac/mouse2020/deg_unfiltered_piggybac_mouse_shRNA_limma_new_pipeline.xlsx")
sh <- sh[ ,c("gene_name", "logFC.shpgbvsctrl", "adj.P.Val.shpgbvsctrl", "logFC.shpgbvsshctrl", "adj.P.Val.shpgbvsshctrl", "logFC.shctrolvsctrl", "adj.P.Val.shctrolvsctrl")]
colnames (sh) <- paste (colnames (sh), "sh.mouse", sep=".")
sh$gene_name.sh.mouse <- toupper (sh$gene_name.sh.mouse)
head (sh)

oe <- read.xlsx ("/Volumes/king/piggybac/human2020/deg_unfiltered_piggybac_overexpression_limma_new_pipeline.xlsx")
oe <- oe[ ,c("gene_name", "logFC", "adj.P.Val")]
colnames (oe) <- paste (colnames (oe), "oe.human", sep=".")
head (oe)

comp <- merge (sh, oe, by.x="gene_name.sh.mouse", by.y="gene_name.oe.human")
comp$significance1 <- "No"
comp$significance1 [comp$adj.P.Val.shpgbvsctrl.sh.mouse < 0.05 | comp$adj.P.Val.shpgbvsshctrl.sh.mouse < 0.05] <- "Yes"
comp$significance2 <- "No"
comp$significance2[comp$significance1 == "Yes" & comp$adj.P.Val.shctrolvsctrl.sh.mouse > 0.05] <- "Yes"
comp$significance.overall <- "No"
comp$significance.overall[comp$significance2 == "Yes" & comp$adj.P.Val.oe.human < 0.05] <- "Yes"

comp$trend.mouse <- "No"
comp$trend.mouse [comp$logFC.shpgbvsctrl.sh.mouse < 0 | comp$logFC.shpgbvsshctrl.sh.mouse <0] <- "Down"
comp$trend.mouse [comp$logFC.shpgbvsctrl.sh.mouse > 0 | comp$logFC.shpgbvsshctrl.sh.mouse >0] <- "Up"
comp$trend.overall <- "No"
comp$trend.overall[comp$trend.mouse == "Down" & comp$logFC.oe.human > 0 & comp$significance.overall == "Yes"] <- "RC1"
comp$trend.overall[comp$trend.mouse == "Up" & comp$logFC.oe.human < 0 & comp$significance.overall == "Yes"] <- "RC2"
table (comp$trend.overall)
#  No  RC1  RC2 
#8541  337   85 

comp.s <- comp[comp$trend.overall != "No", ]
write.xlsx (comp.s, "trend.xlsx", rowNames=F)




