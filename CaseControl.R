#R
options(warn=-1)

cases <- read.table("DDDChildren+BWS/DDDChildren+BWS.metaScores.genes.txt", header=TRUE)
controls <- read.table("1958BC/DDDChildren+BWS+1958BC.metaScores.genes.txt", header=TRUE)
cases$Individuals <- NULL
controls$Individuals <- NULL
cases$GeneCountsRef <- 302 - cases$GeneCounts
controls$GeneCountsRef <- 999 - controls$GeneCounts
data1 <- merge(cases, controls, by.x="Gene", by.y="Gene", all.x=T, all.y=T)
data2 <- data1[complete.cases(data1), ]
data2$P <- 1
for (i in 1:nrow(data2)) {
    data2$P[i] <- (chisq.test(matrix(c(data2$GeneCounts.x[i], data2$GeneCounts.y[i], data2$GeneCountsRef.x[i], data2$GeneCountsRef.y[i]), byrow=T, 2, 2)))$p.value
}
for (i in 1:nrow(data2)) {
    data2$Pfdr[i] <- p.adjust(data2$P[i], method="fdr", n=nrow(data2))
}
data3 <- data2[order(data2$Pfdr),]
write.table(data3, file="Gene.CaseControl.txt", sep="\t", quote=F, row.names=F, na="-9")

######################################################################
cases <- read.table("DDDChildren+BWS/DDDChildren+BWS.metaScores.variants.txt", header=TRUE)
controls <- read.table("1958BC/DDDChildren+BWS+1958BC.metaScores.variants.txt", header=TRUE)
cases$Individuals <- NULL
controls$Individuals <- NULL
cases$VariantCountsRef <- 302 - cases$VariantCounts
controls$VariantCountsRef <- 999 - controls$VariantCounts
data1 <- merge(cases, controls, by.x="Pos", by.y="Pos", all.x=T, all.y=T)
data2 <- data1[complete.cases(data1), ]
data2$P <- 1
for (i in 1:nrow(data2)) {
    data2$P[i] <- (chisq.test(matrix(c(data2$Alt_Cases[i], data2$Alt_Controls[i], data2$Ref_Cases[i], data2$Ref_Controls[i]), byrow=T, 2, 2)))$p.value
}
for (i in 1:nrow(data2)) {
    data2$Pfdr[i] <- p.adjust(data2$P[i], method="fdr", n=nrow(data2))
}
data3 <- data2[order(data2$Pfdr),]
write.table(data3, file="Variant.CaseControl.txt", sep="\t", quote=F, row.names=F, na="-9")


