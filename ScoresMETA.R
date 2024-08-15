
arguments = commandArgs(trailingOnly = TRUE)
args <- commandArgs(TRUE)
data0 <- read.table(file=arguments[1], header=TRUE, stringsAsFactors=FALSE, na.strings="-9", sep="\t", row.names=NULL, quote=NULL, fill=T)
variants <- read.table(file=arguments[2], header=T)
genes <- read.table(file=arguments[3], header=T)
total <- as.numeric(arguments[4])
genesbed <- read.table("genesbed.bed", header=F)

data0$Chr_Start_End <- paste(data0$Chr, data0$Start, data0$End, sep="_")
data1 <- merge(data0, variants, by.x="Chr_Start_End", by.y="Pos", all.x=T)
data2 <- merge(data1, genes, by.x="Gene.refGene", by.y="Gene", all.x=T)
data2$GenePerct <- data2$GeneCounts / data2$TotalGeneCounts
data2$VariantPerct <- data2$VariantCounts / data2$VariantTotalCounts
data2$metaScore <- (as.numeric(data2$Score) * (1 + data2$GenePerct) * (1 + data2$VariantPerct)) / 4

data3 <- data2[,c("Chr.x", "Start.x", "End.x", "cytoBand", "Gene.refGene", "GeneCounts", "TotalGeneCounts", "GenePerct", "Func.refGene", "ExonicFunc.refGene", "avsnp144.x", "VariantCounts", "VariantTotalCounts", "VariantPerct", "Ref", "Alt", "AA1", "AA2", "AA_Pos", "Lenght", "TranscriptID", "Exon", "Interpro_domain", "X1000g2015aug_all", "ExAC_ALL", "esp6500siv2_all", "metaScore", "Score", "TruncScore", "HomozygScore", "XAcorrect", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "VEST3_score", "MutationAssessor_score_Converted", "GERP.._RS_Converted", "phyloP7way_vertebrate_score", "phyloP20way_mammalian_score", "SIFT_scoreInv", "fathmm.MKL_coding_score", "LRT_score", "MutationTaster_score", "FATHMM_score_Converted", "PROVEAN_score_Converted", "CADD_raw_Converted", "DANN_score", "MetaSVM_score_Converted", "MetaLR_score", "integrated_fitCons_score_Converted", "integrated_confidence_value_Converted", "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds_Converted", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "dpsi_max_tissue_Converted", "X1000g2015aug_afr", "X1000g2015aug_eas", "X1000g2015aug_eur", "X1000g2015aug_amr", "X1000g2015aug_sas", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "esp6500siv2_aa", "esp6500siv2_ea", "nci60", "cosmic70", "CLINSIG", "CLNDBN", "CLNACC", "CLNDSDB", "CLNDSDBID")]
colnames(data3)[1] <- "Chr"
colnames(data3)[2] <- "Start"
colnames(data3)[3] <- "End"
colnames(data3)[11] <- "avsnp144"

data4 <- data3[order(-data3$metaScore),]

genes2 <- merge(genes, genesbed, by.x="Gene", by.y="V4", all.x=T, all.y=F)
colnames(genes2)[5] <- "Chr"
colnames(genes2)[6] <- "Start"
colnames(genes2)[7] <- "End"
genes2$V5 <- NULL
genes3 <- genes2[,c("Chr", "Start", "End", "Gene", "GeneCounts", "TotalGeneCounts", "Individuals")]

write.table(data4, file="metaScores.txt", sep="\t", quote=F, row.names=F, na="-9")
write.table(genes3, file="genes3.txt", sep="\t", quote=F, row.names=F, na="-9")
