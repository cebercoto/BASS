# Load data
arguments = commandArgs(trailingOnly = TRUE)
args <- commandArgs(TRUE)
data0 <- read.table(file=arguments[1], header=TRUE, stringsAsFactors=FALSE, na.strings=".", sep="\t", row.names=NULL, quote=NULL, fill=T)
Transcripts <- read.table(file=arguments[2], header=T, fill=T)
TranscriptSize <- read.table("/data/Resources/Software/ceberSUITE/BASS/All_Transcripts_Lenght.txt", header=T)
# Prepare Data
data0$SIFT_scoreInv <- 1-data0$SIFT_score
data0$MutationAssessor_score_Converted <- ((data0$MutationAssessor_score+5.545)/(5.975+5.545))
data0$GERP.._RS_Converted <- (data0$GERP.._RS+12.3)/(6.17+12.3)
data0$phyloP7way_vertebrate_score <- (data0$phyloP7way_vertebrate+5.172)/(1.062+5.172)
data0$phyloP20way_mammalian_score <- (data0$phyloP20way_mammalian+13.282)/(1.199+13.282)
data0$FATHMM_score_Converted <- (data0$FATHMM_score+16.13)/(10.64+16.13)
data0$PROVEAN_score_Converted <- (data0$PROVEAN_score+14)/(14+14)
data0$CADD_raw_Converted <- (data0$CADD_raw+7.535)/(35.789+7.535)
data0$MetaSVM_score_Converted <- (data0$MetaSVM_score+2.006)/(3.04+2.006)
data0$integrated_fitCons_score_Converted <- (data0$integrated_fitCons_score+0)/(0.84+0)
data0$integrated_confidence_value_Converted <- (data0$integrated_confidence_value+0)/(3+0)
data0$SiPhy_29way_logOdds_Converted <- (data0$SiPhy_29way_logOdds-0.0003)/(37.972-0.0003)
data0$dpsi_max_tissue_Converted <- (data0$dpsi_max_tissue+94.8856)/(88.3205+94.8856)
# Merge Data
data0$Chr_Start_End_Ref_Alt <- paste(data0$Chr, data0$Start, data0$End, data0$Ref, data0$Alt, sep="_")
data1 <- merge(data0, Transcripts, by.x="Chr_Start_End_Ref_Alt", by.y="Chr_Start_End_Ref_Alt", all.x=T)
data2 <- merge(data1, TranscriptSize, by.x="TranscriptID", by.y="RefSeq", all.x=TRUE)
# Truncation score
data2$AA_Pos <- as.numeric(as.character(data2$AA_Pos))
data2$TruncScoreTemp <- 1-(data2$AA_Pos/data2$Lenght)
data2$TruncScoreTemp[data2$TruncScoreTemp<0] <- 0
data2$TruncScore <- data2$ExonicFunc.refGene
data2$TruncScore[data2$TruncScore=="stopgain" & !is.na(data2$TruncScore)] <- (data2$TruncScoreTemp[data2$TruncScore=="stopgain" & !is.na(data2$TruncScore)])
data2$TruncScore[data2$TruncScore=="frameshift deletion" & !is.na(data2$TruncScore)] <- data2$TruncScoreTemp[data2$TruncScore=="frameshift deletion" & !is.na(data2$TruncScore)]
data2$TruncScore[data2$TruncScore=="frameshift insertion" & !is.na(data2$TruncScore)] <- data2$TruncScoreTemp[data2$TruncScore=="frameshift insertion" & !is.na(data2$TruncScore)]
data2$TruncScore <- gsub("nonframeshift insertion;nonframeshift insertion", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonframeshift deletion;nonframeshift deletion", 0, data2$TruncScore)
data2$TruncScore <- gsub("unknown;unknown", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonsynonymous SNV;nonsynonymous SNV", 0, data2$TruncScore)
data2$TruncScore <- gsub("stoploss;stoploss", 0, data2$TruncScore)
data2$TruncScore <- gsub("synonymous SNV;synonymous SNV", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonframeshift substitution;nonframeshift substitution", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonframeshift insertion", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonframeshift deletion", 0, data2$TruncScore)
data2$TruncScore <- gsub("unknown", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonsynonymous SNV", 0, data2$TruncScore)
data2$TruncScore <- gsub("stoploss", 0, data2$TruncScore)
data2$TruncScore <- gsub("synonymous SNV", 0, data2$TruncScore)
data2$TruncScore <- gsub("nonframeshift substitution", 0, data2$TruncScore)
data2$TruncScore <- as.numeric(data2$TruncScore)
data2$TruncScore[is.na(data2$TruncScore)] <- 0
# Mean Score
data2$Score <- rowMeans(subset(data2, select=c(Polyphen2_HDIV_score, Polyphen2_HVAR_score, VEST3_score, fathmm.MKL_coding_score, SIFT_scoreInv, MutationAssessor_score_Converted,  GERP.._RS_Converted, phyloP7way_vertebrate_score, phyloP20way_mammalian_score, TruncScore, LRT_score, MutationTaster_score, FATHMM_score_Converted, PROVEAN_score_Converted, CADD_raw_Converted, DANN_score, MetaSVM_score_Converted, MetaLR_score, integrated_fitCons_score_Converted, integrated_confidence_value_Converted, phastCons7way_vertebrate, phastCons20way_mammalian, SiPhy_29way_logOdds_Converted, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE, dpsi_max_tissue_Converted)), na.rm=T)


data3 <- data2[,c("Chr", "Start", "End", "cytoBand", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "avsnp144", "Ref", "Alt", "AA1", "AA2", "AA_Pos", "Lenght", "TranscriptID", "Exon", "Interpro_domain", "X1000g2015aug_all", "ExAC_ALL", "esp6500siv2_all", "FreqCases", "FreqControls", "CHISQ", "P", "OR", "SE", "L95", "U95", "Score", "TruncScore","Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "VEST3_score", "MutationAssessor_score_Converted", "GERP.._RS_Converted", "phyloP7way_vertebrate_score", "phyloP20way_mammalian_score", "SIFT_scoreInv", "fathmm.MKL_coding_score", "LRT_score", "MutationTaster_score", "FATHMM_score_Converted", "PROVEAN_score_Converted", "CADD_raw_Converted", "DANN_score", "MetaSVM_score_Converted", "MetaLR_score", "integrated_fitCons_score_Converted", "integrated_confidence_value_Converted", "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds_Converted", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "dpsi_max_tissue_Converted", "X1000g2015aug_afr", "X1000g2015aug_eas", "X1000g2015aug_eur", "X1000g2015aug_amr", "X1000g2015aug_sas", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "esp6500siv2_aa", "esp6500siv2_ea", "nci60", "cosmic70", "CLINSIG", "CLNDBN", "CLNACC", "CLNDSDB", "CLNDSDBID")]

