# Load data
arguments = commandArgs(trailingOnly = TRUE)
args <- commandArgs(TRUE)
data0 <- read.table(file=arguments[1], header=TRUE, stringsAsFactors=FALSE, na.strings=".", sep="\t", row.names=NULL, quote=NULL, fill=T)
HomCorr <- read.table(file=arguments[2], header=TRUE)
MAFFilter <- as.numeric(arguments[3])
Transcripts <- read.table(file=arguments[4], header=T, fill=T, na.strings=".", sep="\t")
TranscriptSize <- read.table("/data/Resources/Software/ceberSUITE/BASS/All_Transcripts_Size.txt", header=T)
ProteinSize <- read.table("/data/Resources/Software/ceberSUITE/BASS/All_Proteins_Size.txt", header=T)
XAcorrect <- read.table(file=arguments[5], header=TRUE)
# Prepare Data
data0$SIFT_scoreInv <- 1-data0$SIFT_score
data0$LRT_scoreInv <- 1-data0$LRT_score
data0$MutationAssessor_score_Converted <- ((data0$MutationAssessor_score+5.135)/(6.49+5.135))
data0$GERP.._RS_Converted <- (data0$GERP.._RS+12.3)/(6.17+12.3)
data0$phyloP20way_mammalian_score <- (data0$phyloP20way_mammalian+13.282)/(1.199+13.282)
data0$phyloP100way_vertebrate_Converted <- (data0$phyloP100way_vertebrate+20)/(10.003+20)
data0$FATHMM_score_Converted <- (data0$FATHMM_score+16.13)/(10.64+16.13)
data0$PROVEAN_score_Converted <- (data0$PROVEAN_score+14)/(14+14)
data0$CADD_raw_Converted <- (data0$CADD_raw+7.535)/(35.789+7.535)
data0$MetaSVM_score_Converted <- (data0$MetaSVM_score+2.006)/(3.04+2.006)
data0$integrated_fitCons_score_Converted <- (data0$integrated_fitCons_score+0)/(0.84+0)
data0$integrated_confidence_value_Converted <- (data0$integrated_confidence_value+0)/(3+0)
data0$SiPhy_29way_logOdds_Converted <- (data0$SiPhy_29way_logOdds-0.0003)/(37.972-0.0003)
data0$Eigen.raw_Converted <- (data0$Eigen.raw+4.154)/(4.065+4.154)
# Merge Data
data0$Chr_Start_End_Ref_Alt <- paste(data0$Chr, data0$Start, data0$End, data0$Ref, data0$Alt, sep="_")
data1 <- merge(data0, Transcripts, by.x="Chr_Start_End_Ref_Alt", by.y="Chr_Start_End_Ref_Alt", all.x=T)
data1$Pos <- paste(data1$Chr, data1$Start, sep="_")
HomCorr$Pos <- paste(HomCorr$ChrHom, HomCorr$BPHom, sep="_")
data2 <- merge(data1, HomCorr, by.x="Pos", by.y="Pos", all.x=TRUE)
XAcorrect$XAcorrect <- (XAcorrect$XAPA / XAcorrect$TOTAL)
XAcorrect$XAcorrect[is.na(XAcorrect$XAcorrect)] <- 0
data3 <- merge(data2, XAcorrect, by.x="Pos" , by.y="POS", all.x=TRUE)
data4temp <- merge(data3, TranscriptSize, by.x="TranscriptID", by.y="RefSeq", all.x=TRUE)
data4 <- merge(data4temp, ProteinSize, by.x="TranscriptID", by.y="RefSeq", all.x=TRUE)
# Truncation score
data4$AA_Pos <- as.numeric(as.character(data4$AA_Pos))
data4$TruncScoreTemp <- (1-(data4$AA_Pos/data4$ProteinSize))/2 + 0.5
data4$TruncScoreTemp[data4$TruncScoreTemp<0] <- 0
data4$TruncScore <- data4$ExonicFunc.refGene
data4$TruncScore[data4$TruncScore=="stopgain" & !is.na(data4$TruncScore)] <- (data4$TruncScoreTemp[data4$TruncScore=="stopgain" & !is.na(data4$TruncScore)])
data4$TruncScore[data4$TruncScore=="frameshift deletion" & !is.na(data4$TruncScore)] <- data4$TruncScoreTemp[data4$TruncScore=="frameshift deletion" & !is.na(data4$TruncScore)]
data4$TruncScore[data4$TruncScore=="frameshift insertion" & !is.na(data4$TruncScore)] <- data4$TruncScoreTemp[data4$TruncScore=="frameshift insertion" & !is.na(data4$TruncScore)]
data4$TruncScore <- gsub("nonframeshift insertion;nonframeshift insertion", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonframeshift deletion;nonframeshift deletion", NA, data4$TruncScore)
data4$TruncScore <- gsub("unknown;unknown", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonsynonymous SNV;nonsynonymous SNV", NA, data4$TruncScore)
data4$TruncScore <- gsub("stoploss;stoploss", NA, data4$TruncScore)
data4$TruncScore <- gsub("synonymous SNV;synonymous SNV", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonframeshift substitution;nonframeshift substitution", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonframeshift insertion", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonframeshift deletion", NA, data4$TruncScore)
data4$TruncScore <- gsub("unknown", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonsynonymous SNV", NA, data4$TruncScore)
data4$TruncScore <- gsub("stoploss", NA, data4$TruncScore)
data4$TruncScore <- gsub("synonymous SNV", NA, data4$TruncScore)
data4$TruncScore <- gsub("nonframeshift substitution", NA, data4$TruncScore)
data4$TruncScore <- as.numeric(data4$TruncScore)
# Runs of Homozygosity score
data4$HomozygScore <- data4$CorrHom/26000
data4$HomozygScore[is.na(data4$HomozygScore)] <- 0
# BASS Score
data4$MutPred_score <- as.numeric(data4$MutPred_score)
data4$ScoreTemp1 <- rowMeans(subset(data4, select=c(TruncScore, Polyphen2_HDIV_score, Polyphen2_HVAR_score, VEST3_score, fathmm.MKL_coding_score, SIFT_scoreInv, MutationAssessor_score_Converted,  GERP.._RS_Converted, phyloP100way_vertebrate_Converted, phyloP20way_mammalian_score, LRT_scoreInv, MutationTaster_score, FATHMM_score_Converted, PROVEAN_score_Converted, CADD_raw_Converted, DANN_score, MetaSVM_score_Converted, MetaLR_score, integrated_fitCons_score_Converted, integrated_confidence_value_Converted, phastCons100way_vertebrate, phastCons20way_mammalian, SiPhy_29way_logOdds_Converted, Eigen.raw_Converted, M.CAP_score, REVEL_score, MutPred_score, GenoCanyon_score, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE)), na.rm=T)
data4$Nscores <- rowSums(!is.na(subset(data4, select=c(TruncScore, Polyphen2_HDIV_score, Polyphen2_HVAR_score, VEST3_score, fathmm.MKL_coding_score, SIFT_scoreInv, MutationAssessor_score_Converted,  GERP.._RS_Converted, phyloP100way_vertebrate_Converted, phyloP20way_mammalian_score, LRT_scoreInv, MutationTaster_score, FATHMM_score_Converted, PROVEAN_score_Converted, CADD_raw_Converted, DANN_score, MetaSVM_score_Converted, MetaLR_score, integrated_fitCons_score_Converted, integrated_confidence_value_Converted, phastCons100way_vertebrate, phastCons20way_mammalian, SiPhy_29way_logOdds_Converted, Eigen.raw_Converted, M.CAP_score, REVEL_score, MutPred_score, GenoCanyon_score, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE))))
data4$ScoreTemp2 <- (data4$ScoreTemp1 - ((30 - data4$Nscores)/30)/10)
data4$ScoreTemp2[is.na(data4$ScoreTemp1)] <- 0
data4$XAcorrect[is.na(data4$XAcorrect)] <- 0
data4$Score <- data4$ScoreTemp2 - data4$XAcorrect + data4$HomozygScore
data4$Score[is.na(data4$Score)] <- 0
data4$Score[ data4$Score < 0 ] <- 0
data4$Score[ data4$Score > 1 ] <- 1
# Filter data for set MAF threshold
data4MAF <- subset(data4, data4$X1000g2015aug_all < MAFFilter | is.na(data4$X1000g2015aug_all))
data4MAFpi <- subset(data4MAF, data4MAF$gnomAD_genome_ALL < MAFFilter | is.na(data4MAF$gnomAD_genome_ALL))
data4MAF2 <- subset(data4MAFpi, data4MAFpi$ExAC_nontcga_ALL < MAFFilter | is.na(data4MAFpi$ExAC_nontcga_ALL))
# Replace the names of the alleles for somthing nicer?
data4MAF2$GT <- gsub(12, "Alt12", data4MAF2$GT)
data4MAF2$GT <- gsub(11, "Alt11", data4MAF2$GT)
data4MAF2$GT <- gsub(10, "Alt10", data4MAF2$GT)
data4MAF2$GT <- gsub(9, "Alt9", data4MAF2$GT)
data4MAF2$GT <- gsub(8, "Alt8", data4MAF2$GT)
data4MAF2$GT <- gsub(7, "Alt7", data4MAF2$GT)
data4MAF2$GT <- gsub(6, "Alt6", data4MAF2$GT)
data4MAF2$GT <- gsub(5, "Alt5", data4MAF2$GT)
data4MAF2$GT <- gsub(4, "Alt4", data4MAF2$GT)
data4MAF2$GT <- gsub(3, "Alt3", data4MAF2$GT)
data4MAF2$GT <- gsub(2, "Alt2", data4MAF2$GT)
data4MAF2$GT <- gsub(1, "Alt1", data4MAF2$GT)
data4MAF2$GT <- gsub(0, "Ref", data4MAF2$GT)
# Calculate VAF for each variant
if("AD" %in% colnames(data4MAF2)) {
    if("DP" %in% colnames(data4MAF2)) {
        data4MAF2$ADlow <- sapply(strsplit(data4MAF2$AD, ","), min)
        data4MAF2$VAF <- as.numeric(data4MAF2$ADlow) / as.numeric(data4MAF2$DP)
    }
}


if("DV" %in% colnames(data4MAF2)) {
    if("DP" %in% colnames(data4MAF2)) {
        data4MAF2$VAF <- data4MAF2$DV/data4MAF2$DP
    }
}
# Reorder and select columns
if("AD" %in% colnames(data4MAF2)) {
    if("DP" %in% colnames(data4MAF2)) {
        if("GT" %in% colnames(data4MAF2)) {
            if("GQ" %in% colnames(data4MAF2)) {
                if("PL" %in% colnames(data4MAF2)) {
                    data5 <- data4MAF2[,c("Chr", "Start", "End", "cytoBand", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "avsnp150", "HGVS_cDNA", "Ref", "Alt", "CDS_Pos", "TranscriptSize", "HGVS_Prot", "AA1", "AA2", "AA_Pos", "ProteinSize", "TranscriptID", "Exon", "Interpro_domain", "X1000g2015aug_all", "ExAC_ALL", "ExAC_nontcga_ALL", "gnomAD_genome_ALL", "esp6500siv2_all", "GT", "AD", "VAF", "Depth", "Score", "CLNSIG", "InterVar_automated", "cosmic87", "Oncogenicity", "civic_score", "TruncScore", "HomozygScore", "XAcorrect", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "VEST3_score", "fathmm.MKL_coding_score", "SIFT_scoreInv", "MutationAssessor_score_Converted", "GERP.._RS_Converted", "phyloP100way_vertebrate_Converted", "phyloP20way_mammalian_score", "LRT_scoreInv", "MutationTaster_score", "FATHMM_score_Converted", "PROVEAN_score_Converted", "CADD_raw_Converted", "DANN_score", "MetaSVM_score_Converted", "MetaLR_score", "integrated_confidence_value_Converted", "integrated_fitCons_score_Converted", "phastCons100way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds_Converted", "Eigen.raw_Converted", "M.CAP_score", "REVEL_score", "MutPred_score" ,"GenoCanyon_score", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "X1000g2015aug_afr", "X1000g2015aug_eas", "X1000g2015aug_eur", "X1000g2015aug_amr", "X1000g2015aug_sas", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "esp6500siv2_aa", "esp6500siv2_ea", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ", "gnomAD_genome_EAS", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH", "nci60", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNALLELEID", "Mutation_CancerEffect", "PMIDs", "Gene_Cancer_Role", "civic_variant_id", "Q", "GQ", "PL")]
                }
            }
        }
    }
}
if("GT" %in% colnames(data4MAF2)) {
    if("DV" %in% colnames(data4MAF2)) {
        if("DP" %in% colnames(data4MAF2)) {
            data5 <- data4MAF2[,c("Chr", "Start", "End", "cytoBand", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "avsnp150", "HGVS_cDNA", "Ref", "Alt", "CDS_Pos", "TranscriptSize", "HGVS_Prot", "AA1", "AA2", "AA_Pos", "ProteinSize", "TranscriptID", "Exon", "Interpro_domain", "X1000g2015aug_all", "ExAC_ALL", "ExAC_nontcga_ALL", "gnomAD_genome_ALL", "esp6500siv2_all", "GT", "DV", "VAF", "Depth", "Score", "CLINSIG", "InterVar_automated", "cosmic87", "Oncogenicity", "civic_score", "TruncScore", "HomozygScore", "XAcorrect", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "VEST3_score", "fathmm.MKL_coding_score", "SIFT_scoreInv", "MutationAssessor_score_Converted", "GERP.._RS_Converted", "phyloP100way_vertebrate_Converted", "phyloP20way_mammalian_score", "LRT_scoreInv", "MutationTaster_score", "FATHMM_score_Converted", "PROVEAN_score_Converted", "CADD_raw_Converted", "DANN_score", "MetaSVM_score_Converted", "MetaLR_score", "integrated_confidence_value_Converted", "integrated_fitCons_score_Converted", "phastCons100way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds_Converted", "Eigen.raw_Converted", "M.CAP_score", "REVEL_score", "MutPred_score" ,"GenoCanyon_score", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "X1000g2015aug_afr", "X1000g2015aug_eas", "X1000g2015aug_eur", "X1000g2015aug_amr", "X1000g2015aug_sas", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "esp6500siv2_aa", "esp6500siv2_ea", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ", "gnomAD_genome_EAS", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH", "nci60", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNALLELEID", "Mutation_CancerEffect", "PMIDs", "Gene_Cancer_Role", "civic_variant_id", "Q")]
        }
    }
}
# Order file for Score
data6 <- data5[order(-data5$Score),]
# Export data
write.table(data6, file="Scores.txt", sep="\t", quote=F, row.names=F, na="-9")
