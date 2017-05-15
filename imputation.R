library('mice')
# Path
bl_data <- read.csv('/data/processed/baseline_clinical_data.csv')

bl_data_2 <- cbind(bl_data)

public_id <- bl_data_2[, c("PUBLIC_ID")]

bl_data_2$PUBLIC_ID <- NULL

to_remove <- c("CMMC_REASONFORPROC",
               "CMMC_REASONCODE",
               "CMMC_RECDY",
               "CMMC",
               "D_CM_ANEUPLOIDYCAT",
               "D_IM_CD38_PC_PERCENT",
               "D_IM_CD138_PC_PERCENT",
               "D_IM_CD45_PC_PERCENT",
               "D_IM_CD56_PC_PERCENT",
               "D_IM_CD117_PC_PERCENT",
               "D_IM_CD138_DETECTED",
               "D_IM_CD38_DETECTED",
               "D_IM_CD45_PC_TYPICAL_DETECTED",
               "D_IM_CD56_DETECTED",
               "D_IM_CD13_DETECTED",
               "D_IM_CD20_DETECTED",
               "D_IM_CD33_DETECTED",
               "D_IM_CD52_DETECTED",
               "D_IM_CD117_DETECTED",
               "DEMOG_NATIVEHAWAIIA",
               "D_CM_T614",
               "D_IM_fgfr3_detected",
               "D_LAB_chem_ldh",
               "D_LAB_urine_24hr_total_protein",
               "D_LAB_urine_24hr_m_protein",
               "SS_DOESTHEPATIEN")

far_better_bl <- bl_data_2[, !names(bl_data_2) %in% to_remove]

imp <- mice(far_better_bl, m=8, defaultMethod = c("pmm", "logreg", "polyreg", "polr"))

imputed <- complete(imp, 3)

imputed_o <- cbind("PUBLIC_ID" = public_id, imputed)

write.csv(imputed_o, "imputed_h_and_w_update.csv")