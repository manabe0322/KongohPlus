dt_exp <- fread("D:/研究プロジェクト管理/research_note_KongohPlus/250311_KMU30_estimate-parameter/KMU30_stutter-ratio_all-data/dt_exp_rm_dropin.csv")
dt_exp <- change_data_type_exp(dt_exp)

dt_gt <- fread("D:/研究プロジェクト管理/research_note_KongohPlus/250311_KMU30_estimate-parameter/KMU30_stutter-ratio_all-data/genotype-data_KMU30.csv")
dt_gt <- change_data_type_gt(dt_gt)

dt_seq <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/seq/sequence_data_japanese.csv")
dt_seq <- change_data_type_seq(dt_seq)

dt_kit <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)

kit_name <- "GlobalFiler"

warmup <- 2000
sampling <- 2000
chains <- 4
pu_remove <- TRUE
pu_dist <- 1
at_sex <- 40
data_est_sr <- NULL
n_parallel_chains = detectCores()

model_sr_gen <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_sr_gen.stan")
model_sr_uni <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_sr_uni.stan")

# 条件1: Allele model
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
dt_setting_est_sr[["BSR model"]] <- "Allele"
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- NA
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 1

# 条件2: LUS model
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
dt_setting_est_sr[["BSR model"]] <- "LUS"
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- NA
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 2

# 条件3: Multi-seq model, X = 3 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 3
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 3

# 条件4: Multi-seq model, X = 4 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 4
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 4

# 条件5: Multi-seq model, X = 5 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 5
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 5

# 条件6: Multi-seq model, X = 6 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 6
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 6

# 条件7: Multi-seq model, X = 7 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 7
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 7

# 条件8: Multi-seq model, X = 8 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 8
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 8

# 条件9: Multi-seq model, X = 9 (D21S11, D2S441, FGA, SE33, D1S1656, D12S391, D2S1338)
dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
bsr_model_name <- rep("Not consider", 21)
bsr_model_name[c(7, 9, 12, 17, 19, 20, 21)] <- "Multi-seq"
dt_setting_est_sr[["BSR model"]] <- bsr_model_name
dt_setting_est_sr[["FSR model"]] <- "Not consider"
dt_setting_est_sr[["DSR model"]] <- "Not consider"
dt_setting_est_sr[["M2SR model"]] <- "Not consider"
x_values <- rep(NA, 21)
x_values[c(7, 9, 12, 17, 19, 20, 21)] <- 9
dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]] <- x_values
col_int <- c("Minimum uninterrupted stretch (for Multi-seq)")
dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]
condition_no <- 9

#### 以下、MCMCシミュレーション

# Extract information on the kit
markers_all <- unique(dt_kit$Marker)
n_mk_all <- length(markers_all)
dye_all_mk <- extract_dye(dt_kit, markers_all)

# Extract information on setting
markers_auto <- dt_setting_est_sr$Marker
n_mk_auto <- length(markers_auto)
at_auto <- dt_setting_est_sr$AT
names(at_auto) <- markers_auto
meth_model_bsr <- dt_setting_est_sr[["BSR model"]]
meth_model_fsr <- dt_setting_est_sr[["FSR model"]]
meth_model_dsr <- dt_setting_est_sr[["DSR model"]]
meth_model_m2sr <- dt_setting_est_sr[["M2SR model"]]

# Make an object for ATs of all markers
at_all_mk <- rep(NA, n_mk_all)
names(at_all_mk) <- markers_all
at_all_mk[match(markers_auto, markers_all)] <- at_auto
at_all_mk[is.na(at_all_mk)] <- at_sex

# Information on considering stutter ratios
sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
rownames(sr_consider_auto) <- markers_auto
colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
sr_consider_auto[, 1] <- meth_model_bsr != "Not consider"
sr_consider_auto[, 2] <- meth_model_fsr != "Not consider"
sr_consider_auto[, 3] <- meth_model_dsr != "Not consider"
sr_consider_auto[, 4] <- meth_model_m2sr != "Not consider"

# Arrange experimental data
tmp <- arrange_exp_data(dt_exp, dt_gt,
                        markers_auto, sr_consider_auto,
                        markers_all, at_all_mk, dye_all_mk,
                        pu_remove, pu_dist)
peak_list <- tmp$peak_list
size_list <- tmp$size_list
height_list <- tmp$height_list
gt_list <- tmp$gt_list

# Prepare data related to the sequence data
min_us_auto <- dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]]
tmp <- correct_al_rep(dt_seq, markers_auto, min_us_auto)
wlus_auto <- tmp$wlus_auto
wca_auto <- tmp$wca_auto

# Create a data.table for the stutter ratios
dt_sr_data <- create_dt_sr_data(peak_list, height_list, gt_list, sr_consider_auto, wlus_auto)

# Perform MCMC for BSR
tmp <- perform_mcmc_sr_all_mk(dt_sr_data, "BS", markers_auto, meth_model_bsr, wlus_auto, wca_auto,
                              model_sr_gen, model_sr_uni,
                              warmup, sampling, chains, n_parallel_chains)
dt_result_est_bsr <- tmp$dt_result_est_sr
list_trace_bsr <- tmp$trace_all_mk
list_dt_average_bsr <- tmp$dt_average_sr_all_mk

write.csv(dt_result_est_bsr, paste0("Condition", condition_no, "_seq-japanese_dt_result_est_bsr.csv"), row.names = FALSE)
for(i in 1:n_mk_auto){
  write.csv(list_dt_average_bsr[[i]], paste0("Condition", condition_no, "_seq-japanese_dt_average_bsr_", i, "_", markers_auto[i], ".csv"), row.names = FALSE)
}
