dt_csp <- fread("D:/研究プロジェクト管理/research_note_KongohPlus/240816_KMU_validation/other_data/KMU_csp/KMU_mix14_2_csp.csv")
dt_csp <- change_data_type_csp(dt_csp)

dt_af <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/Allele-count_GF_Japanese.csv")
dt_af <- change_data_type_af(dt_af)

dt_kit <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)
dt_kit_auto <- extract_kit_auto(dt_kit)
markers_auto <- unique(dt_kit_auto[, Marker])

dt_anameth <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/anameth/GF_KMU30.csv")

load("D:/RStudio_GitHub/KongohPlus/inst/extdata/param/Param_GF_KMU30.Rdata")
load("D:/RStudio_GitHub/KongohPlus/inst/extdata/stutter/Stutter_GF_KMU30.Rdata")

dt_ignore_mk <- data.table(Ignore = rep(FALSE, length(markers_auto)),
                           Marker = markers_auto)

n_mk <- length(median_alpha)
expect_ae <- n_mk * (median_alpha + 1) / (sum(median_alpha) + n_mk)
cv_allele <- median_cv_allele
cv_stutter <- median_cv_stutter
hnc_from <- 2
hnc_to <- 2
hnc_all <- hnc_from:hnc_to

model_main <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_main.stan")
model_main_rep1 <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_main_rep1.stan")
