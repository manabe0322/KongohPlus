#dt_csp <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/mix15_CSP.csv")
dt_csp <- fread("D:/研究プロジェクト管理/research_note_KongohPlus/240816_KMU_validation/other_data/KMU_csp/KMU_mix106_csp.csv")
dt_csp <- change_data_type_csp(dt_csp)

dt_af <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/Allele-count_GF_Japanese.csv")
dt_af <- change_data_type_af(dt_af)

dt_kit <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)
dt_kit_auto <- extract_kit_auto(dt_kit)
markers_auto <- unique(dt_kit_auto[, Marker])

dt_anameth <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/anameth/240816_KMU_GF30.csv")

# Create a data.table for information on markers
dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                    AT = as.numeric(dt_anameth$Value[dt_anameth$Category == "AT"]))
dt_info_mk <- copy(dt_at)
dt_info_mk <- cbind(
  Ignore = rep(FALSE, nrow(dt_info_mk)),
  dt_info_mk
)

load("D:/RStudio_GitHub/KongohPlus/inst/extdata/param/240816_KMU_edit-cv.Rdata")
dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Allele_probability"],
                         Value = dt_anameth$Value[dt_anameth$Category == "Allele_probability"])
dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype_filter"],
                         Value = dt_anameth$Value[dt_anameth$Category == "Genotype_filter"])
dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                        Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

hnc <- 4
warmup <- 1200
sampling <- 300
chains <- 4

model_main <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_main.stan")
