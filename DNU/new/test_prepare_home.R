dt_csp <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/csp_example.csv")
dt_csp <- change_data_type_csp(dt_csp)

dt_af <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/Allele-count_GF_Japanese.csv")
dt_af <- change_data_type_af(dt_af)

dt_kit <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)
dt_kit_auto <- extract_kit_auto(dt_kit)
markers_auto <- unique(dt_kit_auto[, Marker])

dt_anameth <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/anameth/GF_default.csv")

load("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/param/param_GF_default.Rdata")
load("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/stutter_file/stutter_GF_default.Rdata")

dt_ignore_mk <- data.table(Ignore = rep(FALSE, length(markers_auto)),
                           Marker = markers_auto)

alpha_auto <- median_alpha
cv <- median_cv
hnc_from <- 2
hnc_to <- 2
hnc_all <- hnc_from:hnc_to

model_main <- cmdstan_model("C:/Users/manab/Documents/RStudio/KongohPlus/src/stan/model_main.stan")
