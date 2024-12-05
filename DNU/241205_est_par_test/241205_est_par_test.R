library(data.table)

dt_exp <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")
dt_exp <- change_data_type_exp(dt_exp)

dt_gt <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")
dt_gt <- change_data_type_gt(dt_gt)

load("D:/RStudio_GitHub/KongohPlus/inst/extdata/stutter_file/241205_stutter_GF_default.RData")
# list_dt_average_sr

dt_kit <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)

dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)

kit_name <- "GlobalFiler"

warmup <- 2000
sampling <- 2000
chains <- 4
pu_remove <- TRUE
pu_dist <- 1
at_sex <- 100
data_est_sr <- NULL
n_parallel_chains = detectCores()

model_ae <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_ae.stan")
model_gamma <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_gamma.stan")

##

peak_use <- peak_list[index_prop_ae_use]
gt_use <- gt_list[index_prop_ae_use]
est_alpha <- dt_result_est_ae$Alpha


# cv = 0.5439125
