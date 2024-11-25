library(data.table)

dt_exp <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")
dt_exp <- change_data_type_exp(dt_exp)

dt_gt <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")
dt_gt <- change_data_type_gt(dt_gt)

dt_seq <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/seq/sequence_data_default.csv")
dt_seq <- change_data_type_seq(dt_seq)

dt_kit <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/kit/GlobalFiler.csv")
dt_kit <- change_data_type_kit(dt_kit)

dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)

kit_name <- "GlobalFiler"

warmup <- 2000
sampling <- 2000
chains <- 4
pu_remove <- TRUE
pu_dist <- 1
at_sex <- 100
data_est_sr <- NULL
n_parallel_chains = detectCores()

model_sr_gen <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_sr_gen.stan")
model_sr_uni <- cmdstan_model("D:/RStudio_GitHub/KongohPlus/src/stan/model_sr_uni.stan")
