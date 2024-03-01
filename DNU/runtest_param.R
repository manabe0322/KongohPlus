path_pack <- "D:/RStudio_GitHub/KongohPlus"

dt_exp <- fread(paste0(path_pack, "/inst/extdata/examples/expData_kakeiken.csv"))
dt_exp <- change_data_type_exp(dt_exp)

dt_gt <- fread(paste0(path_pack, "/inst/extdata/examples/gtData_kakeiken.csv"))
dt_gt <- change_data_type_gt(dt_gt)

dt_seq <- fread(paste0(path_pack, "/inst/extdata/examples/Sequence-data_GF_Example.csv"))
dt_seq <- change_data_type_seq(dt_seq)

kit <- "GlobalFiler"
dt_kit <- fread(paste0(path_pack, "/inst/extdata/kit/", kit, ".csv"))
dt_kit_auto <- extract_kit_auto(dt_kit)
markers_par <- unique(dt_kit_auto[, Marker])
dyes_par <- extract_dye(dt_kit_auto, markers_par)

dt_kit_set <- data.table(Marker = markers_par, Dye = dyes_par)
dt_kit_set$AT <- c(rep(12, 5), rep(20, 3), rep(13, 4), rep(16, 5), rep(19, 4))
dt_kit_set$Model_BSR <- "Best"
dt_kit_set$Model_FSR <- c(rep("Multiple loci together", 12), "Best", rep("Multiple loci together", 8))
dt_kit_set$Model_DSR <- "Multiple loci together"
dt_kit_set$Model_M2SR <- c(rep("Not consider", 16), "Uniform", "Not consider", "Uniform", "Not consider", "Not consider")

n_warmup <- 2000
n_sampling <- 2000
n_chains <- 4
pu_auto <- TRUE
pu_dist <- 1
min_n_motif_lus <- 2
min_n_motif_multi <- 3

error_message <- check_est_par(dt_exp, dt_gt, dt_seq, dt_kit_set)

model_ae <- readRDS(paste0(path_pack, "/inst/model/model_ae.rds"))
model_sr_gen <- readRDS(paste0(path_pack, "/inst/model/model_sr_gen.rds"))
model_sr_together <- readRDS(paste0(path_pack, "/inst/model/model_sr_together.rds"))
model_sr_uni <- readRDS(paste0(path_pack, "/inst/model/model_sr_uni.rds"))
model_gamma <- readRDS(paste0(path_pack, "/inst/model/model_gamma.rds"))
