# Msglm fit chunk for the Solene Denolly CANX AP-MS

Sys.setenv(MKL_NUM_THREADS = 1)

#job.args <- c("sdenolly_canxapms", "sdenolly_canxapms_msglm", "20211125", "17")
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message("Project ID=", project_id)
job_name <- as.character(job.args[[2]])
job_version <- job.args[[3]]
data_version <- job.args[[3]]
fit_version <- job_version
job_chunk <- as.integer(job.args[[4]])

message('Job ', job_name, '(id=',job_chunk,
        " data_version=", data_version, " fit_version=", fit_version,
        " running on ", Sys.info()["nodename"], ")")

source('~/R/config.R')
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id, "_msglm_data_", data_version, ".RData"))
message("Loading data from ", rdata_filepath)
load(rdata_filepath)

if (Sys.getenv('SLURM_CPUS_PER_TASK') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) # SLURM way
} else if (Sys.getenv('NSLOTS') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('NSLOTS')) # SGE way
} else {
  mcmc_nchains <- 8
}

require(rlang)
require(dplyr)
require(msglm)
require(tidyr)
require(stringr)

sel_object_ids <- msdata$objects$object_id[[job_chunk]]
message(paste0(sel_object_ids, collapse=" "), " ", msdata$msentities[['object']], " ID(s): ",
        paste0(sort(unique(dplyr::filter(msdata$objects, object_id %in% sel_object_ids) %>% dplyr::pull(object_label))), collapse=' '))
model_data <- msglm_data(msglm_def, msdata, sel_object_ids,
                         max_quantobjects = 50L,
                         eager_msprotocols = FALSE,
                         specificity_msexp_group_cols = "bait",
                         cooccurrence_msexp_group_cols = "bait",
                         quantobject_fdr = 0.01)

# remove unneeded data to free some memory
msdata <- NULL
gc()

message('Running STAN in NUTS mode...')
options(mc.cores=mcmc_nchains)

msglm.stan_fit <- fit_model(model_data, adapt_delta=0.9,
                            iter=4000L, chains=mcmc_nchains,
                            stanmodel_options = list(
                               object_labu_min_scale = 3,
                               batch_effect_sigma = 1.0,
                               object_msprobe_shift_tau = 0.05, object_msprobe_shift_df = 2,
                               effect_slab_scale = 4.0, effect_slab_df = 6,
                               quant_batch_tau = 1.0, quant_batch_df = 2,
                               quantobject_shift_sigma = 1.0, quantobject_shift_df = 4
                            ))

msglm_results <- process.stan_fit(msglm.stan_fit, model_data)

res_prefix <- paste0(project_id, "_msglm")
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rfit_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', fit_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rfit_filepath, '...')
results_info <- list(project_id = project_id, msfolder = data_info$msfolder,
                     data_version = data_version, fit_version = fit_version,
                     job_name = job_name, job_chunk = job_chunk)
save(data_info, results_info, model_data, msglm_results,
     file = rfit_filepath)
message('Done.')
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)
