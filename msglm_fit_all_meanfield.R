# Msglm fit chunk for the Nadine Hinkel PCP hits AP-MS data

Sys.setenv(MKL_NUM_THREADS = 1)

project_id <- "sdenolly_canxapms"
msfolder <- "mq20211122"
fit_version <- "20211126"

message("Project ID=", project_id,
        " (fit_version=", fit_version, ") running on ", Sys.info()["nodename"], ")")

source('~/R/config.R')
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id, "_msglm_data_", fit_version, ".RData"))
message("Loading data from ", rdata_filepath)
load(rdata_filepath)

mcmc_nchains <- 8

require(rlang)
require(dplyr)
require(msglm)
require(stringr)
require(furrr)

plan(multicore, workers = 32)

try_fit_model <- function(stan_data, ntries, ...) {
  Reduce(function(last, i){
    if (!is.null(last) && !inherits_any(last, "try-error")) {
      last
    } else {
      message("Try #", i, " to fit the model")
      try({
        res <- fit_model(stan_data, ...)
        # check if res actually contains the fit, return error if not
        if(inherits_any(try(res$metadata()), "try-error")) res$metadata() else res
      })
    }
  }, seq_len(ntries), init=NULL)
}

opts <- furrr_options(packages = c("dplyr", "msglm", "rlang"),
              globals = c("msglm_def", "msdata", "try_fit",
                          "project_id", "fit_version",
                          "data_info", "mcmc_nchains"),
              seed = TRUE, stdout = FALSE, scheduling = 1)#Inf)

fit_chunks <- rowwise(msdata$objects) %>% group_map(function(r, ...) r) %>%
furrr::future_map(.progress = TRUE, .options = opts,
#purrr::map(
function(selobj_df) {
  message("#", selobj_df$chunk, " ", msdata$msentities[['object']],
          " ID=", selobj_df$object_id, ": ", selobj_df$object_label)

  model_data <- msglm_data(msglm_def, msdata, selobj_df$object_id,
                    specificity_msexp_group_cols = "bait",
                    cooccurrence_msexp_group_cols = "bait",
                    max_quantobjects = 80L,
                    quantobject_fdr = 0.01)
  stan_fit <- try_fit_model(model_data, 10L,
                  method="variational", algorithm="meanfield",
                  stanmodel_options=list(
                           object_labu_min_scale = 3,
                           batch_effect_sigma = 5.0,
                           object_msprobe_shift_tau = 0.05, object_msprobe_shift_df = 2,
                           effect_slab_scale = 4.0, effect_slab_df = 6,
                           quant_batch_tau = 1.0, quant_batch_df = 2,
                           quantobject_shift_sigma = 1.0, quantobject_shift_df = 4),
                  iter=2000L, chains=mcmc_nchains,
                  advi_iter=50000L, tol_rel_obj=0.001)
  if (inherits_any(stan_fit, "try-error")) {
    warning("Fitting chunk #", selobj_df$chunk,
            " (ID=", selobj_df$object_id, ") failed: ",
            attr(stan_fit, "condition"))
    msglm_results <- NULL
  } else if (inherits_any(stan_fit, "CmdStanFit")) {
    msglm_results <- process.stan_fit(stan_fit, model_data, contrast_method="normal")
  } else {
    warning("Unexpected stan_fit result")
    msglm_results <- NULL
  }
  gc()
  list(msglm_results = msglm_results)
})
names(fit_chunks) <- purrr::map_chr(fit_chunks, ~paste0(fit_version, '_',
                                                        .$msglm_results$objects$stats$object_id[1]))
fit_stats <- combine_fit_chunks(fit_chunks, 'stats')
fit_contrasts <- combine_fit_chunks(fit_chunks, 'contrast_stats')

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_meanfield_', fit_version, '.RData'))
results_info <- list(project_id = project_id, fit_version = fit_version)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts, file = rfit_filepath)
message('Done.')
