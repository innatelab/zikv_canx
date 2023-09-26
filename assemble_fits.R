# assembly and analysis of Solene Denolly CANX AP-MS MSGLM analysis
#
# Author: Alexey Stukalov
###############################################################################

project_id <- "sdenolly_canxapms"
message('Project ID=', project_id)
data_version <- "20211125"
fit_version <- "20211126"
message("Assembling fit results for project ", project_id,
        " (dataset v", data_version, ", fit v", fit_version, ")")

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

require(msglm)
require(dplyr)
require(tidyr)
require(stringr)
require(furrr)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', fit_version, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_msdata_full_', data_version, '.RData')))

message('Loading MSGLM model fit results...')

fit_path <- file.path(scratch_path, paste0(project_id, '_msglm'))#_', fit_version))
fit_files <- list.files(fit_path, paste0(project_id, '_msglm_', fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
    tidyr::extract(filename, "chunk", ".+_(\\d+).RData$", convert=TRUE, remove=FALSE) %>%
    dplyr::left_join(dplyr::select(msdata$objects, chunk, object_id), by="chunk") %>%
    dplyr::arrange(chunk)
require(RMySQL)
chunk_dispatcher_conn <- dbConnect(RMySQL::MySQL(),
                                    dbname="inlab_computing", user="inlab_dispatcher",
                                    host="tumevi4-websrv1.srv.mwn.de",
                                    password=Sys.getenv("INLAB_DISPATCHER_PASSWD"),
                                    timeout=300)
chunk_statuses.df <- dbGetQuery(chunk_dispatcher_conn,
                                str_c("SELECT * FROM jobchunks WHERE user='ge68wan2' AND job_id='",
                                      project_id, "_", fit_version, "'")) %>%
  dplyr::mutate(start_time = as.POSIXlt(start_time),
                end_time = as.POSIXlt(end_time),
                fit_time = end_time - start_time)
dbDisconnect(chunk_dispatcher_conn)
table(chunk_statuses.df$status)
fit_files.df <- dplyr::inner_join(fit_files.df,
                                  dplyr::filter(chunk_statuses.df, status == "complete" & fit_time > 0) %>%
                    dplyr::select(chunk, status, fit_time), by="chunk")
# load fit results in parallel
plan(multisession, workers = 32)
fit_chunks <- seq_len(nrow(fit_files.df)) %>%
  furrr::future_map(.progress = TRUE, .options = furrr_options(stdout=FALSE, packages=c("msglm"), globals=c("fit_path", "fit_files.df")),
                    ~load_fit_chunk(.x))
names(fit_chunks) <- purrr::map_chr(fit_chunks, ~paste0(.$results_info$fit_version, '_',
                                                        .$msglm_results$objects$stats$object_id[1]))

fit_stats <- combine_fit_chunks(fit_chunks, 'stats')
fit_contrasts <- combine_fit_chunks(fit_chunks, 'contrast_stats')

rm(fit_chunks)

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_', fit_version, '.RData'))
results_info <- list(project_id = project_id, data_version = data_version,
                     fit_version = fit_version)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts, file = rfit_filepath)
message('Done.')

contrasts.df <- dplyr::ungroup(msglm_def$contrasts) %>%
  dplyr::inner_join(tidyr::pivot_wider(dplyr::mutate(as.data.frame.table(msglm_def$metaconditionXcontrast, responseName="w"),
                                   side = if_else(w > 0, "lhs", "rhs")) %>% dplyr::filter(w != 0),
                        c(contrast), names_from = "side", values_from = "metacondition",
                        names_glue = "{.value}_{side}"), by="contrast") %>%
  dplyr::mutate(offset = 0, offset_prior = 0) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, bait_lhs = bait, treatment_lhs = treatment)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, bait_rhs = bait, treatment_rhs = treatment))

object_contrasts_thresholds.df <- dplyr::select(contrasts.df, offset, offset_prior, contrast, contrast_type) %>%
  dplyr::mutate(
    p_value_threshold = case_when(contrast_type=="filter" ~ 1E-3,
                                  contrast_type=="comparison" ~ 1E-3,
                                  TRUE ~ NA_real_),
    median_threshold = case_when(contrast_type=="filter" ~ pmax(2.0, 2.0 + abs(offset - offset_prior)),
                                 contrast_type=="comparison" ~ pmax(1.0, 0.25 + abs(offset - offset_prior)),
                                 TRUE ~ NA_real_))

object_contrasts.df <- fit_contrasts$object_conditions %>%
  dplyr::filter(var %in% c('obj_cond_labu')) %>%
  dplyr::inner_join(dplyr::select(object_contrasts_thresholds.df, contrast, p_value_threshold, median_threshold)) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_reverse & !is_contaminant,
                is_hit = is_hit_nomschecks,
                change = if_else(is_signif, if_else(median > 0, "+", "-"), NA_character_)) %>%
  dplyr::select(object_id, contrast, mean, median, sd, rhat, p_value, change)

object_effects_thresholds.df <- dplyr::select(msglm_def$effects, effect, prior_mean) %>%
  dplyr::mutate(p_value_threshold = 1E-3,
                median_threshold = c(0.25)
  )

object_effects.df <- fit_stats$object_effects %>%
  dplyr::filter(var %in% c('obj_effect')) %>%
  dplyr::inner_join(dplyr::select(object_effects_thresholds.df, effect, p_value_threshold, median_threshold)) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - prior_mean) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_reverse & !is_contaminant,
                is_hit = is_hit_nomschecks,
                change = if_else(is_signif, if_else(median > 0, "+", "-"), NA_character_)) %>%
  dplyr::select(object_id, bait, treatment, effect_label, mean, median, sd, rhat, p_value, change)

object_bait_treatment_effects.df <-
  dplyr::filter(object_effects.df, !is.na(bait) & !is.na(treatment)) %>%
  pivot_wider(id_cols = "object_id", names_from="effect_label", names_sep=".",
              values_from = c("median", "mean", "sd", "p_value", "change"))

object_contrasts_report.df <- dplyr::select(msdata$objects,
    object_id, gene_names, majority_protein_acs, protein_description, npepmods, npepmods_unique, starts_with("is_")) %>%
  dplyr::left_join(
  pivot_wider_spec(object_contrasts.df,
  build_wider_spec(object_contrasts.df,
      names_from = "contrast", names_sep = ".", names_sort=TRUE,
      values_from = c("median", "mean", "sd", "p_value", "change")) %>%
      dplyr::arrange(desc(.value == "change"), contrast),
    id_cols = c("object_id")), by="object_id") %>%
  dplyr::left_join(object_bait_treatment_effects.df, by="object_id") %>%
  dplyr::relocate(`change.CANX+ZikV`, .before=`change.CANX_vs_empty@Mock`) %>%
  dplyr::mutate(is_hit = !is.na(`change.CANX+ZikV`) &
                          (coalesce(`change.CANX_vs_empty@Mock`, "") == "+" |
                           coalesce(`change.CANX_vs_empty@ZikV`, "") == "+") &
                         !is_reverse & !is_contaminant) %>%
  dplyr::relocate(is_hit, .before=`change.CANX+ZikV`)

writexl::write_xlsx(object_contrasts_report.df,
                    file.path(analysis_path, "results", paste0(project_id, "_msglm_report_", fit_version, "_3.xlsx")))
