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
require(stringr)
require(furrr)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', fit_version, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_msdata_full_', data_version, '.RData')))

message('Loading MSGLM model fit results...')

modelobj <- msdata$msentities['modelobject']
quantobj <- msdata$msentities['quantobject']
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

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
