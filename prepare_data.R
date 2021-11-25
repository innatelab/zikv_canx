# Preparing data for the analysis of proteins relocation upon ZIKV infection (Solene Denolly)
#
# Author: Alexey Stukalov
###############################################################################

project_id <- "sdenolly_canxapms"
msfolder <- "mq20211122"
data_version <- "20211125"
message('Project ID=', project_id, " data version=", data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(msimportr)
require(msglm)
require(dplyr)
require(jsonlite)
require(stringr)
require(readr)
require(pheatmap)

data_info <- list(project_id = project_id, data_ver=data_version, fit_ver=data_version,
                  pepmodstate_mscalib_filename = "mscalib_QEP5_intensity_pepmodstate_cov2_20211108.json",
                  msfolder = msfolder)
msdata_path <- file.path(data_path, msfolder)

message('Loading MS instrument protein calibration data from ', data_info$pepmodstate_mscalib_filename, '...')
pepmodstate_mscalib <- read_mscalib_json(file.path(data_path, data_info$pepmodstate_mscalib_filename))# %>% convert_logintensityBase(new_base=2)

msruns.df <- read_tsv(file.path(msdata_path, "msruns.txt")) %>%
  tidyr::extract("msrun", c('bait', 'treatment', 'batch'),
                 "^(.+)_(.+)_(.+)$", remove = FALSE, convert=TRUE) %>%
  dplyr::mutate(msrun_ix = row_number(),
                rawfile = str_remove(rawfile, "\\.raw$"),
                bait = factor(bait) %>% relevel("empty"),
                treatment = factor(treatment) %>% relevel("Mock")) %>%
  dplyr::arrange(bait, treatment, batch) %>%
  dplyr::mutate(condition = paste0(bait, '_', treatment),
                condition = factor(condition, levels=unique(condition)) %>% relevel("empty_Mock"))

fasta.dfs <- list(
  HAtag = read_innate_uniprot_fasta(file.path(data_path, "fasta", "hatag.fasta")),
  ZikV = read_innate_uniprot_fasta(file.path(data_path, "fasta", "A0A024B7W1_ZIKVF.fasta")),
  human = bind_rows(
      read_innate_uniprot_fasta(file.path(data_path, "fasta", "UP000005640_9606.fasta")),
      read_innate_uniprot_fasta(file.path(data_path, "fasta", "UP000005640_9606_additional.fasta"))
  )
)

fix_viral_protein_acs <- function(acs) {
  str_remove_all(str_remove_all(acs, "(?<=^|;)(?:sp|tr)\\|"), "\\|[^|]+(?=;|$)")
}

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(msdata_path, 'combined/txt'), import_data = c(data_info$quant_type, "ident_type"))
# fix protein acs
msdata.wide <- dplyr::mutate(msdata.wide,
                             across(matches("_acs?$"), fix_viral_protein_acs))
msdata_colgroups <- attr(msdata.wide, "column_groups")

mqevidence <- read.MaxQuant.Evidence(file.path(msdata_path, 'combined', 'txt'),
                                     mschannel_annotate.f = function(mschans_df) {
                                       res <- dplyr::inner_join(dplyr::select(mschans_df, mstag, rawfile),
                                                                dplyr::select(msruns.df, rawfile, condition, treatment, bait,
                                                                              msrun=msrun, msexperiment=msrun, batch),
                                                                by="rawfile")
                                                                print(res)
                                       attr(res, "column_scopes") <- c(treatment = "msexperiment",
                                                                       bait = "msexperiment",
                                                                       condition = "msexperiment",
                                                                       batch = "msexperiment")
                                       return(res)
                                     })
mqevidence$peptides <- read.MaxQuant.Peptides(file.path(msdata_path, 'combined', 'txt'), file_name='peptides.txt',
                                              import_data='ident_type')
mqevidence$peaks <- NULL # exclude big data frame

for (dfname in c("pepmods", "pepmodstates", "peptides")) {
  mqevidence[[dfname]] <- dplyr::mutate(mqevidence[[dfname]], across(matches("_acs?$"), fix_viral_protein_acs))
}

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

# all ms data
msdata_full <- list(msexperiments = dplyr::select(mqevidence$msexperiments, -starts_with("msfraction")),
                    msruns = mqevidence$msruns)

msdata_full <- append_protgroups_info(msdata_full, dplyr::mutate(msdata.wide, organism=NULL),
                                      proteins_info = dplyr::bind_rows(fasta.dfs) %>%
                                        dplyr::mutate(protein_ac_noiso = str_remove(protein_ac, "-\\d+(?:#.+)*$"),
                                                      protein_isoform_ix = replace_na(as.integer(str_match(protein_ac, "-(\\d+)$")[, 2]), 1L)),
                                      import_columns = c("organism"))

msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(dplyr::bind_rows(fasta.dfs), lead_razor_protein_ac = protein_ac, organism)) %>%
  dplyr::mutate(peptide_rank = 1L)
msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::mutate(pepmod_rank = 1L)
msdata_full$pepmodstates <- mqevidence$pepmodstates

# redefine protein groups (protregroups)
peptides.df <- dplyr::select(msdata_full$peptides, peptide_id, protgroup_ids, protein_acs, lead_razor_protein_ac,
                             peptide_seq, is_reverse, peptide_rank)
proteins.df <- msdata_full$proteins
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, '_', data_version, "_peptides.RData")),
     peptides.df, proteins.df)
# .. run protregroup.jl
msdata_full$protregroups <- read_tsv(file.path(data_path, msfolder,
                                               str_c(project_id, "_", msfolder, '_', data_version, "_protregroups.txt")),
                                     col_types = list(protregroup_id = "i"))

msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
  separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
  dplyr::mutate(is_majority = TRUE) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::mutate(protein_ac_rank = row_number()) %>%
  dplyr::ungroup()

msdata_full$protregroup2peptide <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, peptide_id=spec_peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, peptide_id=peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, peptide_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()
msdata_full$protregroup2pepmod <- dplyr::inner_join(msdata_full$protregroup2peptide,
                                                    dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id)) %>%
  dplyr::select(-peptide_id)
msdata_full$protregroup2pepmodstate <- dplyr::semi_join(msdata_full$protregroup2pepmod,
                                                        dplyr::select(msdata_full$pepmods, pepmod_id)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmod_id, pepmodstate_id))

msdata_full$protregroups <- msdata_full$protregroups %>%
  dplyr::mutate(gene_label = strlist_label2(gene_names),
                protein_label = strlist_label2(protein_names),
                protein_description = strlist_label2(protein_descriptions),
                is_viral = replace_na(str_detect(organism, "ZIKV"), FALSE),
                protac_label = strlist_label2(majority_protein_acs),
                protregroup_label = case_when(is_viral ~ protein_label,
                                              !is.na(gene_label) ~ gene_label,
                                              !is.na(protac_label) ~ protac_label,
                                              TRUE ~ str_c('#', protregroup_id))) %>%
  dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                   dplyr::group_by(protregroup_id) %>%
                   dplyr::summarise(npeptides = n_distinct(peptide_id),
                                    npepmods = n_distinct(pepmod_id),
                                    npeptides_unique = n_distinct(peptide_id[is_specific]),
                                    npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                 npeptides_razor = 0L,
                                 npepmods_unique_razor = npepmods_unique,
                                 npepmods_razor = 0L))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  dplyr::mutate(is_idented = str_detect(ident_type, "MSMS"))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod,
                                                    msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(msexperiment, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# setup experimental design matrices
conditions.df <- bind_rows(dplyr::select(msdata_full$msexperiments, condition, treatment, bait) %>%
                           dplyr::distinct() %>%
                           dplyr::mutate(is_virtual = FALSE)) %>%
  dplyr::arrange(is_virtual, condition) %>%
  dplyr::mutate(is_negctrl = bait == 'empty' & treatment == 'Mock')

conditionXeffect.mtx <- model.matrix(~ 1 + treatment*bait, conditions.df)
conditionXeffect.mtx <- conditionXeffect.mtx[, setdiff(colnames(conditionXeffect.mtx), '(Intercept)')] # handled separately
dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))
all_conditions <- as.character(conditions.df$condition)
plot_path <- file.path(analysis_path, 'plots', str_c(msfolder, "_", data_info$fit_ver))
if (!dir.exists(plot_path)) dir.create(plot_path)
pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, na_col = "gray", border_color = "white",
         filename = file.path(plot_path,
                              paste0(project_id, "_conditionXeffect_", data_info$msfolder, "_", data_info$fit_ver, ".pdf")),
         width = 6, height = 5)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  dplyr::mutate(treatment = effect_factor(effect, "treatment", levels(conditions.df$treatment), NA_character_),
                bait = effect_factor(effect, "bait", levels(conditions.df$bait), NA_character_),
                is_positive = FALSE) %>%
  dplyr::mutate(effect_type = case_when(
                              !is.na(bait) & !is.na(treatment) ~ "bait_treatment",
                              !is.na(bait) ~ "bait",
                              !is.na(treatment) ~ "treatment")) %>%
  dplyr::mutate(prior_tau = case_when(effect_type == "bait" ~ 2.0,
                                      effect_type == "treatment" ~ 1.0,
                                      effect_type == "bait_treatment" ~ 0.5,
                                      TRUE ~ 1.0),
                prior_mean = 0.0,
                prior_df1 = case_when(effect_type == "bait" ~ 1.0,
                                      effect_type == "treatment" ~ 2.0,
                                      effect_type == "bait_treatment" ~ 4.0,
                                      TRUE ~ 1.0),
                prior_df2 = prior_df1)
effects.df$effect_label <- sapply(1:nrow(effects.df), function(i) {
  bait <- as.character(effects.df$bait[[i]])
  trt <- as.character(effects.df$treatment[[i]])
  chunks <- if (!is.na(bait)) { bait } else { character(0) }
  if (!is.na(trt)) chunks <- append(chunks, trt)
  paste0(chunks, collapse="+")
})

all_metaconditions <- c(all_conditions)
conditionXmetacondition.mtx <- msglm::constant_matrix(FALSE, list(condition = all_conditions,
                                                                  metacondition = all_metaconditions))
for (cname in conditions.df$condition) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}
conditionXmetacondition.df <- as.data.frame.table(conditionXmetacondition.mtx, responseName = "part_of",
                                                  stringsAsFactors = FALSE) %>% as_tibble() %>%
  dplyr::filter(part_of) %>% dplyr::select(-part_of)
pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plot_path,
                              paste0(project_id, "_metaconditions_", data_info$msfolder, "_", data_info$fit_ver, ".pdf")),
         width = 6, height = 5)

canx_vs_empty <- paste0("CANX_vs_empty@", levels(conditions.df$treatment))
zikv_vs_mock <- paste0("ZikV_vs_Mock@", levels(conditions.df$bait))
all_contrasts <- c(canx_vs_empty, zikv_vs_mock)
metaconditionXcontrast.mtx <- msglm::constant_matrix(0, list(metacondition = all_metaconditions,
                                                             contrast = all_contrasts))
for (bait in levels(conditions.df$bait)) {
  metaconditionXcontrast.mtx[paste0(bait, '_ZikV'), paste0('ZikV_vs_Mock@', bait)] <- 1.0
  metaconditionXcontrast.mtx[paste0(bait, '_Mock'), paste0('ZikV_vs_Mock@', bait)] <- -1.0
}
for (trt in levels(conditions.df$treatment)) {
  metaconditionXcontrast.mtx[paste0('CANX_', trt), paste0('CANX_vs_empty@', trt)] <- 1.0
  metaconditionXcontrast.mtx[paste0('empty_', trt), paste0('CANX_vs_empty@', trt)] <- -1.0
}
pheatmap(metaconditionXcontrast.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plot_path,
                              paste0(project_id, "_contrasts_", data_info$msfolder, "_", data_info$fit_ver, ".pdf")),
         width = 6, height = 5)

metaconditionXcontrast.df <- as.data.frame.table(metaconditionXcontrast.mtx, responseName="weight",
                                                 stringsAsFactors = FALSE) %>% as_tibble() %>%
  dplyr::filter(weight != 0) %>%
  dplyr::mutate(contrast_type = if_else(contrast %in% canx_vs_empty, "filter", "comparison"),
                condition_role = if_else(contrast_type == 'filter',
                                         if_else(weight > 0, 'signal', 'background'),
                                         'signal'))

contrasts.df <- dplyr::select(metaconditionXcontrast.df, contrast, contrast_type) %>%
    dplyr::distinct()

msrunXbatchEffect.mtx <- model.matrix(~ 1 + batch, msdata_full$msruns %>% dplyr::mutate(batch = factor(batch)))
msrunXbatchEffect.mtx <- msrunXbatchEffect.mtx[, -1, drop=FALSE]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata_full$msruns$msrun,
                                        batch_effect = colnames(msrunXbatchEffect.mtx))
pheatmap(msrunXbatchEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         file=file.path(analysis_path, "plots", str_c(msfolder, "_", data_info$fit_ver),
                        paste0(project_id, "_", msfolder, "_", data_info$fit_ver, "_batch_effects.pdf")),
         width=4, height=10)

batch_effects.df <- tibble(batch_effect = colnames(msrunXbatchEffect.mtx),
                           is_positive = FALSE,
                           prior_mean = 0.0)

msglm_def <- msglm_model(conditionXeffect.mtx, conditions.df, effects.df, verbose=TRUE) %>%
  msglm::set_contrasts(metaconditionXcontrast.mtx, conditionXmetacondition.mtx, contrasts.df) %>%
  msglm::set_batch_effects(msrunXbatchEffect.mtx,
                           batch_effects.df, applies_to="object")

msdata <- import_msglm_data(msdata_full, msglm_def,
                            object="protregroup", quantobject="pepmodstate",
                            verbose = TRUE)

msdata4norm.df <- msdata_full$pepmodstate_intensities %>%
    dplyr::semi_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id)) %>%
    dplyr::semi_join(dplyr::semi_join(msdata_full$pepmods,
        dplyr::filter(msdata_full$protregroups %>% dplyr::mutate(is_viral=FALSE),
                      !is_reverse & !is_contaminant & !is_viral) %>%
        dplyr::inner_join(msdata_full$protregroup2pepmod) %>% dplyr::select(pepmod_id)) %>%
        dplyr::filter(!is_reverse & !is_contaminant) %>%
        dplyr::select(pepmod_id))

require(cmdstanr)
options(mc.cores=8L)
Sys.setenv(MKL_NUM_THREADS=1L)

# normalize experiments:
msruns_hnorm <- multilevel_normalize_experiments(msdata$pepmodstate_mscalib,
    dplyr::select(msdata$msruns, msrun, treatment, condition, bait),
    msdata4norm.df,
    quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
    #mcmc.iter = 100L,
    verbose=TRUE,
    stan_method = "mcmc",
    quant_ratio.max = 3.0,
    norm_levels = list(msrun = list(cond_col="msrun", nmschan_ratio.min=1.0, max_objs=1000L),
                       condition = list(cond_col="condition", max_objs=1000L),
                       bait = list(cond_col = "bait", nmschan_ratio.min=1.0)))
# skip bait shift because it's too big
msdata$msrun_shifts <- msruns_hnorm$mschannel_shifts

# estimate the effect scales per each msrun
pepmodstate_intensities4eff_scale.df <- msdata_full$pepmodstate_intensities %>%
  # fitler for CANX interactors that are known ER-localized proteins
  dplyr::semi_join(dplyr::filter(msdata_full$proteins, str_detect(gene_name, "^(CANX|EMC(1|2|3|4|10)|ERLIN(1|2)|GOLGA4|ERP29|PGRMC1|PDIA6|ESYT(1|2)|POR|TXNDC(5|12)|MLEC|MESD|SRPR[AB]|UGGT1|NOMO1|MOGS|SPTLC1|ASPH)$")) %>%
  dplyr::inner_join(msdata_full$protein2protregroup, by="protein_ac") %>%
    dplyr::inner_join(dplyr::filter(msdata_full$protregroup2pepmodstate, is_specific), by="protregroup_id") %>%
    dplyr::select(pepmodstate_id) %>% dplyr::distinct()) %>%
  dplyr::inner_join(msdata$msrun_shifts) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::mutate(intensity_norm = intensity*2^(-total_msrun_shift),
                intensity_norm_log2 = log2(intensity_norm)) %>%
  dplyr::group_by(pepmodstate_id, condition, batch) %>%
  dplyr::mutate(intensity_cond = median(intensity_norm, na.rm=TRUE)) %>%
  dplyr::group_by(pepmodstate_id, batch) %>%
  dplyr::mutate(intensity_empty_mock = intensity_cond[[which.max(condition == "empty_Mock")]],
                intensity_empty_zikv = intensity_cond[[which.max(condition == "empty_ZikV")]],
                intensity_CANX_mock = intensity_cond[[which.max(condition == "CANX_Mock")]],
                intensity_CANX_zikv = intensity_cond[[which.max(condition == "CANX_ZikV")]]) %>%
  dplyr::ungroup() %>%
  dplyr::filter(bait == "CANX") %>%
  dplyr::mutate(is_infected = treatment != "Mock",
                intensity_cond_log2 = log2(intensity_cond),
                intensity_cond_vs_empty_log2 = intensity_cond_log2 - log2(if_else(is_infected, intensity_empty_zikv, intensity_empty_mock)),
                effect = if_else(intensity_cond_vs_empty_log2 == 0, NA_real_, intensity_cond_vs_empty_log2)) %>%
  #dplyr::group_by(is_infected, batch) %>%
  #dplyr::mutate(effect_bin = as.integer(cut(effect, 50))) %>%
  #dplyr::group_by(is_infected, batch, effect_bin) %>%
  #dplyr::slice_sample(n = 200) %>%
  dplyr::ungroup()

pepmodstate_intensities_pairs4eff_scale.df <- dplyr::left_join(pepmodstate_intensities4eff_scale.df,
                                                               pepmodstate_intensities4eff_scale.df,
                                                          by=c("bait", "batch", "pepmodstate_id")) %>%
  dplyr::filter(msrun.x != msrun.y & (effect.x > 1.0) & (effect.y > 1.0))

ggplot(dplyr::filter(pepmodstate_intensities_pairs4eff_scale.df, is_infected.y & msrun.x != msrun.y)) +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dotted") +
  geom_point(aes(x = effect.x, y = effect.y, color = factor(batch)))

require(cmdstanr)
options(mc.cores=8)
effect_scales.stanmodel <- cmdstanr::cmdstan_model(file.path(project_scripts_path, "effect_scales_normalize.stan"))

fit_effect_scales <- function(data_df, sigma_t_offset=1E-3, scale_sigma_df=2, data_sigma_df=16) {
  data_df <- mutate(data_df,
                    msrun.x = factor(msrun.x),
                    msrun.y = factor(msrun.y, levels=levels(msrun.x)),
                    msrun_ix.x = as.integer(msrun.x),
                    msrun_ix.y = as.integer(msrun.y))

  effect_scales.data <- list(Nmsruns = n_distinct(data_df$msrun.x),
                             NeffectPairs = nrow(data_df),
                             outlierProb = 1E-3,
                             effect_x = data_df$effect.x,
                             effect_y = data_df$effect.y,
                             msrun_x = data_df$msrun_ix.x,
                             msrun_y = data_df$msrun_ix.y,
                             sigma_t_offset = sigma_t_offset,
                             scale_sigma_df = scale_sigma_df,
                             data_sigma_df = data_sigma_df)
  effect_scales.stanfit <- effect_scales.stanmodel$sample(
                              data = effect_scales.data, chains = 8,
                              iter_warmup=1000, iter_sampling=1000)
  effect_scales_stats.df <- effect_scales.stanfit$summary(
        variables = c("data_sigma", "scale_sigma", "scale_log2", "scale_mult"),
        msglm:::posterior_summary_metrics) %>%
    dplyr::rename(varspec = variable) %>%
    tidyr::extract(varspec, c("var", "index_msrun"), "^(\\w+)(?:\\[(\\d+)\\])?$",
                   remove=FALSE, convert=TRUE) %>%
    dplyr::mutate(msrun = levels(data_df$msrun.x)[index_msrun])
}

effect_scales_fit.df <- fit_effect_scales(pepmodstate_intensities_pairs4eff_scale.df)
#saveRDS(effect_scales_fit.df, file.path(scratch_path, str_c(project_id, "_", msfolder, "_", data_info$fit_ver, "_effect_scales.rds")))

effect_scales.df <- dplyr::transmute(dplyr::filter(effect_scales_fit.df, var=="scale_log2"),
                                     msrun,
                                     msrun_scale_log2=median, msrun_scale_log2_sd=sd,
                                     ess_bulk, rhat) %>%
  dplyr::left_join(dplyr::select(msdata_full$msruns, msrun, bait, condition, treatment, batch)) %>%
  dplyr::group_by(bait, batch) %>%
  dplyr::mutate(msrun_scale_log2 = msrun_scale_log2 - max(msrun_scale_log2),
                msrun_scale=2^msrun_scale_log2,
                msrun_scale_log2_median = median(msrun_scale_log2),
                msrun_scale_avg = mean(msrun_scale)) %>%
  dplyr::ungroup()

msrunXeffect.df <- inner_join(as.data.frame.table(conditionXeffect.mtx, responseName="w"), msdata$msruns) %>%
  dplyr::filter(w != 0) %>%
  dplyr::left_join(dplyr::select(effect_scales.df, msrun, msrun_scale)) %>%
  dplyr::left_join(dplyr::select(effects.df, effect, eff_bait = bait)) %>%
  dplyr::mutate(effect_scale = if_else(is.na(eff_bait), 1.0, msrun_scale),
                w = w * effect_scale)
msexpXeffect.mtx <- msglm::frame2matrix(msrunXeffect.df, row_col="msexperiment", col_col="effect",
                                        cols = colnames(conditionXeffect.mtx),
                                        rows = dplyr::arrange(msdata$msexperiments, bait, treatment, batch) %>% .$msexperiment)

pheatmap(msexpXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers = TRUE, number_format = "%.3f",
         filename = file.path(analysis_path, "plots", paste0(msfolder, "_", data_info$fit_ver),
                              paste0(project_id, "_msexperimentXeffect_", msfolder, "_", data_info$fit_ver, ".pdf")),
         width = 8, height = 12)

# redefine model with scaled effects
msglm_def <- msglm_model(conditionXeffect.mtx, conditions.df, effects.df,
                         msprobeXeffect = msexpXeffect.mtx, verbose=TRUE) %>%
  msglm::set_contrasts(metaconditionXcontrast.mtx, conditionXmetacondition.mtx, contrasts.df) %>%
  msglm::set_batch_effects(msrunXbatchEffect.mtx,
                           batch_effects.df, applies_to="object")

rmsglmdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', data_info$fit_ver, '.RData'))
message('Saving data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata, msglm_def, msruns_hnorm, effect_scales.df,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', data_info$data_ver, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     file = rfulldata_filepath)

message('Done.')

tmp.env <- new.env(parent=baseenv())
load("/pool/analysis/astukalov/aherrmann_cov2omap/scratch/aherrmann_cov2omap_msglm_data_20211003.RData", envir = tmp.env)
ggplot(inner_join(msdata$msrun_shifts, tmp.env$msdata$msrun_shifts, by=c("conditionXmsfraction", "compartmentXmsfraction", "msrun")),
       aes(x=total_msrun_shift.x, y=total_msrun_shift.y, label=msrun, color=compartmentXmsfraction)) +
  geom_point() +
  geom_text(vjust=-1.1, size=3) +
  geom_abline(intercept = 0, slope = 1, color="firebrick", linetype="dotted") +
  theme_bw()

ggplot(total_msrun_shifts.df %>%
       dplyr::inner_join(total_msrun_shifts_old.df %>%
       dplyr::select(mschannel_ix, total_mschannel_shift_old = total_mschannel_shift))) +
  geom_point()
