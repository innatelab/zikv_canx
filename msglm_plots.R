# plotting the volcanos for contrasts/effects
# peptide heatmaps
# how model fits the data
# Alexey Stukalov 2021.05.03

#The following part is only needed if starting from a fresh environment----
project_id <- "sdenolly_canxapms"
fit_version <- "20211126"
message('Project ID=', project_id, " fit version=", fit_version)

require(tidyverse)
require(rlang)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

message('Loading data...')
load(file.path(scratch_path, str_c(project_id, '_msglm_fit_', fit_version, '.RData')))
#load(file.path(scratch_path, str_c(project_id, '_msglm_fit_meanfield_', fit_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', fit_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', results_info$data_ver, '.RData')))

modelobj <- msdata$msentities['object']
quantobj <- msdata$msentities['quantobject']
obj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")
obj_labu_shift <- msdata[[str_c(quantobj, "_mscalib")]]$zShift

# The plotting starts from here

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)
require(ggpubr)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(misc_scripts_path, 'furrr_utils.R'))

treatment_palette <- c("Mock"="gray", "ZikV" = "#F4982A")
hit_palette <- c("non-hit"="grey", hit="black", viral="#F4982A", "viral hit"="#F4982A")
base_font_family <- "Segoe UI Symbol"
base_plot_path <- file.path(analysis_path, 'plots', str_c(data_info$msfolder, "_", fit_version))#, "_meanfield"))
#sel_ci_target <- "replicate"
sel_ci_target <- "average"

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
                                 TRUE ~ NA_real_),
    median_max = case_when(contrast_type=="filter" ~ 5,
                           contrast_type=="comparison" ~ 4,
                           TRUE ~ NA_real_)
  )

# Volcano plots of all contrasts
object_contrasts_4show.df <- fit_contrasts$object_conditions %>%
  dplyr::filter(var %in% c('obj_cond_labu', 'obj_cond_labu_replCI')) %>%
  dplyr::ungroup() %>%
  select(-contains('threshold')) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_sel = FALSE,
                is_hit_nomschecks = is_signif & !is_reverse & !is_contaminant,
                is_hit = is_hit_nomschecks,
                hit_type = case_when(is_hit & is_viral ~ "viral",
                                     is_hit ~ "hit", TRUE ~ "non-hit"),
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset,
                truncation = volcano_truncation(median, median_trunc, p_value, is_hit_nomschecks, !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                show_label = coalesce(is_hit_nomschecks, FALSE))

require(furrr)
plan(multicore, workers=16)
plot_furrr_opts <- furrr_options(globals = c("base_plot_path", "base_font_family",
                                             "project_id", "data_version", "fit_version", "obj_labu_shift",
                                             "mlog10_trans", "mlog_pow_trans", "mlog_breaks",
                                             "theme_bw_ast",
                                             "point_truncation_shape_palette", "point_truncation_size_palette",
                                             "treatment_palette",
                                             "sel_ci_target", "modelobj", "quantobj", "modelobj_idcol", "quantobj_idcol",
                                             "msglm_def", "fit_stats", "msdata", "msdata_full", "modelobject_intensities.df"),
                                 packages = c("dplyr", "ggplot2", "Cairo", "ggrepel", "ggrastr", "stringr"),
                                 stdout=TRUE)

group_by(object_contrasts_4show.df, ci_target, contrast,
         offset, median_threshold, p_value_threshold) %>%
#group_walk(.keep=TRUE,
future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
                  function(sel_object_contrast.df, contrast_info) {
  message("Plotting ", contrast_info$contrast, " ci_target=", contrast_info$ci_target)

  p <- ggplot(sel_object_contrast.df,
              aes(x=median_trunc, y=p_value, shape=truncation, size=truncation_type, color=hit_type)) +
    geom_hline(data=contrast_info, aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
    #geom_hline(data=contrast_info, aes(yintercept = p_value_max), linetype=1, color="darkgray") +
    geom_vline(data=contrast_info, aes(xintercept = offset), linetype=1, color="darkgray") +
    geom_vline(data=contrast_info, aes(xintercept = offset + median_threshold), linetype=2, color="darkgray") +
    geom_vline(data=contrast_info, aes(xintercept = offset - median_threshold), linetype=2, color="darkgray") +
    geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                    alpha=0.1, size=0.5, color="darkgray") +
    geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif)) +
    geom_text_repel(data=dplyr::filter(sel_object_contrast.df, is_hit_nomschecks) %>% dplyr::arrange(is_hit),
                    aes(label = object_label),
                    nudge_y = -0.12,
                    size=2.5, point.padding=0.1, box.padding = 0.1,
                    max.overlaps = 20,
                    show.legend = FALSE, segment.color = "gray") +
    scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
    #scale_fill_gradient(low="gray75", high="black") +
    #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
    scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
    scale_size_manual(values=point_truncation_size_palette, guide="none") +
    scale_color_manual(values=hit_palette, na.value="magenta") +
    #scale_color_manual(values=orgcode_palette, guide="none") +
    #facet_grid(p_value_range ~ contrast, scales = "free_y") +
    ggtitle(contrast_info$contrast, subtitle=str_c("ci_target=", contrast_info$ci_target)) +
    theme_bw_ast(base_family = base_font_family)
  plot_path <- file.path(base_plot_path, str_c("volcanos_contrasts_", contrast_info$ci_target))
  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

  ggsave(filename = file.path(plot_path,
                              str_c(project_id, '_', fit_version, '_volcano_',
                                    str_replace_all(contrast_info$contrast, ":|@", "_"), '.pdf')),
         plot=p, width=15, height=18, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
})

#### volcano plots for effects

object_effects_thresholds.df <- dplyr::select(msglm_def$effects, effect, prior_mean) %>%
  dplyr::mutate(p_value_threshold = 1E-3,
                median_threshold = c(1.0),
                median_max = c(5)
  )

object_effects_4show.df <- fit_stats$object_effects %>%
  dplyr::filter(var %in% c('obj_effect', 'obj_effect_replCI')) %>%
  select(-contains('threshold')) %>%
  dplyr::inner_join(dplyr::select(msglm_def$effects, effect, treatment, bait)) %>%
  #dplyr::left_join(dplyr::select(fit_contrasts$object_conditions, ci_target, contrast, object_id,
  #                               contrast_median = median, contrast_pvalue = p_value) %>%
  #                  dplyr::inner_join(dplyr::select(contrasts.df, contrast, bait,
  #                                                  treatment = treatment_lhs))) %>%
  dplyr::inner_join(object_effects_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - prior_mean) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_reverse & !is_contaminant,
                is_hit = is_hit_nomschecks,# & (is.na(treatment) | (contrast_median * median < 0 | contrast_pvalue >= 0.05)),
                hit_type = case_when(is_hit & is_viral ~ "viral",
                                     is_hit ~ "hit", TRUE ~ "non-hit"),
                show_label = coalesce(is_hit_nomschecks, FALSE),
                mean_trunc = pmax(-median_max, pmin(median_max, mean - prior_mean)) + prior_mean,
                median_trunc = pmax(-median_max, pmin(median_max, median - prior_mean)) + prior_mean,
                truncation = volcano_truncation(median, median_trunc, p_value, is_hit, !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif))

group_by(object_effects_4show.df, ci_target, effect, effect_label,
         prior_mean, median_threshold, p_value_threshold) %>%
future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
#group_walk(.keep=TRUE,
                   function(sel_object_effect.df, effect_info) {
  message("Plotting ", effect_info$effect, " ci_target=", effect_info$ci_target)

  p <- ggplot(sel_object_effect.df,
              aes(x=median_trunc, y=p_value, shape=truncation, size=truncation_type, color=hit_type)) +
    geom_hline(data=effect_info, aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
    #geom_hline(data=effect_info, aes(yintercept = p_value_max), linetype=1, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean), linetype=1, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean + median_threshold), linetype=2, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean - median_threshold), linetype=2, color="darkgray") +
    geom_point_rast(data=dplyr::filter(sel_object_effect.df, !is_signif),
                    alpha=0.1, size=0.5, color="darkgray") +
    geom_point(data=dplyr::filter(sel_object_effect.df, is_signif)) +
    geom_text_repel(data=dplyr::filter(sel_object_effect.df, is_hit_nomschecks) %>% dplyr::arrange(is_hit),
                    aes(label = object_label),
                    nudge_y = -0.12, force=0.25, max.overlaps=20,
                    size=2.5, point.padding=0.2, box.padding=0.15,
                    show.legend = FALSE, segment.color = "gray") +
    scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
    #scale_fill_gradient(low="gray75", high="black") +
    #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
    scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
    scale_size_manual(values=point_truncation_size_palette, guide="none") +
    scale_color_manual(values=hit_palette, na.value="magenta") +
    #scale_color_manual(values=orgcode_palette, guide="none") +
    #facet_grid(p_value_range ~ contrast, scales = "free_y") +
    ggtitle(effect_info$effect_label, subtitle=str_c("ci_target=", effect_info$ci_target)) +
    theme_bw_ast(base_family = base_font_family)
  plot_path <- file.path(base_plot_path,
                         str_c("volcanos_effects_", effect_info$ci_target))
  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

  ggsave(filename = file.path(plot_path,
                              str_c(project_id, '_', fit_version, '_effect_volcano_',
                                    str_replace_all(effect_info$effect, ":|@", "_"), '.pdf')),
         plot = p, width=15, height=18, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
})

#### volcano plots for batch effects

object_batch_effects_thresholds.df <- dplyr::select(msglm_def$batch_effects, batch_effect, prior_mean) %>%
  dplyr::mutate(p_value_threshold = 1E-3,
                median_threshold = c(0.5),
                median_max = c(5)
  )

object_batch_effects_4show.df <- fit_stats$object_batch_effects %>%
  dplyr::filter(var %in% c('obj_batch_effect')) %>%
  select(-contains('threshold')) %>%
  dplyr::inner_join(object_batch_effects_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - prior_mean) >= median_threshold),
                hit_type = case_when(is_signif & is_viral ~ "viral hit",
                                     is_signif & !is_contaminant & !is_reverse ~ "hit",
                                     is_signif ~ "non-hit",
                                     TRUE ~ "non-hit"),
                mean_trunc = pmax(-median_max, pmin(median_max, mean - prior_mean)) + prior_mean,
                median_trunc = pmax(-median_max, pmin(median_max, median - prior_mean)) + prior_mean,
                truncation = volcano_truncation(median, median_trunc, p_value, is_signif, !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif))

group_by(object_batch_effects_4show.df, batch_effect,
         prior_mean, median_threshold, p_value_threshold) %>%
future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
#group_walk(.keep=TRUE,
                   function(sel_object_effect.df, effect_info) {
  message("Plotting ", effect_info$batch_effect)

  p <- ggplot(sel_object_effect.df,
              aes(x=median_trunc, y=p_value, color=hit_type,
                  shape=truncation, size=truncation_type)) +
    geom_hline(data=effect_info, aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
    #geom_hline(data=effect_info, aes(yintercept = p_value_max), linetype=1, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean), linetype=1, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean + median_threshold), linetype=2, color="darkgray") +
    geom_vline(data=effect_info, aes(xintercept = prior_mean - median_threshold), linetype=2, color="darkgray") +
    geom_point_rast(data=dplyr::filter(sel_object_effect.df, !is_signif),
                    alpha=0.1, size=0.5, color="darkgray") +
    geom_point(data=dplyr::filter(sel_object_effect.df, is_signif)) +
    geom_text_repel(data=dplyr::filter(sel_object_effect.df, is_signif),
                    aes(label = object_label),
                    nudge_y = -0.12, force=0.25, max.overlaps=20,
                    size=2.5, point.padding=0.2, box.padding=0.15,
                    show.legend = FALSE, segment.color = "gray") +
    scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
    scale_color_manual(values=hit_palette, na.value="magenta") +
    scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
    scale_size_manual(values=point_truncation_size_palette, guide="none") +
    ggtitle(effect_info$batch_effect) +
    theme_bw_ast(base_family = base_font_family)
  plot_path <- file.path(base_plot_path, "volcanos_batch_effects")
  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

  ggsave(filename = file.path(plot_path,
                              str_c(project_id, '_', fit_version, '_batch_effect_volcano_',
                                    str_replace_all(effect_info$batch_effect, ":|@", "_"), '.pdf')),
         plot = p, width=15, height=18, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
})

object_treatment_specificity.df <-
  dplyr::filter(fit_contrasts$object_conditions, var %in% c('obj_cond_labu', 'obj_cond_labu_replCI')) %>%
  dplyr::mutate(contrast_type = if_else(contrast_type == "filtering", "filter", contrast_type)) %>%
  dplyr::inner_join(dplyr::select(dplyr::filter(contrasts.df, contrast_type == "filter" & treatment_lhs == "Mock"), contrast)) %>%
  dplyr::select(-contains('threshold')) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (median - offset >= median_threshold),
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset) %>%
  dplyr::rename_with(~paste0(.x, "_mock"), matches("^q\\d|^ess_|^prob_|_threshold$|^offset|contrast|median|mean") |
                      any_of(c("mean", "median", "sd", "mad", "p_value", "is_signif",
                              "rhat", "bandwidth", "bin_width", "nperms"))) %>%
  dplyr::inner_join(
    dplyr::filter(fit_contrasts$object_conditions, var %in% c('obj_cond_labu', 'obj_cond_labu_replCI')) %>%
    dplyr::mutate(contrast_type = if_else(contrast_type == "filtering", "filter", contrast_type)) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(contrasts.df, contrast_type == "filter" & treatment_lhs == "ZikV"), contrast)) %>%
    dplyr::select(-contains('threshold')) %>%
    dplyr::inner_join(object_contrasts_thresholds.df) %>%
    dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (median - offset >= median_threshold),
                  mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                  median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset) %>%
    dplyr::rename_with(~paste0(.x, "_zikv"), matches("^q\\d|^ess_|^prob_|_threshold$|^offset|contrast|median|mean") |
                       any_of(c("mean", "median", "sd", "mad", "p_value", "is_signif",
                                "rhat", "bandwidth", "bin_width", "nperms")))) %>%
  dplyr::inner_join(
    dplyr::filter(fit_stats$object_effects, var %in% c('obj_effect', 'obj_effect_replCI')
                  & !is.na(bait) & !is.na(treatment)) %>%
    dplyr::select(-contains('threshold')) %>%
    dplyr::inner_join(object_effects_thresholds.df) %>%
    dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - prior_mean) >= median_threshold)) %>%
    dplyr::rename_with(~paste0(.x, "_eff"), matches("^q\\d|^ess_|^prob_|_threshold$|^offset|contrast|median") |
                       any_of(c("var", "mean", "median", "sd", "mad", "p_value", "is_signif",
                                "rhat", "bandwidth", "bin_width", "nperms")))) %>%
  dplyr::mutate(is_signif = is_signif_eff | is_signif_mock | is_signif_zikv,
                is_sel = FALSE,
                is_hit_nomschecks = is_signif & !is_reverse & !is_contaminant,
                is_hit = is_signif_eff & (is_signif_mock | is_signif_zikv) & !is_reverse & !is_contaminant,
                is_viral = object_id %in% dplyr::filter(msdata$objects, is_viral)$object_id,
                hit_type = case_when(is_hit & is_viral ~ "viral hit",
                                     is_hit ~ "hit",
                                     is_viral ~ "viral", TRUE ~ "non-hit"),
                truncation = scatter_truncation(median_mock, median_trunc_mock,
                                                median_zikv, median_trunc_zikv,
                                                is_hit | !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                show_label = coalesce(is_hit_nomschecks, FALSE))

group_by(object_treatment_specificity.df, ci_target, median_threshold_zikv, median_threshold_mock, median_threshold_eff) %>%
#future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts
group_walk(.keep=TRUE,
                   function(plot_data.df, plot_info) {
  message("Plotting ci_target=", plot_info$ci_target)
  p <- ggplot(plot_data.df,
              aes(x=median_trunc_mock, y=median_trunc_zikv, shape=truncation,
                  size=truncation_type, color=hit_type)) +
    geom_hline(yintercept = 0, linetype=1, color="darkgray") +
    geom_vline(xintercept = 0, linetype=1, color="darkgray") +
    geom_hline(yintercept = c(-plot_info$median_threshold_mock, +plot_info$median_threshold_mock), linetype=2, color="darkgray") +
    geom_vline(xintercept = c(-plot_info$median_threshold_zikv, +plot_info$median_threshold_zikv), linetype=2, color="darkgray") +
    geom_abline(intercept = 0, slope = 1, linetype=1, color="darkgray") +
    geom_abline(intercept = c(-plot_info$median_threshold_eff, +plot_info$median_threshold_eff), slope = 1, linetype=2, color="darkgray") +
    geom_point_rast(data=dplyr::filter(plot_data.df, !is_signif),
                    alpha=0.1, size=0.5, color="darkgray") +
    geom_point(data=dplyr::filter(plot_data.df, is_signif)) +
    geom_text_repel(data=dplyr::filter(plot_data.df, is_hit_nomschecks),
                    aes(label = object_label),
                    nudge_y = -0.05,
                    size=2.5, force=1.0,
                    point.padding=0.01, box.padding=0.1, max.overlaps = 25,
                    show.legend = FALSE, segment.color = "gray") +
    scale_x_continuous("log2(fold-change) Mock",  limits = c(-1.15, 5.15)) +
    scale_y_continuous("log2(fold-change) ZikV",  limits = c(-1.15, 5.15)) +
    scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
    scale_size_manual(values=point_truncation_size_palette, guide="none") +
    scale_color_manual(values=hit_palette, na.value="magenta") +
    ggtitle(paste0("ZikV vs Mock of CANX AP-MS"),
            subtitle=str_c("ci_target=", plot_info$ci_target)) +
    theme_bw_ast(base_family = base_font_family)
  plot_path <- file.path(base_plot_path,
                         str_c("scatter_compartment_", plot_info$ci_target))
  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

  ggsave(filename = file.path(plot_path,
                              str_c(project_id, '_', fit_version, '_zikv_vs_mock_scatter.pdf')),
         plot = p, width=15, height=14, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
})

#Making boxplots for all proteins (the original script makes timecourse for each protein). The dots are from the original LFQ values of MQ.----
sel_objects.df <- dplyr::filter(msdata$objects, str_detect(gene_names, "^VPS\\d+|STAT\\d+|ZC3HAV")) #This is used for debugging
sel_objects.df <- dplyr::semi_join(msdata$objects, dplyr::select(fit_stats$objects, object_id), by="object_id") #This is all!

# boxplots of model fit for each protein group
rowwise(sel_objects.df) %>%
#group_walk(
future_group_walk(.progress=TRUE, .options=plot_furrr_opts,
                  function(sel_obj.df, ignore.df) {
  obj_label <- sel_obj.df$object_label
  obj_label_safe <- str_replace_all(str_remove(coalesce(obj_label, "noname"), "\\.\\.\\.$"), "[/.]+", "_")
  message("Plotting ", obj_label, " box plots")
  sel_obj_conds.df <- dplyr::semi_join(dplyr::filter(fit_stats$object_conditions, ci_target == sel_ci_target & str_starts(var, "obj_cond_labu")),
                                       sel_obj.df, by="object_id") %>%
    dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, bait, treatment), by="condition") %>%
    dplyr::mutate(dplyr::across(c(mean, median, starts_with("q")),
                  ~2^(.x + obj_labu_shift)))%>%
    dplyr::mutate(q97.5_limit = max(median + (q75-q25)*5)) %>%
    #dplyr::mutate(bait = factor(bait) %>% relevel("empty")) %>%
    dplyr::arrange(bait, treatment) %>%
    dplyr::mutate(condition = factor(condition, levels=condition))
  sel_obj_contrasts.df <- dplyr::semi_join(dplyr::filter(object_contrasts_4show.df, ci_target == sel_ci_target & str_starts(var, "obj_cond_labu")),
                                                         sel_obj.df, by="object_id") %>%
                        dplyr::left_join(contrasts.df) %>%
                        dplyr::left_join(dplyr::select(sel_obj_conds.df, obj_abu_lhs = q97.5, metacondition_lhs = condition)) %>%
                        dplyr::left_join(dplyr::select(sel_obj_conds.df, obj_abu_rhs = q97.5, metacondition_rhs = condition)) %>%
                        mutate(obj_abu = pmax(obj_abu_lhs, obj_abu_rhs),
                               group1 = metacondition_lhs, group2 = metacondition_rhs,
                               y.position = log10(obj_abu) + 0.1, treatment = NA_character_,
                               p_value = sprintf("%.3f p=%.2e", median, p_value))
  #sel_obj_msdata.df <- dplyr::semi_join(modelobject_intensities.df, sel_obj.df, by="object_id")

  if (nrow(sel_obj_conds.df) > 0) {
    obj_plot <-
      ggplot(data=sel_obj_conds.df,
             aes(x=condition, color=treatment, fill=treatment)) +
      geom_boxplot(aes(middle=median, lower=q25, upper=q75, ymin=q2.5, ymax=pmin(q97.5, q97.5_limit)),
                   alpha=0.5, stat = "identity") +
      stat_pvalue_manual(data = dplyr::filter(sel_obj_contrasts.df, !is_signif),
                         #aes(y.position = obj_abu,
                         #    xmin = metacondition_lhs, xmax = metacondition_rhs),
                         color="darkgray", label = "p_value", size=2,
                         tip.length = 0.005, step.increase = 0.01) +
      stat_pvalue_manual(data = dplyr::filter(sel_obj_contrasts.df, is_signif),
                         #aes(y.position = obj_abu,
                         #    xmin = metacondition_lhs, xmax = metacondition_rhs),
                         color="black", label = "p_value", size=2,
                         tip.length = 0.005, step.increase = 0.01) +
      scale_colour_manual(values=treatment_palette) +
      scale_fill_manual(values=treatment_palette) +
      scale_linetype_manual(values=c("AP"="solid", "proteome"="dotted")) +
      #geom_point(data=sel_obj_msdata.df, aes(y = intensity_norm),
      #           position = position_jitterdodge(jitter.width = 0.25, dodge.width=1.0), size=0.5) +
      theme_bw_ast(base_family = base_font_family, base_size = 8) +
      theme(legend.position = "bottom") +
      scale_y_log10() +
      ggtitle(str_c(obj_label, " fit"),
              subtitle=sel_obj.df$protein_description)

    plot_path <- file.path(base_plot_path, str_c(modelobj, "s_fit_", sel_ci_target))
    if (!dir.exists(plot_path)) dir.create(plot_path)
    ggsave(obj_plot, file = file.path(plot_path,
                               str_c(project_id, "_", fit_version, "_", obj_label_safe,
                                     "_", sel_obj.df$object_id[[1]], ".pdf")),
           width=4, height=6, device = cairo_pdf, family=base_font_family)
  }
})

#Making peptide heatmaps for every protein group----
#sel_objects.df <- dplyr::filter(msdata$objects, is_viral)
sel_objects.df <- msdata$objects # all!!!

sel_pepmodstates.df <- dplyr::inner_join(sel_objects.df, msdata_full[[str_c(modelobj, "2pepmod")]]) %>%
  #dplyr::filter(is_specific) %>%
  dplyr::inner_join(msdata_full$pepmodstates) %>%
  dplyr::inner_join(dplyr::transmute(msdata_full$pepmodstate_intensities, pepmodstate_id,
                                     is_idented, is_quanted = !is.na(intensity)) %>%
                    dplyr::group_by(pepmodstate_id) %>%
                    dplyr::summarise(nidented = sum(is_idented, na.rm=TRUE),
                                     nquanted = sum(is_quanted, na.rm=TRUE),
                                     .groups="drop")) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id, pepmod_seq, peptide_seq)) %>%
  dplyr::select(object_id, pepmod_id, majority_protein_acs, protac_label, gene_label, gene_names,
                nidented, nquanted, pepmod_seq, peptide_seq, charge, pepmodstate_id, is_specific)
if (exists("fit_stats") && rlang::has_name(fit_stats, "quantobjects")) {
  sel_pepmodstates.df <- dplyr::left_join(sel_pepmodstates.df,
                                          dplyr::select(dplyr::filter(fit_stats$quantobjects, var=="qobj_shift"),
                                                        pepmodstate_id=quantobject_id, pms_median = median),
                                          by="pepmodstate_id") %>%
    dplyr::arrange(object_id, desc(is_specific), desc(coalesce(pms_median, -1000.0)), desc(nidented), desc(nquanted),
                   peptide_seq, pepmod_seq, charge)
} else {
  sel_pepmodstates.df <- dplyr::arrange(sel_pepmodstates.df, object_id, desc(is_specific), desc(nidented), desc(nquanted),
                                        peptide_seq, pepmod_seq, charge) %>%
    mutate(pms_median = NA_real_)
}
sel_pepmodstates.df <- dplyr::mutate(sel_pepmodstates.df,
                                     pepmodstate_ext = paste0(str_trunc(pepmod_seq, 60), ".", charge,
                                                              " (", pepmodstate_id, ") ",
                                                              if_else(is.na(pms_median),
                                                                      if_else(is_specific, "+", ""), "*"))
                                     %>% factor(., levels=unique(.)))

sel_pepmod_intens.df <- tidyr::crossing(pepmodstate_id = sel_pepmodstates.df$pepmodstate_id,
                                        msrun = unique(msdata_full$pepmodstate_intensities$msrun)) %>%
  dplyr::inner_join(dplyr::select(sel_pepmodstates.df, object_id, protac_label, gene_label, gene_names,
                                  pepmod_id, pepmodstate_id, charge, pepmodstate_ext, pms_median) %>%
                      dplyr::distinct()) %>%
  dplyr::left_join(dplyr::select(msdata_full$pepmodstate_intensities, pepmodstate_id, msrun, intensity, ident_type)) %>%
  dplyr::inner_join(distinct(select(msdata$msruns, msexperiment, msrun, condition, treatment, bait, batch)) %>%
                      dplyr::left_join(dplyr::select(msdata$msrun_shifts, msrun, total_msrun_shift), by="msrun") %>%
                      dplyr::arrange(bait, treatment, batch, msexperiment, msrun) %>%
                      dplyr::mutate(msrun.2 = factor(msrun, levels=unique(msrun)),
                                    msexperiment.2 = factor(msexperiment, levels=unique(msexperiment)))) %>%
  dplyr::mutate(intensity_norm = 2^(-total_msrun_shift)*intensity,
                msrun = msrun.2, msrun.2 = NULL,
                msexperiment = msexperiment.2, msexperiment.2 = NULL) %>%
  dplyr::group_by(pepmodstate_id) %>%
  dplyr::mutate(intensity_norm_trunc = pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE))) %>%
  dplyr::ungroup()

sel_pepmods.df <- dplyr::group_by(sel_pepmod_intens.df, pepmod_id) %>%
  dplyr::summarize(n_pepmod_quants = sum(!is.na(intensity)),
                   median_quant = median(intensity, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_pepmod_quants), desc(median_quant), pepmod_id)

group_by(sel_pepmod_intens.df, object_id) %>%
#group_walk(.keep=TRUE,
future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
                  function(obj_pepmod_intens.df, obj_row) {
  obj_row <- dplyr::semi_join(msdata$objects, obj_row, by="object_id")
  shown_pepmod_intens.df <- dplyr::mutate(obj_pepmod_intens.df, !!obj_idcol := object_id) %>%
                            #dplyr::group_by(msprotocol, instrument, compartment, pepmodstate_id) %>%
                            #dplyr::mutate(intensity_norm = intensity_norm / median(intensity_norm, na.rm=TRUE)) %>%
                            #dplyr::ungroup() %>%
                            dplyr::mutate(intensity_norm =
                                           pmax(pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE)),
                                                         quantile(intensity_norm, 0.05, na.rm=TRUE))) #sel_pepmod_intens.df
  obj_label <- obj_row$object_label
  obj_label_safe <- str_replace_all(str_remove(coalesce(obj_label, "noname"), "\\.\\.\\.$"), "[/.]+", "_")
  message("Plotting ", obj_label, "...")
  p <- ggplot(shown_pepmod_intens.df) +
    geom_tile(aes(x=msexperiment, y=pepmodstate_ext,
                  fill=intensity_norm, color=ident_type),
              na.rm = FALSE, size=0.5, width=0.85, height=0.85) +
    theme_bw_ast(base_family = base_font_family, base_size = 10) +
    theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0)) +
    guides(color=guide_legend("ident_type", override.aes = list(fill=NA, size=2))) +
    ggtitle(str_c(obj_label,  " (", obj_idcol, "=", obj_row$object_id,
                  ", ac=", obj_row$protac_label, ") peptide map"),
            subtitle = obj_row$protein_description) +
    scale_fill_distiller(na.value="#00000000", type="div", palette = "Spectral") +
    scale_y_discrete(breaks = levels(obj_pepmod_intens.df$pepmodstate_ext)) +
    scale_color_manual(na.value="#00000000",
                       values=c("MULTI-MSMS"="black", "MULTI-MATCH-MSMS"="khaki",
                                "ISO-MSMS"="slateblue",
                                "MSMS"="cornflowerblue", "MULTI-SECPEP"="firebrick",
                                "MULTI-MATCH"="gray"))

  plot_path <- file.path(base_plot_path, "pepmodstate_heatmaps")
  if (!dir.exists(plot_path)) dir.create(plot_path)
  ggsave(p, file = file.path(plot_path,
                             paste0(project_id, "_", results_info$data_version, "_pepmodstate_heatmap_",
                                    obj_label_safe, "_", obj_row$object_id, ".pdf")),
         width=2 + 0.25*n_distinct(shown_pepmod_intens.df$msrun),
         height=2 + 2*n_distinct(shown_pepmod_intens.df$object_id) +
                    2*min(20, 0.05*n_distinct(shown_pepmod_intens.df$pepmodstate_id)),
         device=cairo_pdf, family=base_font_family, limitsize=FALSE)
  gc()
})
