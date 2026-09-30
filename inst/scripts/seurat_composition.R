#!/usr/bin/env Rscript

# Purpose: Report how the composition of one categorical metadata column
#          (--comp_var, e.g. celltype) changes across the levels of another
#          (--group_var, e.g. tissue_type). Statistics use the sample
#          (--sample_var) as the unit of replication: speckle's propeller
#          transform of per-sample proportions, then limma on every pairwise
#          contrast of --group_var. When a sample contributes cells to more
#          than one group, limma::duplicateCorrelation() blocks on the sample.
#          No test is run on pooled cell counts (pseudoreplication).

pacman::p_load(optparse)
pdf(NULL)


##################### Options ########################
option_list <- list(
  make_option("--input", default = "", type = "character",
              help = "Seurat object (.rds or .qs2), a metadata data.frame saved as .rds or .qs2 (e.g. from seurat_save_metadata.R), or a metadata .tsv/.tsv.gz [required]"),
  make_option("--group_var", default = "", type = "character",
              help = "Metadata column whose levels are compared, e.g. tissue_type [required]"),
  make_option("--comp_var", default = "celltype", type = "character",
              help = "Metadata column whose composition is measured [default: %default]"),
  make_option("--sample_var", default = "patient", type = "character",
              help = "Metadata column that identifies the sample, the unit of replication [default: %default]"),
  make_option("--levels", default = "", type = "character",
              help = "Comma-separated level order for --group_var. Contrasts are all pairs in this order, earlier level is the numerator. Default: existing factor levels, else natural sort"),
  make_option("--min_cells", default = 0, type = "integer",
              help = "Drop sample x group combinations with fewer cells than this [default: %default]"),
  make_option("--transform", default = "logit", type = "character",
              help = "Transform of the proportions before limma: logit or asin [default: %default]"),
  make_option("--outdir", default = "plots", type = "character",
              help = "Output directory for PNGs [default: %default]"),
  make_option("--resultsdir", default = "results", type = "character",
              help = "Output directory for TSV and XLSX tables [default: %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

inputfile  = opts$input
group_var  = opts$group_var
comp_var   = opts$comp_var
sample_var = opts$sample_var
min_cells  = opts$min_cells
transform  = opts$transform
plotsdir   = opts$outdir
resultsdir = opts$resultsdir

if (inputfile == "" || !file.exists(inputfile)) {
  print_help(opt_parser)
  stop("--input must be an existing .rds, .qs2, .tsv or .tsv.gz file")
}
input_ext = tolower(tools::file_ext(sub("\\.gz$", "", inputfile, ignore.case = TRUE)))
if (!input_ext %in% c("rds", "qs2", "tsv")) {
  stop("--input extension must be .rds, .qs2, .tsv or .tsv.gz")
}
if (group_var == "") {
  print_help(opt_parser)
  stop("--group_var is required")
}
if (length(unique(c(group_var, comp_var, sample_var))) < 3) {
  stop("--group_var, --comp_var and --sample_var must be three different columns")
}
# Output tables name their columns after these three, so a metadata column
# with one of the script's own column names would collide.
reserved = c("n", "total", "prop", "transformed_prop", "observed", "expected",
             "comp_total", "group_total", "log2_enrichment", "contrast", "numerator",
             "denominator", "mean_prop_numerator", "mean_prop_denominator", "prop_ratio",
             "estimate", "t", "p_value", "fdr", "model", "transform", "Total_n")
clash = intersect(c(group_var, comp_var, sample_var), reserved)
if (length(clash) > 0) {
  stop("Column name(s) clash with output column names: ", paste(clash, collapse = ", "),
       ". Rename the column(s) in the metadata")
}
if (!transform %in% c("logit", "asin")) {
  stop("--transform must be logit or asin")
}
# optparse returns the raw string, with only a warning, when an integer option
# gets a non-number.
if (!is.numeric(min_cells) || is.na(min_cells) || min_cells < 0) {
  stop("--min_cells must be a whole number, 0 or more")
}

group_levels = trimws(strsplit(opts$levels, ",", fixed = TRUE)[[1]])
group_levels = group_levels[nzchar(group_levels)]
if (anyDuplicated(group_levels)) {
  stop("--levels has duplicate values")
}

message("\nArguments:")
print(opts)
message("")
#################################################


pacman::p_load(nvutils, sct2, speckle, limma, gtools, tidyverse)


# Slug a value for use in a file or directory name: labels like CAF/MSC exist.
slug = function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}


################## Read metadata ###################
# ReadSeurat() stops on anything that is not a Seurat object, so the .rds/.qs2
# branch reads the object directly and lets its class decide: a Seurat object
# gives its meta.data, a data.frame (seurat_save_metadata.R output) is used as is.
message("\nReading ", inputfile)
if (input_ext == "tsv") {
  metadata = read_tsv(inputfile, show_col_types = FALSE)
} else {
  if (input_ext == "qs2") {
    obj = qs2::qs_read(inputfile)
  } else {
    obj = readRDS(inputfile)
  }
  if (inherits(obj, "Seurat")) {
    metadata = obj[[]]
  } else if (is.data.frame(obj)) {
    metadata = obj
  } else {
    stop("--input holds an object of class ", paste(class(obj), collapse = ", "),
         "; expected a Seurat object or a data.frame")
  }
  rm(obj)
  invisible(gc())
}
message("Read ", nrow(metadata), " cells")

missing_cols = setdiff(c(group_var, comp_var, sample_var), colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Column(s) not found: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(colnames(metadata), collapse = ", "))
}

# Level order for comp_var and group_var: keep an existing factor's levels,
# else sort naturally. Keep only levels that occur in the data.
level_order = function(x) {
  if (is.factor(x)) {
    lv = levels(droplevels(x))
  } else {
    lv = gtools::mixedsort(unique(as.character(x[!is.na(x)])))
  }
  lv
}

md = tibble(sample = metadata[[sample_var]],
            group  = metadata[[group_var]],
            comp   = metadata[[comp_var]])
rm(metadata)

n_na = sum(!complete.cases(md))
if (n_na > 0) {
  message("Dropping ", n_na, " cells with NA in ", group_var, ", ", comp_var, " or ", sample_var)
  md = md[complete.cases(md), ]
}

# Quote values in messages so stray spaces in the data are visible.
quoted = function(x) {
  paste0('"', x, '"', collapse = ", ")
}

data_group_levels = level_order(md$group)
if (length(group_levels) == 0) {
  group_levels = data_group_levels
} else {
  not_in_levels = setdiff(data_group_levels, group_levels)
  if (length(not_in_levels) > 0) {
    stop(group_var, " value(s) not in --levels: ", quoted(not_in_levels))
  }
  not_in_data = setdiff(group_levels, data_group_levels)
  if (length(not_in_data) > 0) {
    stop("--levels value(s) not found in ", group_var, ": ", quoted(not_in_data),
         ". Values in the data: ", quoted(data_group_levels))
  }
}
if (length(group_levels) < 2) {
  stop(group_var, " needs at least 2 levels; found: ", paste(group_levels, collapse = ", "))
}

md = md %>%
  mutate(sample = factor(as.character(sample), levels = level_order(sample)),
         group = factor(as.character(group), levels = group_levels),
         comp  = factor(as.character(comp), levels = level_order(comp)))

message("\n", group_var, " levels: ", paste(group_levels, collapse = ", "))
message(comp_var, " levels: ", nlevels(md$comp))


################## Filter small sample x group combinations ###################
unit_sizes = md %>% count(sample, group, name = "total")
if (min_cells > 0) {
  dropped = unit_sizes %>% filter(total < min_cells)
  if (nrow(dropped) > 0) {
    message("\nDropping ", nrow(dropped), " ", sample_var, " x ", group_var,
            " combination(s) with fewer than ", min_cells, " cells:")
    print(as.data.frame(dropped))
    md = md %>% anti_join(dropped, by = c("sample", "group"))
    unit_sizes = unit_sizes %>% filter(total >= min_cells)
    # A comp level whose cells were all in dropped combinations would otherwise
    # stay as an all-zero row in every table and in the limma fit.
    lost_levels = setdiff(levels(md$comp), unique(as.character(md$comp)))
    if (length(lost_levels) > 0) {
      message("Removing ", comp_var, " level(s) with no cells left: ", quoted(lost_levels))
      md = md %>% mutate(comp = droplevels(comp))
    }
  }
}
if (nlevels(md$comp) < 2) {
  stop(comp_var, " needs at least 2 levels; found: ", quoted(levels(md$comp)))
}

units_per_group = unit_sizes %>% count(group, name = "n_samples")
message("\n", sample_var, "s per ", group_var, ":")
print(as.data.frame(units_per_group))
single_groups = as.character(units_per_group$group[units_per_group$n_samples == 1])
if (length(single_groups) > 0) {
  message("WARNING: only one ", sample_var, " in ", group_var, " level(s) ", quoted(single_groups),
          ". Contrasts with these levels rest on one sample; interpret with care")
}
empty_groups = setdiff(group_levels, as.character(units_per_group$group))
if (length(empty_groups) > 0) {
  stop("No samples left in ", group_var, " level(s): ", paste(empty_groups, collapse = ", "),
       " (check --min_cells)")
}
if (nrow(unit_sizes) <= length(group_levels)) {
  stop("Need more ", sample_var, " x ", group_var, " combinations (", nrow(unit_sizes),
       ") than ", group_var, " levels (", length(group_levels), ") to estimate variance")
}

outname = paste0("composition_", slug(group_var), "_", slug(comp_var))
dir.create(plotsdir, recursive = TRUE, showWarnings = FALSE)
dir.create(resultsdir, recursive = TRUE, showWarnings = FALSE)

write_table = function(df, name) {
  path = file.path(resultsdir, paste0(outname, "_", name))
  write_tsv(df, paste0(path, ".tsv"))
  write_xlsx_pretty(df, paste0(path, ".xlsx"))
  message("Wrote ", path, ".tsv/.xlsx")
}


################## Table 1: pooled counts and percentages ###################
# One wide table: per group, the cell count (n), the share of the group's cells
# (pct_of_group) and the share of the comp level's cells (pct_of_comp). A
# "Total" row gives group sizes; a Total_n column gives comp level sizes.
message("\nPooled composition")
pooled_long = md %>%
  count(comp, group, name = "n", .drop = FALSE) %>%
  group_by(group) %>% mutate(pct_of_group = n / sum(n)) %>%
  group_by(comp) %>% mutate(pct_of_comp = n / sum(n)) %>%
  ungroup() %>%
  mutate(comp = as.character(comp))
pooled_totals = md %>%
  count(group, name = "n", .drop = FALSE) %>%
  mutate(comp = "Total", pct_of_group = 1, pct_of_comp = n / sum(n))
pooled = bind_rows(pooled_long, pooled_totals) %>%
  pivot_wider(names_from = group, values_from = c(n, pct_of_group, pct_of_comp),
              names_glue = "{group}_{.value}", names_vary = "slowest") %>%
  left_join(bind_rows(pooled_long, pooled_totals) %>%
              group_by(comp) %>% summarise(Total_n = sum(n)),
            by = "comp")
write_table(pooled %>% rename(!!comp_var := comp), "pooled")


################## Table 2: per-sample proportions ###################
# One row per sample x group x comp level. complete() fills comp levels absent
# from a sample with n = 0 so they count as 0, not missing.
per_sample = md %>%
  count(sample, group, comp, name = "n", .drop = TRUE) %>%
  complete(nesting(sample, group), comp, fill = list(n = 0L)) %>%
  left_join(unit_sizes, by = c("sample", "group")) %>%
  mutate(prop = n / total) %>%
  arrange(group, sample, comp)
# Written after the statistics, once the transformed proportions exist.


################## Table 3: pooled enrichment ###################
# log2(observed / expected) from the pooled table. Descriptive only; -Inf when
# a comp level has no cells in a group.
n_total = nrow(md)
enrich = md %>%
  count(comp, group, name = "observed", .drop = FALSE) %>%
  group_by(comp) %>% mutate(comp_total = sum(observed)) %>%
  group_by(group) %>% mutate(group_total = sum(observed)) %>%
  ungroup() %>%
  mutate(expected = as.numeric(comp_total) * group_total / n_total,  # integer product overflows past ~46k x 46k
         log2_enrichment = log2(observed / expected))
write_table(enrich %>% rename(!!comp_var := comp, !!group_var := group), "enrichment")


################## Table 4: propeller statistics ###################
message("\nStatistics: ", transform, " transform, limma")
md = md %>% mutate(unit = paste(sample, group, sep = "__"))
# The key must be one-to-one with sample x group; a name containing "__" could
# make two combinations share a key.
if (n_distinct(md$unit) != nrow(unit_sizes)) {
  stop(sample_var, " x ", group_var, " keys collide; a ", sample_var, " or ", group_var,
       " value contains \"__\"")
}
props = getTransformedProps(clusters = md$comp, sample = md$unit, transform = transform)

unit_info = unit_sizes %>%
  mutate(unit = paste(sample, group, sep = "__")) %>%
  slice(match(colnames(props$TransformedProps), unit))
stopifnot(identical(unit_info$unit, colnames(props$TransformedProps)))

# The per-sample table carries the exact values limma fits (transformed_prop):
# one limma "gene" per comp level, one column per sample x group, design on
# group, block on sample. The logit transform adds 0.5 to every count, so
# transformed_prop is not exactly logit(prop).
tp = props$TransformedProps
per_sample = per_sample %>%
  mutate(unit = paste(sample, group, sep = "__"), comp_chr = as.character(comp)) %>%
  left_join(tibble(unit = rep(colnames(tp), each = nrow(tp)),
                   comp_chr = rep(rownames(tp), times = ncol(tp)),
                   transformed_prop = as.vector(tp)),
            by = c("unit", "comp_chr")) %>%
  select(-unit, -comp_chr)
stopifnot(!anyNA(per_sample$transformed_prop))
write_table(per_sample %>% rename(!!sample_var := sample, !!group_var := group, !!comp_var := comp),
            "per_sample")

# Syntactic names for the design so makeContrasts() accepts levels like
# "Tumor core"; the original labels go back into the output table.
design_names = make.names(group_levels, unique = TRUE)
design = model.matrix(~ 0 + unit_info$group)
colnames(design) = design_names

pairs = combn(seq_along(group_levels), 2)
contrast_labels = paste0(group_levels[pairs[1, ]], "_vs_", group_levels[pairs[2, ]])
contrast_defs = paste(design_names[pairs[1, ]], "-", design_names[pairs[2, ]])
contrasts = makeContrasts(contrasts = contrast_defs, levels = design)
colnames(contrasts) = contrast_labels

repeated = any(duplicated(unit_info$sample))
model_name = "limma"
if (repeated) {
  model_name = "limma_dupcor"
  n_shared = sum(table(unit_info$sample) > 1)
  message(n_shared, " ", sample_var, "(s) contribute to more than one ", group_var,
          " level: fitting with duplicateCorrelation(block = ", sample_var, ")")
  dupcor = duplicateCorrelation(props$TransformedProps, design, block = unit_info$sample)
  message("Consensus within-", sample_var, " correlation: ", round(dupcor$consensus.correlation, 3))
  fit = lmFit(props$TransformedProps, design, block = unit_info$sample,
              correlation = dupcor$consensus.correlation)
} else {
  message("Each ", sample_var, " is in one ", group_var, " level: fitting without blocking")
  fit = lmFit(props$TransformedProps, design)
}
fit2 = eBayes(contrasts.fit(fit, contrasts), robust = TRUE)

# Mean per-sample proportion per group, for the proportion ratio.
mean_props = per_sample %>%
  group_by(comp, group) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

stats = map_dfr(seq_along(contrast_labels), function(i) {
  num = group_levels[pairs[1, i]]
  den = group_levels[pairs[2, i]]
  topTable(fit2, coef = i, number = Inf, sort.by = "none") %>%
    rownames_to_column("comp") %>%
    transmute(contrast = contrast_labels[i], numerator = num, denominator = den,
              comp, estimate = logFC, t, p_value = P.Value, fdr = adj.P.Val) %>%
    left_join(mean_props %>% filter(group == num) %>% select(comp, mean_prop_numerator = mean_prop),
              by = "comp") %>%
    left_join(mean_props %>% filter(group == den) %>% select(comp, mean_prop_denominator = mean_prop),
              by = "comp") %>%
    mutate(prop_ratio = mean_prop_numerator / mean_prop_denominator)
}) %>%
  mutate(comp = as.character(comp)) %>%
  relocate(estimate, t, p_value, fdr, .after = last_col()) %>%
  mutate(model = model_name, transform = transform)
write_table(stats %>% rename(!!comp_var := comp), "stats")

message("\nComp levels with FDR < 0.05:")
print(as.data.frame(stats %>% filter(fdr < 0.05) %>% select(contrast, comp, prop_ratio, fdr)))


################## Plots ###################
save_plot = function(p, name, width, height) {
  path = file.path(plotsdir, paste0(outname, "_", name, ".png"))
  ggsave(path, plot = p, width = width, height = height, bg = "white")
  message("Wrote ", path)
}

n_comp = nlevels(md$comp)
n_units = nrow(unit_sizes)

# Pooled stacked bar per group
p = two_category_barplot(md, category = "group", subcategory = "comp",
                         title = paste(comp_var, "composition by", group_var),
                         legend_title = comp_var) +
  labs(x = group_var)
save_plot(p, "pooled_stackedbar", width = 3 + 0.8 * length(group_levels), height = 7)

# Stacked bar per sample, faceted by group
p = two_category_barplot(md, category = "sample", subcategory = "comp",
                         title = paste(comp_var, "composition per", sample_var),
                         legend_title = comp_var) +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
  labs(x = sample_var)
save_plot(p, "per_sample_stackedbar", width = 4 + 0.3 * n_units, height = 7)

# Same bars faceted by sample, so one sample's groups sit side by side
p = two_category_barplot(md, category = "group", subcategory = "comp",
                         title = paste(comp_var, "composition per", sample_var, "and", group_var),
                         legend_title = comp_var) +
  facet_grid(cols = vars(sample), scales = "free_x", space = "free_x") +
  labs(x = group_var)
save_plot(p, "per_sample_stackedbar_by_sample", width = 4 + 0.3 * n_units, height = 7)

# Per-comp-level boxplots of sample-level proportions; points of one sample
# joined across groups.
ncol_box = ceiling(sqrt(n_comp))
p = ggplot(per_sample, aes(x = group, y = prop)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_line(aes(group = sample), color = "grey60", linewidth = 0.3) +
  geom_point(size = 1) +
  facet_wrap(vars(comp), scales = "free_y", ncol = ncol_box) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = group_var, y = paste("Proportion of", sample_var, "cells"),
       title = paste(comp_var, "proportion per", sample_var, "by", group_var)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_plot(p, "per_sample_boxplot",
          width = max(7, 1 + ncol_box * (0.6 + 0.4 * length(group_levels))),
          height = 1 + ceiling(n_comp / ncol_box) * 2.2)

# Heatmap of pooled log2 enrichment; -Inf (no cells) drawn as NA.
# The 0.1 floor keeps the color scale valid when every value is 0 or non-finite.
limit = max(c(0.1, abs(enrich$log2_enrichment[is.finite(enrich$log2_enrichment)])))
p = enrich %>%
  mutate(log2_enrichment = if_else(is.finite(log2_enrichment), log2_enrichment, NA_real_),
         comp = fct_rev(comp)) %>%
  ggplot(aes(x = group, y = comp, fill = log2_enrichment)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(log2_enrichment, 1)), size = 3, na.rm = TRUE) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "indianred",
                       limits = c(-limit, limit), na.value = "grey80",
                       name = "log2(obs/exp)") +
  labs(x = group_var, y = comp_var, title = paste(comp_var, "enrichment by", group_var)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_plot(p, "enrichment_heatmap", width = 3 + 0.8 * length(group_levels), height = 2 + 0.25 * n_comp)


message("\nDone.")


cat("\n\n")
devtools::session_info()
