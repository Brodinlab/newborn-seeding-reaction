#!/usr/bin/env Rscript
## Bulk RNA-seq (S3A / Fig4D / Fig4E / Fig6D)
## Input: input/bulkRNAseq/  ->  Output: output/bulkRNAseq/

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(DESeq2)
  library(stringr)
  library(purrr)
  library(sva)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")

input_dir <- file.path(repo_root, "input", "bulkRNAseq")
output_dir <- file.path(repo_root, "output", "bulkRNAseq")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path_meta <- file.path(input_dir, "merged_meta_all.csv")
path_counts <- file.path(input_dir, "merged_data_all_raw_counts.csv.gz")
path_tcr_meta <- file.path(input_dir, "tcr_sample_metadata.csv")
path_trust4_dir <- file.path(input_dir, "trust4")

ifng_score <- c(
  "CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "IRF1",
  "GBP1", "GBP2", "GBP4", "GBP5", "TAP1", "TAP2", "PSMB8", "PSMB9"
)
ifn1_isg <- c(
  "OAS1", "OAS2", "OAS3", "OASL", "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3",
  "ISG15", "RSAD2", "IFI6", "IFI27", "IFI44", "IFI44L", "LY6E", "HERC5",
  "CMPK2", "EPSTI1", "STAT1", "STAT2", "IRF7", "IRF9", "USP18"
)
nfkb_score <- c(
  "IL6", "TNF", "IL1B", "CCL2", "CCL3", "CCL4", "CXCL8", "CXCL2",
  "ICAM1", "VCAM1", "PTGS2"
)

gini_coef <- function(x) {
  x <- sort(x[x > 0])
  n <- length(x)
  if (n < 2) return(NA_real_)
  2 * sum((seq_len(n)) * x) / (n * sum(x)) - (n + 1) / n
}

validate_inputs <- function() {
  missing <- c(path_meta, path_counts)[!file.exists(c(path_meta, path_counts))]
  if (length(missing)) {
    stop("Missing required bulk RNA-seq inputs in input/bulkRNAseq/: ",
         paste(basename(missing), collapse = ", "))
  }
}

read_trust4_report <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  lines <- readLines(filepath, warn = FALSE)
  if (length(lines) < 2) return(NULL)
  cols <- str_split(str_remove(lines[1], "^#\\s*"), "\t")[[1]]
  dt <- read.delim(text = paste(lines[-1], collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  colnames(dt) <- cols[seq_len(ncol(dt))]
  dt
}

cat("=== bulk_mrnaseq ===\n")
validate_inputs()

meta_out <- read_csv(path_meta, show_col_types = FALSE)
counts_df <- read_csv(path_counts, show_col_types = FALSE)
counts_mat <- as.matrix(counts_df[, -1])
rownames(counts_mat) <- counts_df$gene
storage.mode(counts_mat) <- "integer"

meta_final <- as.data.frame(meta_out)
rownames(meta_final) <- meta_final$sample_id
counts_mat <- counts_mat[, meta_final$sample_id, drop = FALSE]
stopifnot(identical(colnames(counts_mat), rownames(meta_final)))
cat("Loaded:", nrow(meta_out), "samples;", nrow(counts_df), "genes\n")

plot_s3a_cxcl10 <- function() {
  if (!"CXCL10" %in% rownames(counts_mat)) {
    stop("CXCL10 not found in count matrix for S3A.")
  }
  meta_cb <- meta_final
  meta_cb$sample_id <- rownames(meta_cb)
  meta_cb <- meta_cb[!duplicated(meta_cb$sample_id), , drop = FALSE]
  rownames(meta_cb) <- meta_cb$sample_id
  counts_cb <- counts_mat[, meta_cb$sample_id, drop = FALSE]
  stopifnot(identical(colnames(counts_cb), rownames(meta_cb)))

  batch_corrected <- as.data.frame(
    ComBat_seq(counts = as.matrix(counts_cb), batch = meta_cb$batch)
  )

  meta_1w <- meta_cb %>%
    dplyr::filter(!is.na(time), time > 0, time <= 7)

  batch_corrected_1w <- batch_corrected[, meta_1w$sample_id, drop = FALSE]
  df_cxcl10 <- data.frame(
    Sample = colnames(batch_corrected_1w),
    CXCL10 = as.numeric(batch_corrected_1w["CXCL10", ]),
    stringsAsFactors = FALSE
  )
  meta_subset <- data.frame(
    Sample = rownames(meta_cb),
    IFNG_category = meta_cb$IFNG_category,
    stringsAsFactors = FALSE
  )
  df_cxcl10 <- df_cxcl10 %>%
    dplyr::left_join(meta_subset, by = "Sample") %>%
    dplyr::mutate(
      IFNG_category = ifelse(is.na(IFNG_category), "Other", IFNG_category),
      Group = ifelse(IFNG_category == "Other", "No information", IFNG_category)
    ) %>%
    dplyr::arrange(CXCL10) %>%
    dplyr::mutate(Index = dplyr::row_number())

  median_value <- median(df_cxcl10$CXCL10, na.rm = TRUE)
  median_position <- which.min(abs(df_cxcl10$CXCL10 - median_value))
  q1_value <- quantile(df_cxcl10$CXCL10, 0.25, na.rm = TRUE)
  q3_value <- quantile(df_cxcl10$CXCL10, 0.75, na.rm = TRUE)
  q4_value <- quantile(df_cxcl10$CXCL10, 1.0, na.rm = TRUE)
  q1_position <- which.min(abs(df_cxcl10$CXCL10 - q1_value))
  q3_position <- which.min(abs(df_cxcl10$CXCL10 - q3_value))
  q4_position <- which.min(abs(df_cxcl10$CXCL10 - q4_value))

  unique_groups <- unique(df_cxcl10$Group)
  n_groups <- length(unique_groups)
  colors <- if (n_groups <= 3) {
    c("steelblue", "red", "forestgreen", "purple", "orange")
  } else {
    rainbow(n_groups, alpha = 0.8)
  }
  colors <- colors[seq_len(n_groups)]
  names(colors) <- unique_groups
  if ("No information" %in% names(colors)) {
    colors["No information"] <- "lightgray"
  }

  y_max <- max(df_cxcl10$CXCL10, na.rm = TRUE)
  n_samp <- nrow(df_cxcl10)

  p <- ggplot(df_cxcl10, aes(x = Index, y = CXCL10 + 1)) +
    geom_col(aes(fill = Group), color = "white", linewidth = 0.2, alpha = 0.8) +
    scale_fill_manual(values = colors, name = "IFNG Category") +
    geom_hline(yintercept = median_value + 1, color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = q1_value + 1, color = "blue", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = q3_value + 1, color = "blue", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = q4_value + 1, color = "darkgreen", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    geom_vline(xintercept = median_position, color = "red", linetype = "dotted", alpha = 0.7) +
    annotate("text",
      x = median_position + n_samp * 0.15,
      y = median_value + 1 + y_max * 0.05,
      label = paste("Median:", round(median_value, 2)),
      color = "red", fontface = "bold", size = 3.5, hjust = 0.5, vjust = 0.5) +
    annotate("text",
      x = q1_position - n_samp * 0.08,
      y = q1_value + 1 + y_max * 0.03,
      label = paste("Q1:", round(q1_value, 2)),
      color = "blue", fontface = "bold", size = 3, hjust = 0.5, vjust = 0.5) +
    annotate("text",
      x = q3_position + n_samp * 0.08,
      y = q3_value + 1 + y_max * 0.03,
      label = paste("Q3:", round(q3_value, 2)),
      color = "blue", fontface = "bold", size = 3, hjust = 0.5, vjust = 0.5) +
    annotate("text",
      x = q4_position - n_samp * 0.08,
      y = q4_value + 1 - y_max * 0.03,
      label = paste("Q4 (Max):", round(q4_value, 2)),
      color = "darkgreen", fontface = "bold", size = 3, hjust = 0.5, vjust = 0.5) +
    annotate("segment",
      x = median_position + n_samp * 0.1,
      y = median_value + 1 + y_max * 0.02,
      xend = median_position, yend = median_value + 1,
      arrow = arrow(length = grid::unit(0.2, "cm")),
      color = "red", linewidth = 0.6) +
    labs(
      title = "CXCL10 Gene Expression Level Distribution (1W)",
      x = "Samples (ordered by CXCL10 expression)",
      y = "CXCL10 Expression Level"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    ) +
    scale_x_continuous(
      breaks = seq(1, n_samp, by = max(1, n_samp %/% 10)),
      labels = function(x) df_cxcl10$Sample[x]
    )

  out_name <- "S3A_cxcl10_threshold_compared_with_ifng_67samples.pdf"
  ggsave(file.path(output_dir, out_name), p, width = 9, height = 6)
  message("Wrote S3A (", n_samp, " 1W samples)")
}

plot_s3a_cxcl10()
cat("Wrote S3A\n")

plot_signature_child_lines <- function(gene_list, out_name) {
  nm <- deparse(substitute(gene_list))
  parts <- strsplit(nm, "_")[[1]]
  cap <- paste0(toupper(parts[1]), "_", paste(parts[-1], collapse = "_"))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta_final,
    design = ~1
  )
  vst_mat <- assay(vst(dds, blind = TRUE))
  avail <- intersect(gene_list, rownames(vst_mat))
  z <- t(scale(t(vst_mat[avail, , drop = FALSE])))
  score <- colMeans(z, na.rm = TRUE)

  pm <- meta_final %>%
    dplyr::mutate(
      signature_score = score[match(sample_id, names(score))],
      IFNG_label = ifelse(IFNG_category == "IFNgHIGH", "High", "Low"),
      time_plot = factor(
        time_group,
        levels = c("CB", "1W", "3_4M", "5_6M"),
        labels = c("CB", "1W", "3-4M", "5-6M")
      )
    ) %>%
    dplyr::filter(!is.na(signature_score), !is.na(IFNG_label), !is.na(time_plot))

  child_traj <- pm %>%
    dplyr::group_by(study_id, IFNG_label, time_plot) %>%
    dplyr::summarise(signature_score = median(signature_score, na.rm = TRUE), .groups = "drop")

  med_traj <- child_traj %>%
    dplyr::group_by(IFNG_label, time_plot) %>%
    dplyr::summarise(median_score = median(signature_score, na.rm = TRUE), .groups = "drop")

  p <- ggplot() +
    geom_line(
      data = child_traj,
      aes(x = time_plot, y = signature_score, group = interaction(IFNG_label, study_id), color = IFNG_label),
      alpha = 0.22, linewidth = 0.45
    ) +
    geom_jitter(
      data = child_traj, aes(x = time_plot, y = signature_score, color = IFNG_label),
      width = 0.12, height = 0, alpha = 0.28, size = 1.2
    ) +
    geom_line(
      data = med_traj,
      aes(x = time_plot, y = median_score, group = IFNG_label, color = IFNG_label),
      linewidth = 1.6
    ) +
    geom_point(data = med_traj, aes(x = time_plot, y = median_score, color = IFNG_label), size = 2.6) +
    scale_color_manual(name = "IFNg", values = c("High" = "#FCBA6F", "Low" = "#5E4C5F")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "Time", y = "Gene Score", title = paste(cap, "Gene Panel expression Over time"))

  ggsave(file.path(output_dir, out_name), p, width = 9, height = 6)
}

plot_signature_child_lines(ifng_score, "Fig4D_IFNG_score_child_lines_median_signature_score_over_time.pdf")
plot_signature_child_lines(ifn1_isg, "Fig4D_IFN1_isg_child_lines_median_signature_score_over_time.pdf")
plot_signature_child_lines(nfkb_score, "Fig4D_NFKB_score_child_lines_median_signature_score_over_time.pdf")
cat("Wrote Fig4D\n")

plot_fig4e_ridgelines <- function() {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Install ggridges: install.packages(\"ggridges\")")
  }
  library(ggridges)

  meta_tbl <- meta_final %>%
    dplyr::mutate(tnum = as.numeric(time)) %>%
    dplyr::filter(!is.na(tnum), !is.na(batch), IFNG_category %in% c("IFNgLOW", "IFNgHIGH"))

  meta_nocb <- meta_tbl %>% dplyr::filter(tnum > 0)
  counts_nocb <- counts_mat[, meta_nocb$sample_id, drop = FALSE]
  rownames(meta_nocb) <- meta_nocb$sample_id
  meta_nocb$batch <- factor(meta_nocb$batch)
  meta_nocb$IFNG_category <- factor(meta_nocb$IFNG_category, levels = c("IFNgLOW", "IFNgHIGH"))

  design_f <- if (nlevels(meta_nocb$batch) < 2) ~1 else ~batch
  dds <- DESeqDataSetFromMatrix(round(counts_nocb), meta_nocb, design_f)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, quiet = TRUE)
  vst_nocb <- assay(vst(dds, blind = FALSE))

  t_x <- log2(meta_nocb[colnames(vst_nocb), "tnum"] + 1)
  r2 <- apply(vst_nocb, 1, function(z) {
    r <- cor(as.numeric(z), t_x, use = "pairwise.complete.obs")
    if (!is.finite(r)) 0 else r^2
  })
  hvg <- names(sort(r2, decreasing = TRUE))[seq_len(min(2000L, length(r2)))]
  pca <- prcomp(t(vst_nocb[hvg, , drop = FALSE]), center = TRUE, scale. = FALSE)
  sgn <- sign(cor(pca$x[, 1], t_x, use = "pairwise.complete.obs"))
  if (is.na(sgn) || sgn == 0) sgn <- 1
  pseu <- pca$x[, 1] * sgn

  df_st <- meta_nocb %>%
    dplyr::mutate(pseu = pseu[match(sample_id, names(pseu))]) %>%
    dplyr::filter(tnum == 7 | tnum >= 90) %>%
    dplyr::mutate(
      time_split4 = factor(
        dplyr::case_when(
          tnum == 7 ~ "1W_7d",
          tnum %in% c(90, 120) ~ "3_4M_90_120d",
          tnum %in% c(150, 180) ~ "5_6M_150_180d",
          TRUE ~ NA_character_
        ),
        levels = c("1W_7d", "3_4M_90_120d", "5_6M_150_180d")
      )
    ) %>%
    dplyr::filter(!is.na(time_split4))

  ref <- df_st %>%
    dplyr::filter(time_split4 == "1W_7d") %>%
    dplyr::group_by(IFNG_category) %>%
    dplyr::summarise(ref = median(pseu, na.rm = TRUE), .groups = "drop")

  fig5e_col_low <- "#67596D"
  fig5e_col_high <- "#F8C48C"
  fig5e_xlim_max <- 100

  ridge_plot <- function(d, title, fname) {
    d <- d %>%
      dplyr::left_join(ref, by = "IFNG_category") %>%
      dplyr::mutate(delta = pseu - ref)
    p <- ggplot(d, aes(x = delta, y = IFNG_category, fill = IFNG_category)) +
      ggridges::geom_density_ridges(alpha = 0.88, scale = 0.95) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
      scale_fill_manual(values = c("IFNgHIGH" = fig5e_col_high, "IFNgLOW" = fig5e_col_low)) +
      coord_cartesian(xlim = c(NA, fig5e_xlim_max), expand = FALSE) +
      ggridges::theme_ridges(center_axis_labels = TRUE) +
      theme(legend.position = "none") +
      labs(title = title, x = "Delta pseudotime (signed PC1 - within-IFN ref 1W)", y = NULL)
    ggsave(file.path(output_dir, fname), p, width = 7, height = 3.6)
  }

  ridge_plot(
    dplyr::filter(df_st, time_split4 == "3_4M_90_120d"),
    "3-4M (90-120 d): pseu change vs within-IFN median 1W_7d",
    "Fig4E_plot_delta_pseu_ridgeline_3_4M_vs_withinIFN_median1W.pdf"
  )
  ridge_plot(
    dplyr::filter(df_st, time_split4 == "5_6M_150_180d"),
    "5-6M (150-180 d): pseu change vs within-IFN median 1W_7d",
    "Fig4E_plot_delta_pseu_ridgeline_5_6M_vs_withinIFN_median1W.pdf"
  )
}
plot_fig4e_ridgelines()
cat("Wrote Fig4E\n")

plot_fig6d_tcr <- function() {
  if (!file.exists(path_tcr_meta) || !dir.exists(path_trust4_dir)) {
    message("Fig6D skipped: add tcr_sample_metadata.csv and trust4/*_report.tsv under input/bulkRNAseq/")
    return(invisible(NULL))
  }
  reps <- list.files(path_trust4_dir, "_report.tsv$", full.names = TRUE)
  if (length(reps) == 0) {
    message("Fig6D skipped: no trust4/*_report.tsv files found")
    return(invisible(NULL))
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install patchwork: install.packages(\"patchwork\")")
  }
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Install ggridges: install.packages(\"ggridges\")")
  }
  library(patchwork)
  library(ggridges)

  col_low <- "#5E4C5F"
  col_high <- "#FCBA6F"

  stars_from_fdr <- function(fdr) {
    if (is.na(fdr) || fdr >= 0.05) return("")
    if (fdr < 0.001) return("***")
    if (fdr < 0.01) return("**")
    if (fdr < 0.05) return("*")
    ""
  }

  meta_t <- read.csv(path_tcr_meta, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      IFNG_category = dplyr::case_when(
        IFNG_category == "High" ~ "IFNgHIGH",
        IFNG_category == "Low" ~ "IFNgLOW",
        TRUE ~ as.character(IFNG_category)
      ),
      IFN_group = factor(
        ifelse(IFNG_category == "IFNgHIGH", "High", "Low"),
        levels = c("Low", "High")
      ),
      batch = str_replace(sample_id, "_.*$", "")
    ) %>%
    dplyr::filter(IFNG_category %in% c("IFNgHIGH", "IFNgLOW"))

  all_tcr <- purrr::map_dfr(reps, function(f) {
    sid <- sub("_report\\.tsv$", "", basename(f))
    dt <- read_trust4_report(f)
    if (is.null(dt)) return(NULL)
    dt$sample_id <- sid
    dt
  })

  tcr <- all_tcr %>%
    dplyr::filter(!str_detect(V, "^IGH|^IGK|^IGL")) %>%
    dplyr::mutate(count = as.numeric(count)) %>%
    dplyr::filter(count >= 2)

  div <- tcr %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      total_clones = sum(count),
      shannon = -sum((count / sum(count)) * log(count / sum(count) + 1e-12)),
      simpson = 1 - sum((count / sum(count))^2),
      clonality = if (dplyr::n() > 1) 1 - shannon / log(dplyr::n()) else NA_real_,
      gini = gini_coef(count),
      .groups = "drop"
    )

  dm <- div %>%
    dplyr::inner_join(meta_t, by = "sample_id") %>%
    dplyr::mutate(log10_depth = log10(total_clones + 1))

  fig6d_metrics <- c("clonality", "gini", "simpson", "shannon")
  fig6d_titles <- c(
    clonality = "TCR clonality",
    gini = "TCR Gini coeff.",
    simpson = "TCR, Simpson diversity",
    shannon = "TCR Shannon diversity"
  )

  reg_rows <- purrr::map_dfr(fig6d_metrics, function(m) {
    fit <- tryCatch(
      lm(as.formula(paste(m, "~ IFN_group + log10_depth + batch")), data = dm),
      error = function(e) NULL
    )
    p_group <- NA_real_
    if (!is.null(fit)) {
      s <- summary(fit)$coefficients
      if ("IFN_groupHigh" %in% rownames(s)) {
        p_group <- s["IFN_groupHigh", "Pr(>|t|)"]
      }
    }
    tibble::tibble(metric = m, p_group = p_group)
  })
  reg_rows$fdr_group <- p.adjust(reg_rows$p_group, method = "fdr")
  fdr_lookup <- setNames(reg_rows$fdr_group, reg_rows$metric)

  ridge_theme <- theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      panel.grid = element_blank(),
      plot.margin = margin(4, 6, 4, 4)
    )

  plots <- purrr::map(fig6d_metrics, function(m) {
    fit <- tryCatch(
      lm(as.formula(paste(m, "~ IFN_group + log10_depth + batch")), data = dm),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    valid <- dm %>% dplyr::filter(is.finite(.data[[m]]))
    valid$predicted <- stats::predict(fit, newdata = valid)
    if (nrow(valid) < 4) return(NULL)

    med_df <- valid %>%
      dplyr::group_by(IFN_group) %>%
      dplyr::summarise(med = median(predicted, na.rm = TRUE), .groups = "drop")

    star <- stars_from_fdr(fdr_lookup[[m]])
    ggplot(valid, aes(x = predicted, y = IFN_group, fill = IFN_group, color = IFN_group)) +
      geom_density_ridges(
        alpha = 0.65,
        scale = 0.92,
        rel_min_height = 0.01,
        color = NA
      ) +
      geom_density_ridges(
        alpha = 0.65,
        scale = 0.92,
        rel_min_height = 0.01,
        fill = NA,
        linewidth = 0.25
      ) +
      geom_vline(
        data = med_df,
        aes(xintercept = med, color = IFN_group),
        linetype = "dashed",
        linewidth = 0.55
      ) +
      scale_fill_manual(values = c("Low" = col_low, "High" = col_high)) +
      scale_color_manual(values = c("Low" = col_low, "High" = col_high)) +
      scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
      labs(
        title = paste0(fig6d_titles[[m]], " ", star),
        x = "Model-predicted value",
        y = NULL
      ) +
      ridge_theme
  })
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0) {
    message("Fig6D skipped: no TCR panels could be built")
    return(invisible(NULL))
  }

  pdf(
    file.path(output_dir, "Fig6D_tcr_repertoire_diversity_cleaned_significant_predicted.pdf"),
    width = 12,
    height = 3.2
  )
  print(patchwork::wrap_plots(plots, ncol = 4, nrow = 1))
  dev.off()
}
plot_fig6d_tcr()
cat("Wrote Fig6D\n")
cat("Done. output -> ", output_dir, "\n")
