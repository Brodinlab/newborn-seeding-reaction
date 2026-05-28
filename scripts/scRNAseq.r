#!/usr/bin/env Rscript
## Xpress scRNA-seq (Experiment 2) -> output/scRNAseq (Fig7C–7E)
##
## Inputs (on GitHub):  input/scRNAseq/{barcode_annotation.txt, plate_metadata.csv, DATA_SOURCES.md}
## Outputs (on GitHub): output/scRNAseq/*.pdf
##
## output/scRNAseq/cache/ — local only (gitignored). Speeds up re-runs on your machine after
## a full rebuild; not required for publication and not on GitHub.
##
## Default:  Rscript scripts/scRNAseq.r
##   Uses cache if present; otherwise see message (use output PDFs or REBUILD_FROM_RAW=TRUE).
##
## Rebuild:  REBUILD_FROM_RAW=TRUE Rscript scripts/scRNAseq.r
##   Builds cache from UMI (DATA_SOURCES.md); first run ~30–60 min.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggridges)
  library(Matrix)
})

parse_rebuild_flag <- function() {
  tr <- commandArgs(trailingOnly = TRUE)
  identical(Sys.getenv("REBUILD_FROM_RAW"), "TRUE") || "--rebuild-from-raw" %in% tr
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
scrna_input <- file.path(repo_root, "input", "scRNAseq")
# Rebuild-only intermediates (gitignored); default run does not use this folder.
cache_dir <- file.path(repo_root, "output", "scRNAseq", "cache")
output_dir <- file.path(repo_root, "output", "scRNAseq")
tables_dir <- file.path(output_dir, "tables")
dir.create(scrna_input, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

REBUILD_FROM_RAW <- parse_rebuild_flag()

path_data_sources <- file.path(scrna_input, "DATA_SOURCES.md")
path_barcode <- file.path(scrna_input, "barcode_annotation.txt")
path_plate <- file.path(scrna_input, "plate_metadata.csv")
path_built <- file.path(cache_dir, "xpress_seurat_object_experiment2_only.rds")
path_processed <- file.path(cache_dir, "xpress_seurat_processed_experiment2_only.rds")
path_gene_cache <- file.path(cache_dir, "xpress_gene_mapping_cache.rds")
path_hvg_cache <- file.path(cache_dir, "xpress_hvg_biotype_cache.rds")

CONDITIONS <- c("NPC", "LPS (E.coli)", "LPS (Bacteroides)")
COND_COL <- c(
  "NPC" = "#4dac26",
  "LPS (E.coli)" = "#d01c8b",
  "LPS (Bacteroides)" = "#0571b0"
)

cluster_labels <- c(
  "0" = "C0 B-lineage-primed LMPP-like",
  "1" = "C1 qHSC",
  "2" = "C2 Stress",
  "3" = "C3 Erythroid progenitor-like",
  "4" = "C4 MPP-like",
  "5" = "C5 GMP",
  "6" = "C6 Basophil/Mast progenitor-like"
)
cluster_cols <- c(
  "C0 B-lineage-primed LMPP-like" = "#F28E2B",
  "C1 qHSC" = "#4E79A7",
  "C2 Stress" = "#59A14F",
  "C3 Erythroid progenitor-like" = "#EDC948",
  "C4 MPP-like" = "#2F5597",
  "C5 GMP" = "#E15759",
  "C6 Basophil/Mast progenitor-like" = "#B07AA1"
)
cluster_labels_pub <- c(
  "0" = "B-primed LMPP", "1" = "qHSC", "2" = "Stress",
  "3" = "Erythroid prog", "4" = "MPP-like", "5" = "GMP", "6" = "Baso/Mast prog"
)

gene_sets <- list(
  HSC_primitive = c("HLF", "MECOM", "AVP"),
  HSC_prog = c("HOXA9", "MEIS1", "CD34"),
  MPP_core = c("KIT", "SPINK2", "GATA2", "CRHBP", "SOX4", "FLT3"),
  LMPP = c("IKZF1", "BCL11A", "TCF3", "IL7R", "TCF7", "GATA3"),
  CLP_early = c("RAG1", "RAG2"),
  B_commitment = c("EBF1", "PAX5"),
  preB_like = c("VPREB1", "IGLL1"),
  T_early = c("BCL11B", "PTCRA", "LCK", "NOTCH1"),
  GMP_granulocyte = c("MPO", "ELANE", "CEBPE"),
  MYELOID_TF = c("CEBPA", "SPI1"),
  MYELOID_INFLAM = c("CYBB", "LYZ")
)
pathway_labels <- c(
  HSC_primitive = "HSC (primitive)",
  HSC_prog = "HSC (prog)",
  MPP_core = "MPP (core)",
  LMPP = "LMPP (priming)",
  CLP_early = "CLP (early)",
  B_commitment = "B commitment",
  preB_like = "preB-like",
  T_early = "T early",
  GMP_granulocyte = "GMP (granulocyte)",
  MYELOID_TF = "Myeloid TF",
  MYELOID_INFLAM = "Myeloid/inflam"
)
cond_short <- c("NPC" = "NPC", "LPS (E.coli)" = "E.coli", "LPS (Bacteroides)" = "Bact.")
hm_low <- "#4D4D4D"
hm_mid <- "white"
hm_high <- "#E69F00"

theme_science <- function(base_size = 8, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

is_cluster2 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x == 2)
  x <- trimws(as.character(x))
  x %in% c("2", "C2", "Stress", "C2 Stress")
}

read_umi_path <- function() {
  env <- Sys.getenv("XPRESS_UMI_COUNTS", "")
  if (nzchar(env)) return(path.expand(env))
  if (!file.exists(path_data_sources)) {
    stop("Missing ", path_data_sources, " (first line must be the UMI matrix path).")
  }
  lines <- readLines(path_data_sources, warn = FALSE)
  lines <- trimws(lines[nzchar(trimws(lines))])
  path_line <- lines[!grepl("^The UMI matrix", lines)][1]
  path_line <- gsub("^`|`$", "", path_line)
  if (!nzchar(path_line)) stop("No UMI path in ", path_data_sources)
  path.expand(path_line)
}

read_xpress_umi_matrix <- function(file_path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table: install.packages(\"data.table\")")
  }
  bc_row <- data.table::fread(file_path, nrows = 1L, header = FALSE, data.table = FALSE)
  barcodes <- as.character(bc_row[1L, ])
  mat <- data.table::fread(file_path, skip = 1L, header = FALSE, data.table = FALSE)
  genes <- as.character(mat[, 1L])
  count_matrix <- as.matrix(mat[, -1L, drop = FALSE])
  if (length(barcodes) != ncol(count_matrix)) {
    stop("Barcode count (", length(barcodes), ") != matrix columns (", ncol(count_matrix), ").")
  }
  rownames(count_matrix) <- genes
  colnames(count_matrix) <- barcodes
  storage.mode(count_matrix) <- "numeric"
  count_matrix
}

ensembl_count_matrix_to_seurat <- function(count_matrix,
                                           project = "xpress_scRNA",
                                           gene_cache_rds = NULL) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Install biomaRt for REBUILD_FROM_RAW=TRUE: BiocManager::install(\"biomaRt\")")
  }
  keep_idx <- grepl("^ENSG", rownames(count_matrix))
  count_matrix <- count_matrix[keep_idx, , drop = FALSE]
  gene_names <- rownames(count_matrix)
  gene_names_clean <- gsub("\\..*$", "", gene_names)

  use_cache <- !is.null(gene_cache_rds) && file.exists(gene_cache_rds)
  if (use_cache) {
    message("Loading gene mapping from cache: ", gene_cache_rds)
    cache <- readRDS(gene_cache_rds)
    gene_info <- cache$gene_info
    symbol_biotype <- cache$symbol_biotype
  } else {
    message("Querying Ensembl via biomaRt (several minutes on first run)...")
    ensembl <- biomaRt::useEnsembl(
      biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast"
    )
    gene_info <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
      filters = "ensembl_gene_id", values = gene_names_clean, mart = ensembl
    )
    symbol_biotype <- NULL
  }

  idx <- match(gene_names_clean, gene_info$ensembl_gene_id)
  biotype <- gene_info$gene_biotype[idx]
  biotype[is.na(biotype)] <- ""
  keep <- !grepl("pseudogene", biotype)
  gene_names_clean <- gene_names_clean[keep]
  count_matrix <- count_matrix[keep, , drop = FALSE]
  gene_mapping <- gene_info$hgnc_symbol[idx[keep]]
  gene_mapping[is.na(gene_mapping) | gene_mapping == ""] <-
    gene_names_clean[is.na(gene_mapping) | gene_mapping == ""]
  rownames(count_matrix) <- gene_mapping
  count_matrix <- rowsum(count_matrix, group = rownames(count_matrix))

  if (!use_cache) {
    ensembl <- biomaRt::useEnsembl(
      biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast"
    )
    symbol_biotype <- biomaRt::getBM(
      attributes = c("hgnc_symbol", "gene_biotype"),
      filters = "hgnc_symbol", values = rownames(count_matrix), mart = ensembl
    )
    if (!is.null(gene_cache_rds)) {
      saveRDS(list(gene_info = gene_info, symbol_biotype = symbol_biotype), gene_cache_rds)
      message("Gene mapping cached to: ", gene_cache_rds)
    }
  }

  pseudo_symbols <- unique(
    symbol_biotype$hgnc_symbol[
      grepl("pseudogene", symbol_biotype$gene_biotype, ignore.case = TRUE)
    ]
  )
  count_matrix <- count_matrix[!rownames(count_matrix) %in% pseudo_symbols, , drop = FALSE]
  obj <- CreateSeuratObject(counts = count_matrix, project = project)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

validate_inputs <- function(path_umi) {
  missing <- c(path_data_sources, path_barcode, path_plate)[
    !file.exists(c(path_data_sources, path_barcode, path_plate))
  ]
  if (length(missing)) {
    stop("Missing required files under input/scRNAseq/: ", paste(basename(missing), collapse = ", "))
  }
  if (REBUILD_FROM_RAW && !file.exists(path_umi)) {
    stop("UMI matrix not found: ", path_umi,
         "\nSet path in input/scRNAseq/DATA_SOURCES.md or export XPRESS_UMI_COUNTS=/path/to/matrix")
  }
}

build_seurat_from_raw <- function(path_umi) {
  meta_long <- read_csv(path_plate, show_col_types = FALSE)
  bc_ann <- read.delim(path_barcode, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
  names(bc_ann)[1:5] <- c("BCset", "plateBC", "XC_DNBPE", "XC_DNBPE_TS", "WellID")
  bc_ann$plateBC <- as.character(as.integer(bc_ann$plateBC))
  meta_long$plateBC <- as.character(meta_long$plateBC)
  bc_ann <- bc_ann %>% dplyr::left_join(meta_long, by = "plateBC")

  message("Reading UMI matrix (large file)...")
  umi <- read_xpress_umi_matrix(path_umi)
  bc_ann <- bc_ann %>%
    dplyr::filter(XC_DNBPE_TS %in% colnames(umi), !is.na(Priming_condition))
  umi <- umi[, bc_ann$XC_DNBPE_TS, drop = FALSE]

  cell_meta <- bc_ann %>%
    dplyr::mutate(cell_barcode = XC_DNBPE_TS) %>%
    as.data.frame()
  rownames(cell_meta) <- cell_meta$cell_barcode
  cell_meta <- cell_meta[colnames(umi), , drop = FALSE]

  seu <- ensembl_count_matrix_to_seurat(umi, gene_cache_rds = path_gene_cache)
  meta_cols <- intersect(
    c("ELN_ID", "Priming_condition", "Cell_sorting_date", "BCset", "plateBC",
      "WellID", "plate_slot", "XC_DNBPE", "XC_DNBPE_TS"),
    names(cell_meta)
  )
  seu <- AddMetaData(seu, metadata = cell_meta[, meta_cols, drop = FALSE])
  saveRDS(seu, path_built)
  message("Built Seurat object: ", path_built, " (", ncol(seu), " cells)")
  seu
}

process_seurat <- function(seu) {
  obj <- seu
  n_before <- ncol(obj)
  obj <- subset(obj, subset = nFeature_RNA > 1000 & percent.mt < 10)
  message(sprintf("QC filter: %d -> %d cells", n_before, ncol(obj)))

  obj <- SCTransform(obj, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
  hvgs <- VariableFeatures(obj)
  hvg_info <- if (file.exists(path_hvg_cache)) {
    tryCatch(readRDS(path_hvg_cache), error = function(e) NULL)
  } else {
    NULL
  }
  if (!is.null(hvg_info) && "gene_biotype" %in% names(hvg_info)) {
    hvgs <- intersect(hvgs, unique(hvg_info$hgnc_symbol[hvg_info$gene_biotype == "protein_coding"]))
  } else {
    message("HVG biotype cache unavailable; using all SCT variable features.")
  }
  VariableFeatures(obj) <- hvgs

  obj <- RunPCA(obj, features = hvgs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, seed.use = 42, verbose = FALSE)
  saveRDS(obj, path_processed)
  message("Processed Seurat object: ", path_processed)
  obj
}

run_pseudotime_no_c2 <- function(obj) {
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("Install slingshot: BiocManager::install(\"slingshot\")")
  }
  local({
    library(slingshot)
    library(SingleCellExperiment)

    obj2 <- subset(obj, subset = seurat_clusters != "2")
    message("Slingshot no-C2 on ", ncol(obj2), " cells")
    pca_emb <- Embeddings(obj2, "pca")[, 1:30, drop = FALSE]
    umap_emb <- Embeddings(obj2, "umap")
    sce <- SingleCellExperiment(
      assays = list(counts = GetAssayData(obj2, assay = "RNA", layer = "counts")),
      reducedDims = list(PCA = pca_emb, UMAP = umap_emb)
    )
    colData(sce)$cluster <- as.character(obj2$seurat_clusters)
    sds <- slingshot(sce, clusterLabels = colData(sce)$cluster,
                     reducedDim = "PCA", start.clus = "1")
    pt <- as.data.frame(slingPseudotime(sds))
    colnames(pt) <- paste0("Lineage_", seq_len(ncol(pt)))

    pt %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::mutate(
        Priming_condition = obj2$Priming_condition[match(cell, colnames(obj2))],
        plateBC = obj2$plateBC[match(cell, colnames(obj2))],
        seurat_clusters = as.character(obj2$seurat_clusters[match(cell, colnames(obj2))])
      )
  })
}

load_processed_object <- function(path_umi) {
  if (REBUILD_FROM_RAW) {
    if (file.exists(path_processed)) {
      message("Loading processed object: ", path_processed)
      return(readRDS(path_processed))
    }
    seu <- if (file.exists(path_built)) readRDS(path_built) else build_seurat_from_raw(path_umi)
    return(process_seurat(seu))
  }
  if (!file.exists(path_processed)) {
    pdf_ok <- all(file.exists(file.path(
      output_dir,
      c(
        "Fig7C_umap_cluster_annotation.pdf",
        "Fig7C_per_cluster_fraction_stress.pdf",
        "Fig7D_overall_pseudotime_ks_no_C2.pdf",
        "Fig7E_heatmap_lineage_panel_pseudobulk_merged_noC2.pdf"
      )
    )))
    msg <- paste0(
      "No local cache at ", path_processed, " (gitignored; not on GitHub).\n",
      if (pdf_ok) {
        "Manuscript figures are already under output/scRNAseq/*.pdf — no re-run needed.\n"
      } else {
        ""
      },
      "To build cache on this machine: REBUILD_FROM_RAW=TRUE Rscript scripts/scRNAseq.r\n",
      "(set UMI path in input/scRNAseq/DATA_SOURCES.md)."
    )
    stop(msg)
  }
  message("Loading processed object: ", path_processed)
  readRDS(path_processed)
}

plot_figures <- function(obj, pt_raw) {
  obj <- subset(obj, subset = Priming_condition %in% CONDITIONS)
  DefaultAssay(obj) <- "RNA"

  obj@meta.data$cluster_label <- factor(
    cluster_labels[as.character(obj$seurat_clusters)],
    levels = unname(cluster_labels[as.character(0:6)])
  )

  emb <- Embeddings(obj, "umap")
  umap_df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                        cluster_label = obj$cluster_label)
  umap_df <- umap_df[sample(nrow(umap_df)), , drop = FALSE]
  centroids <- umap_df %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
    dplyr::mutate(lab = sub("^C([0-6]).*", "C\\1", as.character(cluster_label)))

  p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = cluster_label)) +
    geom_point(size = 0.35, alpha = 0.9, stroke = 0) +
    geom_text(data = centroids, aes(UMAP_1, UMAP_2, label = lab),
              colour = "black", size = 2.2, inherit.aes = FALSE) +
    scale_colour_manual(values = cluster_cols, name = NULL) +
    coord_fixed() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme_science(8) + theme(legend.position = "right")
  ggsave(file.path(output_dir, "Fig7C_umap_cluster_annotation.pdf"), p_umap,
         width = 89 / 25.4, height = 70 / 25.4, device = cairo_pdf)

  meta_comp <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(
      Priming_condition = factor(.data$Priming_condition, levels = CONDITIONS),
      cluster_label = factor(cluster_labels_pub[as.character(.data$seurat_clusters)],
        levels = c("qHSC", "MPP-like", "B-primed LMPP", "GMP", "Erythroid prog",
                   "Baso/Mast prog", "Stress"))
    )
  comp_df <- meta_comp %>%
    dplyr::count(plateBC, Priming_condition, cluster_label, name = "n") %>%
    dplyr::group_by(plateBC, Priming_condition) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  comp_stress <- comp_df %>% dplyr::filter(cluster_label == "Stress")
  comp_summary_stress <- comp_stress %>%
    dplyr::group_by(Priming_condition) %>%
    dplyr::summarise(mean_prop = mean(prop), sd_prop = sd(prop), .groups = "drop")
  library(scales)
  p_stress <- ggplot(comp_stress, aes(x = Priming_condition, y = prop,
                                      colour = Priming_condition, fill = Priming_condition)) +
    geom_crossbar(data = comp_summary_stress,
      aes(x = Priming_condition, y = mean_prop, ymin = mean_prop, ymax = mean_prop),
      width = 0.55, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    geom_errorbar(data = comp_summary_stress,
      aes(x = Priming_condition, ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop),
      width = 0.18, colour = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
    geom_point(size = 2.4, shape = 21, colour = "black", stroke = 0.35,
               position = position_jitter(width = 0.1, height = 0, seed = 2)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = COND_COL, guide = "none") +
    labs(title = "Per-cluster fraction: C2 Stress", x = NULL,
         y = "Fraction within plate") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 18, hjust = 1))
  ggsave(file.path(output_dir, "Fig7C_per_cluster_fraction_stress.pdf"),
         p_stress, width = 5.5, height = 4.2)
  write_csv(comp_stress, file.path(tables_dir, "cluster_composition_stress_per_plate.csv"))

  lin_cols <- grep("^Lineage_", names(pt_raw), value = TRUE)
  pt_df <- pt_raw %>%
    dplyr::mutate(
      overall_pseudotime = rowMeans(dplyr::across(dplyr::all_of(lin_cols)), na.rm = TRUE),
      overall_pseudotime = ifelse(is.nan(overall_pseudotime), NA_real_, overall_pseudotime),
      condition_plot = factor(.data$Priming_condition, levels = rev(CONDITIONS))
    ) %>%
    dplyr::filter(
      !is.na(.data$overall_pseudotime),
      .data$Priming_condition %in% CONDITIONS
    )

  p_pt <- ggplot(pt_df, aes(overall_pseudotime, condition_plot,
                            fill = condition_plot, colour = condition_plot)) +
    geom_density_ridges(alpha = 0.45, scale = 1.15, linewidth = 0.6) +
    scale_fill_manual(values = COND_COL) +
    scale_colour_manual(values = COND_COL) +
    labs(title = "Overall pseudotime distribution by condition (no C2)",
         x = "Overall pseudotime (mean across available Slingshot lineages)", y = NULL) +
    theme_ridges(center_axis_labels = TRUE) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  ggsave(file.path(output_dir, "Fig7D_overall_pseudotime_ks_no_C2.pdf"),
         p_pt, width = 8, height = 4.2)

  build_gene_annot <- function() {
    do.call(rbind, lapply(names(gene_sets), function(k) {
      data.frame(gene = unique(gene_sets[[k]]), module_key = k, stringsAsFactors = FALSE)
    }))
  }
  obj_hm <- obj[, !is_cluster2(obj$seurat_clusters)]
  gene_annot <- build_gene_annot()
  gene_annot <- gene_annot[gene_annot$gene %in% rownames(obj_hm), , drop = FALSE]
  raw_counts <- GetAssayData(obj_hm, assay = "RNA", layer = "counts")
  plate_cond <- obj_hm@meta.data %>%
    dplyr::select(plateBC, Priming_condition) %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$Priming_condition %in% CONDITIONS) %>%
    dplyr::rename(condition = Priming_condition) %>%
    dplyr::mutate(condition = factor(.data$condition, levels = CONDITIONS)) %>%
    dplyr::arrange(condition, plateBC)
  pb_mat <- sapply(plate_cond$plateBC, function(pl) {
    cells <- colnames(obj_hm)[obj_hm$plateBC == pl]
    Matrix::rowSums(raw_counts[, cells, drop = FALSE])
  })
  hm_genes <- gene_annot$gene
  hm_mat <- pb_mat[hm_genes, , drop = FALSE]
  hm_cpm <- sweep(hm_mat, 2, colSums(pb_mat) / 1e6, "/")
  hm_z <- t(scale(t(log1p(hm_cpm))))
  hm_z[is.nan(hm_z)] <- 0
  cond_mean_z <- sapply(CONDITIONS, function(cond) {
    pls <- plate_cond$plateBC[plate_cond$condition == cond]
    rowMeans(hm_z[hm_genes, pls, drop = FALSE])
  })
  hm_m_df <- as.data.frame(cond_mean_z) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "condition", values_to = "z") %>%
    dplyr::left_join(gene_annot, by = "gene") %>%
    dplyr::mutate(
      module = factor(pathway_labels[module_key], levels = unname(pathway_labels)),
      gene = factor(.data$gene, levels = rev(gene_annot$gene)),
      condition = factor(.data$condition, levels = CONDITIONS)
    )
  write_csv(hm_m_df, file.path(tables_dir, "lineage_panel_pseudobulk_z_merged_noC2.csv"))
  p_hm <- ggplot(hm_m_df, aes(condition, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high,
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
      name = "Mean Z-score\n(across plates)") +
    scale_x_discrete(labels = cond_short) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(size = 8, angle = 0, hjust = 0),
          axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 7, face = "italic"))
  ggsave(file.path(output_dir, "Fig7E_heatmap_lineage_panel_pseudobulk_merged_noC2.pdf"),
         p_hm, width = 5.6, height = nrow(gene_annot) * 0.23 + 3.1)
}

# ── Main ──────────────────────────────────────────────────────────────────────
cat("=== scRNAseq ===\n")
cat("REBUILD_FROM_RAW:", REBUILD_FROM_RAW, "\n")
path_umi <- read_umi_path()
validate_inputs(path_umi)

obj <- load_processed_object(path_umi)
pt_raw <- run_pseudotime_no_c2(obj)
plot_figures(obj, pt_raw)
cat("Done. output -> ", output_dir, "\n")
