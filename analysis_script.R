# according to https://stackoverflow.com/questions/68620638/how-to-execute-r-inside-snakemake
if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(
      input = "list",
      output = "list",
      params = "list",
      wildcards = "list",
      threads = "numeric",
      log = "list",
      resources = "list",
      config = "list",
      rule = "character",
      bench_iteration = "numeric",
      scriptdir = "character",
      source = "function"
    )
  )
  # si exécution manuelle modifier les paths et les faire correspondre à l'analyse actuelle
  snakemake <- Snakemake(
    input = list(
      quants = "./quants_human/ALL_SAMPLES/ALL_SAMPLES.transcript_model_grouped_counts.tsv",
      gtf = "./quants_human/ALL_SAMPLES/ALL_SAMPLES.transcript_models.gtf"
    ),
    output = list(),
    params = list(
      design = "./design_matrix.tsv",
      comparison = "./contrasts.csv",
      output_dir = "analysis",
      batch = "run"
    ),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list(),
    rule = "",
    bench_iteration = 1,
    scriptdir = "",
    source = function(...) {{ wd <- getwd()
      setwd(snakemake@scriptdir)
      source(...)
      setwd(wd) }}
  )
} else {
  snakemake@source(".Rprofile")
}


library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(aggregation)
library(hexbin)
library(crayon)

library(GenomicFeatures)
library(vsn)
library(DESeq2)
library(DEXSeq)
library(DRIMSeq)
library(biomaRt)
library(IsoformSwitchAnalyzeR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stageR)

rld_pca <- function(rld, config, ntop = 500) {
  rv <- matrixStats::rowVars(assay(rld))
  selected_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t(assay(rld)[selected_genes, ])
  pc <- prcomp(mat)
  eig <- (pc$sdev)^2
  variance <- eig * 100 / sum(eig)

  PCAdata <- as.data.frame(pc$x)
# Join with condition, on name, to be sure of matches between condition and sample
  PCAdata <- PCAdata %>%
    tibble::rownames_to_column(var = "sampleID") %>%
    inner_join(config, by = "sampleID") %>%
    tibble::column_to_rownames(var = "sampleID")
  list("data" = PCAdata, "variance" = variance)
}

# TODO find a way to integrate the conditions automatically
plot_isoforms <- function(matrix_tx, txdf, design_matrix, genes, sig_iso = c()) {
  tmp_data <- matrix_tx %>%
    left_join(
      txdf %>%
      dplyr::rename(
               "isoform_id" = "TXNAME",
               "gene_id" = "GENEID"
             ),
      by = "isoform_id"
    ) %>%
    dplyr::filter(gene_id %in% genes) %>%
    mutate(
      transcript_level = ifelse(
        isoform_id %in% sig_iso,
        "significant",
        "not_significant"
      )
    )
  diff_conditions <- unique(design_matrix$condition)

  for(cond in diff_conditions) {
    samples_to_mean <- design_matrix %>%
                    dplyr::filter(condition == cond) %>%
                    pull(sampleID)
    tmp_data[[cond]] <- rowMeans(tmp_data[samples_to_mean])
  }
  plot_data <- tmp_data %>%
    dplyr::select(all_of(diff_conditions), isoform_id, gene_id, transcript_level) %>%
    pivot_longer(!c(isoform_id, gene_id, transcript_level), names_to = "condition", values_to = "mean")


  plot_tx <- ggplot(plot_data, aes(x = isoform_id, y = mean, fill = condition, color = transcript_level)) +
    geom_bar(stat = "identity", position = "dodge", linewidth = 3) +
    scale_color_manual(values = c("significant" = "black", "not_significant" = "white"))
  if(length(genes) == 1) {
    plot_tx <- plot_tx +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      facet_wrap(~gene_id, scales = "free_x")
  } else if (length(genes) < 8) {
    plot_tx <- plot_tx +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme(axis.text.x = element_text(size = 5)) +
      facet_wrap(~gene_id, scales = "free")
  } else {
    plot_tx <- plot_tx +
      theme(axis.text.x = element_blank()) +
      facet_wrap(~gene_id, scales = "free")
  }

  plot_tx
}


design_matrix <- read_tsv(snakemake@params[["design"]], col_types = list(.default = col_character()))
matrix_tx <- read_tsv(snakemake@input[["quants"]])
comparison_df <- read.table(snakemake@params[["comparison"]], sep = ",")

if (!(all(design_matrix$sampleID == colnames(matrix_tx[, -1]))))
  stop("The sample names defined in the design matrix do not match the sample names in the quants file")
if(!(unique(comparison_df[, 1]) %in% colnames(design_matrix)))
  stop("The specified variable is not present in the design_matrix")

colnames(matrix_tx)[1] <- "isoform_id"


# PCA ------------------------------
exploratory_dds <- DESeqDataSetFromMatrix(
  round(matrix_tx %>% column_to_rownames("isoform_id") %>% as.matrix()),
  colData = design_matrix %>% column_to_rownames("sampleID"),
  ~condition
)
exploratory_dds <- estimateSizeFactors(exploratory_dds)
rld <- rlog(exploratory_dds)
pca_data <- rld_pca(rld, design_matrix)#, nrow(log_tr))
pca_plot <- ggplot(
  pca_data$data,
  aes(
    x = PC1,
    y = PC2,
    label = rownames(pca_data$data)
  )
) +
  ylab(paste0("PC2: ", round(pca_data$variance[2], 1), "% variance")) +
  xlab(paste0("PC1: ", round(pca_data$variance[1], 1), "% variance")) +
  coord_fixed() +
  geom_point(aes(color = condition), size = 5) +
  geom_label_repel()

ggsave(
  file.path(
    snakemake@params[["output_dir"]],
    "pca_plot.png"
  ),
  pca_plot
)




# Filter with DRIMseq -------------------------------
gtf <- snakemake@input[["gtf"]]
txdb_filename <- "transcript_models.sqlite"
txdb <- makeTxDbFromGFF(gtf)
invisible(saveDb(txdb, txdb_filename))
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")

for (i in nrow(comparison_df)) {
  cat(bold(paste(
    "Studying",
    comparison_df[i, 1],
    ":",
    comparison_df[i, 2],
    "vs",
    comparison_df[i, 3],
    "\n"
  )))
  batch_list <- stringr::str_split_1(snakemake@params[["batch"]], ",")
  if (comparison_df[i, 1] %in% batch_list) {
    stop("One batch variable cannot also be the studied variable")
  }

  cur_samples <- design_matrix %>%
    dplyr::filter(.data[[comparison_df[i, 1]]] %in% comparison_df[i, 2:3]) %>%
    pull(sampleID)
  matrix_tx_cur <- matrix_tx[, c("isoform_id", cur_samples)]
  counts <- left_join(
    matrix_tx_cur %>% dplyr::rename("feature_id" = "isoform_id"),
    txdf %>% dplyr::rename("feature_id" = "TXNAME"),
    by = "feature_id"
  ) %>%
    dplyr::rename("gene_id" = "GENEID") %>%
    relocate(gene_id, .after = feature_id)

  # TODO à paramétriser depuis le fichier config.yaml ?
  drim <- dmDSdata(
    counts = counts %>% as.data.frame(),
    samples = design_matrix %>%
      dplyr::filter(sampleID %in% cur_samples) %>%
      dplyr::rename("sample_id" = "sampleID") %>%
      as.data.frame()
  )
  drim <- dmFilter(
    drim,
    min_samps_feature_expr = 3, min_feature_expr = 3,
    min_samps_feature_prop = 3, min_feature_prop = 0.1,
    min_samps_gene_expr = 6, min_gene_expr = 3
  )

  if (batch_list != c("")) {
    cur_full_model_wo_batch <- paste0("~ sample + exon + ", comparison_df[i, 1], ":exon")
    batch_f <- paste0(batch_list, ":exon", collapse = " + ")
    cur_full_model <- formula(paste0(cur_full_model_wo_batch, " + ", batch_f))

    cur_reduced_model_wo_batch <- "~ sample + exon"
    cur_reduced_model <- formula(paste0(cur_reduced_model_wo_batch, " + ", batch_f))
  } else {
    cur_full_model <- formula(paste0("~ sample + exon + ", comparison_df[i, 1], ":exon"))
    cur_reduced_model <- formula("~ sample + exon")
  }

  cat(bold(paste("Full formula :", deparse1(cur_full_model), "\n")))
  cat(bold(paste("Reduced formula :", deparse1(cur_reduced_model), "\n")))



  # Inférence DEXSeq -----------------------------------
  sample_data <- DRIMSeq::samples(drim)
  count_data <- round(as.matrix(counts(drim)[, -c(1, 2)]))
  dex <- DEXSeqDataSet(
    countData = count_data,
    sampleData = sample_data,
    design = cur_full_model,
    featureID = counts(drim)$feature_id,
    groupID = counts(drim)$gene_id
  )
  dex <- DEXSeq::estimateSizeFactors(dex)
  dex <- DEXSeq::estimateDispersions(dex)
  dex <- DEXSeq::testForDEU(dex, reducedModel = cur_reduced_model)
  dex_res <- DEXSeqResults(dex, independentFiltering = FALSE)
  qval <- perGeneQValue(dex_res)
  dex_res_g <- data.frame(gene = names(qval), qval)

  # Descente niveau transcrit avec stageR ------------------

  p_confirmation <- matrix(dex_res$pvalue, ncol = 1)
  dimnames(p_confirmation) <- list(dex_res$featureID, "transcript")
  tx2gene <- as.data.frame(dex_res[, c("featureID", "groupID")])
  p_screen <- qval

  stager_obj <- stageRTx(
    pScreen = p_screen,
    pConfirmation = p_confirmation,
    pScreenAdjusted = TRUE,
    tx2gene = tx2gene
  )
  stager_obj <- stageWiseAdjustment(stager_obj, method = "dtu", alpha = 0.05)
  dex_padj <- getAdjustedPValues(
    stager_obj,
    order = FALSE,
    onlySignificantGenes = TRUE
  )

  # Plots et sauvegarde
  if (is.null(dex_padj)) {
    cat(bold(red("No DTU genes were found\n")))
    file.create(
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "signif_genes_and_isoforms.png")
      )
    )
    file.create(
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "signif_genes_and_isoforms.svg")
      )
    )
    file.create(
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "adjusted_gene_transcript_pval_005.csv")
      )
    )
  } else {
    top_genes <- dex_padj %>%
      pull(geneID) %>%
      unique()
    sig_iso <- dex_padj %>%
      filter(transcript < 0.05) %>%
      pull(txID)
    plot_iso <- plot_isoforms(
      matrix_tx_cur,
      txdf,
      design_matrix %>%
        filter(sampleID %in% cur_samples),
      top_genes,
      sig_iso
    )
    ggsave(
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "signif_genes_and_isoforms.png")
      ),
      plot_iso,
      dpi = 300
    )
    ggsave(
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "signif_genes_and_isoforms.svg")
      ),
      plot_iso
    )

    write_csv(
      dex_padj,
      file.path(
        snakemake@params[["output_dir"]],
        paste0(comparison_df[i, 1], "_", comparison_df[i, 2], "-vs-", comparison_df[i, 3], "_", "adjusted_gene_transcript_pval_005.csv")
      )
    )
  }
}
# Aucun transcrit n'est significatif 🥲
