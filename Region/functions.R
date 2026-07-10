auto_read_effect <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if ("...1" %in% names(df)) names(df)[names(df) == "...1"] <- "gene_id"
  if ("Unnamed: 0" %in% names(df)) names(df)[names(df) == "Unnamed: 0"] <- "gene_id"
  if (!"gene_id" %in% names(df)) df <- tibble::rownames_to_column(as.data.frame(df), "gene_id")
  tibble::as_tibble(df)
}

classify_deg <- function(df, fdr_thresh = FDR_THRESH, lfc_thresh = LFC_THRESH) {
  df %>%
    mutate(
      sig = !is.na(FDR) & !is.na(logFC) & FDR <= fdr_thresh & abs(logFC) >= lfc_thresh,
      direction = case_when(
        sig & logFC > 0 ~ "Up",
        sig & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
}

get_symbol_set <- function(df, dir_label) {
  df %>% filter(direction == dir_label, !is.na(symbol), symbol != "") %>% pull(symbol) %>% unique()
}

get_gene_set <- function(df, dir_label) {
  df %>% filter(direction == dir_label) %>% pull(gene_id) %>% unique()
}

parse_ratio <- function(x) {
  parts <- str_split(x, "/", simplify = TRUE)
  as.numeric(parts[, 1]) / as.numeric(parts[, 2])
}

prep_dotplot_df <- function(df, top_n = 10) {
  if (!nrow(df)) return(tibble())
  df %>%
    mutate(
      gene_ratio = parse_ratio(GeneRatio),
      term = str_trunc(Description, 55),
      score = -log10(pmax(p.adjust, 1e-300))
    ) %>%
    slice_head(n = top_n)
}

make_dotplot <- function(df, title, top_n = 10) {
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate(
          "text", x = 0.5, y = 0.5,
          label = "No significant results",
          colour = "grey50", size = 5
        ) +
        xlim(0, 1) + ylim(0, 1) +
        theme_void() +
        ggtitle(title) +
        theme(
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
        )
    )
  }
  
  df <- df %>%
    head(top_n) %>%
    mutate(
      ratio = parse_gene_ratio(GeneRatio),
      neglog10 = -log10(p.adjust),
      Description = stringr::str_wrap(Description, width = 38),
      Description = factor(Description, levels = rev(unique(Description)))
    )
  
  ggplot(df, aes(x = ratio, y = Description)) +
    geom_point(
      aes(size = Count, fill = neglog10),
      shape = 21, colour = "grey50", stroke = 0.3
    ) +
    scale_size(range = c(3, 10), guide = "none") +
    scale_fill_gradientn(
      colours = rev(brewer.pal(11, "RdYlBu")),
      name = expression(-log[10](p.adj))
    ) +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10, lineheight = 0.95),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.4, colour = "grey85"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
}
save_venn_plot <- function(sets, file, title, fills) {
  grob <- VennDiagram::venn.diagram(
    x = sets,
    filename = NULL,
    fill = fills,
    alpha = 0.55,
    cex = 1.1,
    cat.cex = 1.15,
    cat.fontface = "bold",
    fontface = "bold",
    main = title,
    main.cex = 1.2,
    margin = 0.08
  )
  png(file, width = 1800, height = 900, res = 180)
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
}

ens_to_entrez <- function(ens_ids) {
  mapped <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = ens_ids,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  unique(unname(mapped[!is.na(mapped)]))
}


run_ora <- function(gene_ens, universe_ens, label = "") {
  gene_entrez <- ens_to_entrez(gene_ens)
  bg_entrez   <- ens_to_entrez(universe_ens)
  
  if (length(gene_entrez) < 5) return(list(KEGG = data.frame(), GOBP = data.frame()))
  
  kegg <- tryCatch(
    enrichKEGG(
      gene = gene_entrez,
      universe = bg_entrez,
      organism = "mmu",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    ),
    error = function(e) NULL
  )
  
  gobp <- tryCatch(
    enrichGO(
      gene = gene_entrez,
      universe = bg_entrez,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  
  list(
    KEGG = if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0)
      as.data.frame(kegg) %>% filter(p.adjust <= 0.05) %>% arrange(p.adjust) else data.frame(),
    GOBP = if (!is.null(gobp) && nrow(as.data.frame(gobp)) > 0)
      as.data.frame(gobp) %>% filter(p.adjust <= 0.05) %>% arrange(p.adjust) else data.frame()
  )
}

prep_deg <- function(df) {
  if (!"gene_id" %in% names(df)) df$gene_id <- rownames(df)
  df %>%
    mutate(
      sig = FDR <= FDR_THRESH & abs(logFC) >= LFC_THRESH,
      direction = case_when(
        sig & logFC > 0 ~ "Up",
        sig & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
}


add_biomart_mouse_anno <- function(res,
                                   id_col = NULL,
                                   mirror = c("www", "useast", "uswest", "asia"),
                                   keep_first = TRUE) {
  
  ## usage:
  # res_anno <- add_biomart_mouse_anno(res)
  # head(res_anno)
  ids <- if (is.null(id_col)) rownames(res) else res[[id_col]]
  if (is.null(ids)) stop("No IDs found: provide rownames(res) or set id_col.")
  ids_clean <- sub("\\..*$", "", as.character(ids))
  
  mart <- tryCatch(
    biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",
                        mirror = if (mirror == "www") NULL else mirror),
    error = function(e) biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",
                                            mirror = "useast")
  )
  
  ann <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    filters = "ensembl_gene_id",
    values  = unique(ids_clean),
    mart    = mart
  )
  if (keep_first) ann <- ann[!duplicated(ann$ensembl_gene_id), ]
  
  sym  <- setNames(ann$mgi_symbol,  ann$ensembl_gene_id)
  desc <- setNames(ann$description, ann$ensembl_gene_id)
  
  out <- res
  out$ensembl_gene_id <- ids_clean
  out$symbol          <- unname(sym[ids_clean])
  out$description     <- unname(desc[ids_clean])
  
  out <- out[, c("ensembl_gene_id", "symbol", "description",
                 setdiff(colnames(out), c("ensembl_gene_id","symbol","description")))]
  out
}

create_volcano <- function(res_obj, FC = "logFC", title) {
  stopifnot(is.data.frame(res_obj))
  if (!FC %in% names(res_obj)) stop(sprintf("FC column '%s' not found in res_obj.", FC))
  if (!"FDR" %in% names(res_obj)) stop("Column 'FDR' not found in res_obj.")
  if (!"symbol" %in% names(res_obj)) stop("Column 'symbol' not found in res_obj.")
  
  volcano_data <- data.frame(
    log2FC    = res_obj[[FC]],
    neglog10p = -log10(res_obj$FDR),
    gene      = res_obj$symbol
  )
  
  volcano_data <- volcano_data[complete.cases(volcano_data[, c("log2FC","neglog10p")]), ]
  
  fc_thr  <- 1
  p_thr   <- 0.05
  sig_thr <- -log10(p_thr)
  
  volcano_data$category <- "not significant"
  volcano_data$category[volcano_data$log2FC >  fc_thr & volcano_data$neglog10p > sig_thr] <- "upregulated"
  volcano_data$category[volcano_data$log2FC < -fc_thr & volcano_data$neglog10p > sig_thr] <- "downregulated"
  
  top_labs <- volcano_data %>%
    dplyr::filter(category != "not significant") %>%
    dplyr::arrange(desc(neglog10p)) %>%
    head(10)
  
  ggplot(volcano_data, aes(x = log2FC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_text_repel(data = top_labs, aes(label = gene),
                    size = 3, show.legend = FALSE) +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed") +
    geom_hline(yintercept = sig_thr, linetype = "dashed") +
    labs(
      title = title,
      x = expression(Log[2]~"fold change"),
      y = expression(-Log[10]~P)
    ) +
    scale_color_manual(values = c("downregulated"="blue",
                                  "not significant"="grey50",
                                  "upregulated"="red")) +
    theme_bw() +
    theme(legend.title = element_blank())
}