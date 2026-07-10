sym2entrez <- function(symbols) {
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  ids <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  unique(as.character(ids[!is.na(ids)]))
}

extract_degs <- function(res, fc_col, fdr_cut = 0.05, lfc_cut = 0.5) {
  stopifnot(is.data.frame(res))
  stopifnot(fc_col %in% names(res))
  stopifnot("FDR" %in% names(res))
  stopifnot("symbol" %in% names(res))
  
  fc  <- suppressWarnings(as.numeric(res[[fc_col]]))
  fdr <- suppressWarnings(as.numeric(res$FDR))
  sym <- as.character(res$symbol)
  
  keep <- !is.na(fdr) & !is.na(fc) & !is.na(sym) & nzchar(sym)
  sel  <- keep & fdr < fdr_cut & abs(fc) > lfc_cut
  
  list(
    up = unique(sym[sel & fc > 0]),
    down = unique(sym[sel & fc < 0]),
    n_up = sum(sel & fc > 0, na.rm = TRUE),
    n_down = sum(sel & fc < 0, na.rm = TRUE)
  )
}

run_ora <- function(symbols, universe) {
  entrez <- sym2entrez(symbols)
  if (length(entrez) < 5) return(NULL)
  enrichGO(
    gene = entrez,
    universe = universe,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE,
    minGSSize = 10,
    maxGSSize = 500
  )
}

shorten_term <- function(x, maxchar = 55) {
  ifelse(nchar(x) > maxchar, paste0(substr(x, 1, maxchar - 3), "..."), x)
}

get_top_terms <- function(ora_res, n = 8) {
  if (is.null(ora_res)) return(NULL)
  df <- as.data.frame(ora_res)
  if (!nrow(df)) return(NULL)
  head(df[order(df$p.adjust), , drop = FALSE], n)
}
get_padj_for_terms <- function(ora_res, terms) {
  if (is.null(ora_res)) return(setNames(rep(NA_real_, length(terms)), terms))
  df <- as.data.frame(ora_res)
  out <- setNames(rep(NA_real_, length(terms)), terms)
  hit <- match(terms, df$Description)
  out[!is.na(hit)] <- df$p.adjust[hit[!is.na(hit)]]
  out
}


make_dotplot <- function(
    df, title_text,
    strain_colors = c(CAST = "#E69F00", WSB = "#56B4E9", NZO = "#009E73"),
    label_width = 58
) {
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(grid)
  
  tag <- c("Higher than B6" = "H", "Lower than B6" = "L")
  
  df <- df %>%
    mutate(
      direction = factor(direction, levels = c("Higher than B6", "Lower than B6")),
      strain    = factor(strain, levels = c("CAST", "WSB", "NZO")),
      y_lab     = paste0("[", tag[as.character(direction)], "] ",
                         str_trunc(Description_short, width = label_width, side = "right"))
    )
  
  ord <- df %>%
    group_by(direction, y_lab) %>%
    summarise(best_padj = min(p.adjust), .groups = "drop") %>%
    arrange(direction, best_padj)
  
  hi <- ord$y_lab[ord$direction == "Higher than B6"]
  lo <- ord$y_lab[ord$direction == "Lower than B6"]
  
  df$y_lab <- factor(df$y_lab, levels = c(rev(lo), rev(hi)))
  
  ggplot(df, aes(x = strain, y = y_lab)) +
    geom_point(aes(size = GR_num, color = strain), alpha = 0.95, shape = 16) +
    facet_grid(
      rows = vars(direction),
      scales = "free_y",
      space = "free_y"
    ) +
    scale_color_manual(values = strain_colors, guide = "none") +
    scale_size_continuous(
      name   = NULL,
      range  = c(3, 10),
      breaks = c(0.02, 0.04, 0.06),
      labels = scales::percent_format(accuracy = 0.1),
      guide  = guide_legend(
        override.aes = list(color = "black", alpha = 1)
      )
    ) +
    scale_x_discrete(expand = expansion(add = c(0.6, 0.8))) +
    scale_y_discrete(expand = expansion(add = c(0.35, 0.35))) +
    labs(title = title_text, x = "Strain (vs B6)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "grey94", color = NA),
      panel.background = element_rect(fill = "grey94", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      
      axis.line.x = element_line(linewidth = 0.9, color = "black"),
      axis.line.y = element_line(linewidth = 0.9, color = "black"),
      axis.ticks.y = element_blank(),
      
      axis.text.y = element_text(size = 9.5, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 15, face = "bold"),
      
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      
      strip.background = element_blank(),
      strip.text.y = element_text(size = 13, face = "bold", angle = -90),
      strip.placement = "outside",
      
      panel.spacing.y = unit(1.2, "lines"),
      
      legend.position = "right",
      legend.key = element_rect(fill = "grey94", color = NA),
      legend.text = element_text(size = 11)
    )
}
make_dotplot <- function(df, title_text,
                         strain_colors = c(CAST = "#E69F00", WSB = "#56B4E9", NZO = "#009E73"),
                         label_width = 58) {
  tag <- c("Higher than B6" = "H", "Lower than B6" = "L")
  
  df <- df %>%
    mutate(
      direction = factor(direction, levels = c("Higher than B6", "Lower than B6")),
      strain    = factor(strain, levels = c("CAST", "WSB", "NZO")),
      y_lab     = paste0("[", tag[as.character(direction)], "] ",
                         str_trunc(Description_short, width = label_width, side = "right"))
    )
  
  ord <- df %>%
    group_by(direction, y_lab) %>%
    summarise(best_padj = min(p.adjust), .groups = "drop") %>%
    arrange(direction, best_padj)
  
  hi <- ord$y_lab[ord$direction == "Higher than B6"]
  lo <- ord$y_lab[ord$direction == "Lower than B6"]
  df$y_lab <- factor(df$y_lab, levels = c(rev(lo), rev(hi)))
  
  ggplot(df, aes(x = strain, y = y_lab)) +
    geom_point(aes(size = GR_num, color = strain), alpha = 0.95, shape = 16) +
    facet_grid(rows = vars(direction), scales = "free_y", space = "free_y") +
    scale_color_manual(values = strain_colors, guide = "none") +
    scale_size_continuous(
      name   = "Gene Ratio",
      range  = c(3, 10),
      breaks = c(0.02, 0.04, 0.06),
      labels = scales::percent_format(accuracy = 0.1),
      guide  = guide_legend(override.aes = list(color = "black", alpha = 1))
    ) +
    scale_x_discrete(expand = expansion(add = c(0.6, 0.8))) +
    scale_y_discrete(expand = expansion(add = c(0.35, 0.35))) +
    labs(title = title_text, x = "Strain (vs B6)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "grey94", color = NA),
      panel.background = element_rect(fill = "grey94", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.line.x      = element_line(linewidth = 0.9, color = "black"),
      axis.line.y      = element_line(linewidth = 0.9, color = "black"),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_text(size = 10, face = "bold"),
      axis.text.x      = element_text(size = 11, face = "bold"),
      axis.title.x     = element_text(size = 13, face = "bold"),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text.y     = element_text(size = 12, face = "bold", angle = 0),
      strip.placement  = "outside",
      panel.spacing.y  = unit(1.2, "lines"),
      legend.position  = "right",
      legend.key       = element_rect(fill = "grey94", color = NA),
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10),
      plot.margin      = margin(t = 15, r = 120, b = 15, l = 15)
    )
}
make_deg_bar <- function(deg_counts_df, tissue_filter, title_text,
                         strain_colors = c(CAST = "#E69F00", WSB = "#56B4E9", NZO = "#009E73"),
                         bg = "grey92") {
  
  df_tissue <- deg_counts_df %>%
    filter(tissue == tissue_filter) %>%
    mutate(strain = factor(strain, levels = c("CAST", "WSB", "NZO")))
  
  # Dynamic y-axis limits based on actual data (rounded to nearest 200)
  ymax_r <- ceiling(max(df_tissue$n_up)   / 200) * 200
  ymin_r <- ceiling(max(df_tissue$n_down) / 200) * 200
  
  deg_long <- df_tissue %>%
    mutate(n_down = -n_down) %>%
    pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "count") %>%
    mutate(direction = recode(direction, n_up = "Higher than B6", n_down = "Lower than B6"))
  
  ggplot(deg_long, aes(strain, count, fill = strain)) +
    geom_col(width = 0.62, position = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "grey20") +
    scale_fill_manual(values = strain_colors, name = "Strain") +
    scale_y_continuous(
      limits = c(-(ymin_r + 50), ymax_r + 50),
      breaks = seq(-ymin_r, ymax_r, 200),
      labels = \(x) comma(abs(x)),
      expand = c(0, 0),
      sec.axis = sec_axis(
        ~ .,
        breaks = c(-ymin_r * 0.6, ymax_r * 0.6),
        labels = c("Lower than B6 \u2193", "Higher than B6 \u2191")
      )
    ) +
    labs(
      title    = title_text,
      subtitle = "FDR \u2264 0.05, |logFC| \u2265 0.5",
      x        = "Strain",
      y        = "Number of DEGs"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.0, "lines"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 11, face = "bold", colour = "grey40", hjust = 0.5, margin = margin(b = 8)),
      strip.text = element_text(size = 13, face = "bold"),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = unit(0.12, "in"),
      plot.margin = margin(t = 15, r = 100, b = 15, l = 20)
    )
}
