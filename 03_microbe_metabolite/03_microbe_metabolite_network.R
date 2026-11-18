# ================================================================
# SS & NSS 微生物-代谢物 Spearman 相关 → 网络+热图
# ================================================================

library(igraph)
library(pheatmap)
library(svglite)

# --- 自动定位到脚本所在目录（替代本机绝对路径，便于他人复现）---
local({
  args <- commandArgs(trailingOnly = FALSE)
  fileflag <- grep("^--file=", args, value = TRUE)
  if (length(fileflag) > 0) {
    setwd(dirname(normalizePath(sub("^--file=", "", fileflag[1]))))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             rstudioapi::isAvailable() &&
             nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    message("请先将工作目录设置为脚本所在目录再运行本脚本")
  }
})
# --- 以上为路径适配，以下为原有分析代码 ---
dir.create("output", showWarnings = FALSE)
MM <- 1/25.4

POS_COLOR <- "#D95F02"; NEG_COLOR <- "#4575B4"
MICROBE_COL <- "#4DBBD5"; METAB_COL <- "#7FC97F"

# ================================================================
run_analysis <- function(soil_label, ss_cols, nss_cols) {
  cat(sprintf("\n========== %s ==========\n", soil_label))

  # ---- 1. 读取微生物 ----
  df <- read.delim("../data/microbial_OTU_table.xls", sep = "\t", header = TRUE, check.names = FALSE)
  data_cols <- grep("^W[0-9]+_[0-9]+$", names(df), value = TRUE)
  df$genus_clean <- gsub("^g__", "", df$genus)
  genus_abund <- aggregate(df[, data_cols], by = list(genus = df$genus_clean), FUN = sum)
  rownames(genus_abund) <- genus_abund$genus; genus_abund$genus <- NULL

  # 取对应土壤的样本
  sample_cols <- if (soil_label == "SS") ss_cols else nss_cols
  genus_sub <- genus_abund[, unlist(sample_cols)]

  # 筛选
  is_uncl <- grepl("unclassified", rownames(genus_sub), ignore.case = TRUE)
  genus_clean <- genus_sub[!is_uncl, ]
  rel <- sweep(genus_clean, 2, colSums(genus_clean), "/") * 100
  keep1 <- rowSums(genus_clean > 0) / ncol(genus_clean) >= 0.5
  keep2 <- rowMeans(rel) > 0.01
  genus_filtered <- genus_clean[keep1 & keep2, ]
  mean_micro <- sort(rowMeans(genus_filtered), decreasing = TRUE)
  top_genus <- names(head(mean_micro, 50))
  micro_mat <- genus_filtered[top_genus, ]
  cat(sprintf("  Microbes: %d genera\n", nrow(micro_mat)))

  # ---- 2. 读取代谢物 (使用已筛选的) ----
  metab_file <- if (soil_label == "SS") {
    "../02_metabolite_physiology/metabolite_SS.csv"
  } else {
    "../02_metabolite_physiology/metabolite_NSS.csv"
  }
  metab_raw <- read.csv(metab_file, header = TRUE, check.names = FALSE)
  rownames(metab_raw) <- metab_raw$Metabolite
  metab_mat <- metab_raw[, unlist(sample_cols), drop = FALSE]
  metab_mat <- as.data.frame(lapply(metab_mat, as.numeric))
  rownames(metab_mat) <- rownames(metab_raw)
  cat(sprintf("  Metabolites: %d\n", nrow(metab_mat)))

  # ---- 3. Spearman 相关 ----
  n_micro <- nrow(micro_mat); n_metab <- nrow(metab_mat)
  r_mat <- matrix(NA, n_micro, n_metab)
  p_mat <- matrix(NA, n_micro, n_metab)
  pb <- txtProgressBar(min = 0, max = n_micro, style = 3)
  for (i in 1:n_micro) {
    for (j in 1:n_metab) {
      test <- cor.test(as.numeric(micro_mat[i,]), as.numeric(metab_mat[j,]),
                       method = "spearman", exact = FALSE)
      r_mat[i, j] <- test$estimate
      p_mat[i, j] <- test$p.value
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rownames(r_mat) <- rownames(micro_mat); colnames(r_mat) <- rownames(metab_mat)
  rownames(p_mat) <- rownames(micro_mat); colnames(p_mat) <- rownames(metab_mat)

  p_raw_mat <- p_mat
  n_sig <- sum(p_raw_mat < 0.05, na.rm = TRUE)
  cat(sprintf("  Tests: %d, Raw P<0.05: %d\n", n_micro * n_metab, n_sig))

  # ---- 4. 热图 ----
  has_sig_m <- rowSums(p_raw_mat < 0.05, na.rm = TRUE) > 0
  has_sig_b <- colSums(p_raw_mat < 0.05, na.rm = TRUE) > 0

  if (sum(has_sig_m) >= 2 && sum(has_sig_b) >= 2) {
    r_plot <- r_mat[has_sig_m, has_sig_b, drop = FALSE]
    p_plot <- p_raw_mat[has_sig_m, has_sig_b, drop = FALSE]
    sig_mark <- matrix("", nrow(r_plot), ncol(r_plot))
    sig_mark[p_plot < 0.01] <- "**"
    sig_mark[p_plot >= 0.01 & p_plot < 0.05] <- "*"

    fig_h <- max(80, nrow(r_plot) * 10) * MM
    svglite(sprintf("output/%s_microbe_metabolite_heatmap.svg", soil_label),
            width = 300 * MM, height = fig_h)
    pheatmap(r_plot,
             cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
             color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
             breaks = seq(-1, 1, length.out = 101),
             display_numbers = sig_mark, number_color = "black",
             fontsize_number = 6, fontsize = 6, fontsize_col = 7, fontsize_row = 6,
             main = sprintf("%s: Microbe-Metabolite (Raw P<0.05)", soil_label),
             angle_col = 45, border_color = NA,
             cellwidth = 14, cellheight = 8,
             legend_breaks = c(-1, -0.5, 0, 0.5, 1),
             legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
    dev.off()
    cat(sprintf("  -> %s_microbe_metabolite_heatmap.svg (%d x %d)\n",
                soil_label, nrow(r_plot), ncol(r_plot)))
  } else {
    cat(sprintf("  %s: not enough significant genera/metabolites for heatmap\n", soil_label))
  }

  # ---- 5. 建网络 ----
  edge_list <- data.frame()
  for (i in 1:nrow(r_mat)) for (j in 1:ncol(r_mat)) {
    if (!is.na(p_raw_mat[i,j]) && p_raw_mat[i,j] < 0.05 && abs(r_mat[i,j]) >= 0.7) {
      edge_list <- rbind(edge_list, data.frame(
        Node1 = rownames(r_mat)[i], Node2 = colnames(r_mat)[j],
        Coefficient = round(r_mat[i,j], 3), P_value = signif(p_raw_mat[i,j], 3),
        Direction = ifelse(r_mat[i,j] > 0, "positive", "negative"),
        stringsAsFactors = FALSE))
    }
  }

  if (nrow(edge_list) == 0) {
    cat(sprintf("  %s: no significant edges\n", soil_label))
    return(list(n_sig = 0))
  }

  cat(sprintf("  Edges: %d (pos=%d, neg=%d)\n", nrow(edge_list),
              sum(edge_list$Direction == "positive"),
              sum(edge_list$Direction == "negative")))

  all_nodes <- unique(c(edge_list$Node1, edge_list$Node2))
  g <- make_empty_graph(n = length(all_nodes), directed = FALSE)
  V(g)$name <- all_nodes
  V(g)$type <- ifelse(all_nodes %in% rownames(r_mat), "Microbe", "Metabolite")
  V(g)$color <- ifelse(V(g)$type == "Microbe", MICROBE_COL, METAB_COL)

  for (k in 1:nrow(edge_list)) {
    g <- add_edges(g, c(edge_list$Node1[k], edge_list$Node2[k]),
                   weight = abs(edge_list$Coefficient[k]),
                   coefficient = edge_list$Coefficient[k],
                   sign = edge_list$Direction[k])
  }

  # ---- 6. 双圈网络绘图 ----
  micro_names <- sort(unique(edge_list$Node1))
  metab_names <- sort(unique(edge_list$Node2))
  n_mic <- length(micro_names); n_met <- length(metab_names)

  r_outer <- 1.0   # 外圈=代谢物
  r_inner <- 0.65  # 内圈=微生物

  angles_met <- seq(0, 2 * pi, length.out = n_met + 1)[1:n_met]
  coords_met <- data.frame(name = metab_names, x = r_outer * cos(angles_met),
                           y = r_outer * sin(angles_met), stringsAsFactors = FALSE)
  angles_mic <- seq(0, 2 * pi, length.out = n_mic + 1)[1:n_mic]
  coords_mic <- data.frame(name = micro_names, x = r_inner * cos(angles_mic),
                           y = r_inner * sin(angles_mic), stringsAsFactors = FALSE)
  all_coords <- rbind(coords_met, coords_mic)
  rownames(all_coords) <- c(coords_met$name, coords_mic$name)

  # 节点度 → 大小
  degs <- degree(g)
  deg_met <- degs[coords_met$name]; deg_mic <- degs[coords_mic$name]
  deg_all <- c(deg_met, deg_mic)
  dbrk <- c(1, 2, 5)  # Small=1, Medium=2-4, Large≥5
  sz_fn <- function(d) ifelse(d >= dbrk[3], 4.0, ifelse(d >= dbrk[2], 2.4, 1.2))
  coords_met$cex <- sz_fn(deg_met); coords_mic$cex <- sz_fn(deg_mic)

  n_pos <- sum(edge_list$Direction == "positive")
  n_neg <- sum(edge_list$Direction == "negative")

  svglite(sprintf("output/%s_microbe_metabolite_network.svg", soil_label),
          width = 100 * MM, height = 100 * MM,
          pointsize = 8, system_fonts = list(serif = "Times New Roman"))

  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white", family = "serif")
  plot(0, 0, type = "n", xlim = c(-1.35, 1.35), ylim = c(-1.25, 1.25),
       xlab = "", ylab = "", axes = FALSE, asp = 1)

  # 边
  for (k in 1:nrow(edge_list)) {
    n1 <- edge_list$Node1[k]; n2 <- edge_list$Node2[k]
    x1 <- all_coords[n1, "x"]; y1 <- all_coords[n1, "y"]
    x2 <- all_coords[n2, "x"]; y2 <- all_coords[n2, "y"]
    col <- ifelse(edge_list$Direction[k] == "positive", POS_COLOR, NEG_COLOR)
    alpha_val <- pmax(0.15, abs(edge_list$Coefficient[k]) * 0.5)
    rgb_vals <- col2rgb(col)[, 1] / 255
    edge_col <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha_val)
    lines(c(x1, x2), c(y1, y2), col = edge_col, lwd = abs(edge_list$Coefficient[k]) * 1.2)
  }

  # 外圈代谢物
  points(coords_met$x, coords_met$y, pch = 21, bg = adjustcolor(METAB_COL, 0.85),
         col = "grey40", cex = coords_met$cex)
  # 内圈微生物
  points(coords_mic$x, coords_mic$y, pch = 21, bg = adjustcolor(MICROBE_COL, 0.85),
         col = "grey40", cex = coords_mic$cex)

  # 代谢物标签（外圈圆心）
  for (i in 1:n_met) {
    ang <- angles_met[i]; tx <- r_outer * cos(ang); ty <- r_outer * sin(ang)
    text(tx, ty, labels = coords_met$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }
  # 微生物标签（内圈圆心）
  for (i in 1:n_mic) {
    ang <- angles_mic[i]; tx <- r_inner * cos(ang); ty <- r_inner * sin(ang)
    text(tx, ty, labels = coords_mic$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }

  # 图例
  n_sm <- sum(c(deg_met, deg_mic) < dbrk[2])
  n_md <- sum(c(deg_met, deg_mic) >= dbrk[2] & c(deg_met, deg_mic) < dbrk[3])
  n_lg <- sum(c(deg_met, deg_mic) >= dbrk[3])
  legend("bottomleft",
         legend = c(
           sprintf("Large (deg≥%.0f, n=%d)", dbrk[3], n_lg),
           sprintf("Medium (deg %.0f–%.0f, n=%d)", dbrk[2], dbrk[3]-1, n_md),
           sprintf("Small (deg<%.0f, n=%d)", dbrk[2], n_sm),
           "Metabolite", "Microbe",
           sprintf("Pos (n=%d)", n_pos), sprintf("Neg (n=%d)", n_neg)),
         pch = c(21, 21, 21, 21, 21, NA, NA),
         pt.bg = c("grey60", "grey60", "grey60", METAB_COL, MICROBE_COL, NA, NA),
         col = c("grey40", "grey40", "grey40", "grey40", "grey40", POS_COLOR, NEG_COLOR),
         pt.cex = c(4.0, 2.4, 1.2, 1.0, 1.0, NA, NA),
         lty = c(NA, NA, NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, NA, NA, 1.2, 1.2),
         cex = 0.55, bty = "n")

  dev.off()
  cat(sprintf("  -> %s_microbe_metabolite_network.svg\n", soil_label))

  # ---- 7. 导出 ----
  edge_list <- edge_list[order(edge_list$P_value), ]
  write.csv(edge_list, sprintf("output/%s_microbe_metabolite_edges.csv", soil_label), row.names = FALSE)
  node_df <- data.frame(Name = V(g)$name, Type = V(g)$type, Degree = degree(g),
                         Betweenness = round(betweenness(g, normalized = TRUE), 4))
  node_df <- node_df[order(node_df$Degree, decreasing = TRUE), ]
  write.csv(node_df, sprintf("output/%s_microbe_metabolite_nodes.csv", soil_label), row.names = FALSE)

  list(n_sig = n_sig, edge_list = edge_list, graph = g)
}

# ================================================================
cat("============================================================\n")
cat("SS & NSS Microbe-Metabolite Correlation\n")
cat("============================================================\n")

ss_cols  <- list(CK=c("W1_1","W1_2","W1_3"), S=c("W2_1","W2_2","W2_3"),
                 SC1=c("W3_1","W3_2","W3_3"), SC5=c("W7_1","W7_2","W7_3"))
nss_cols <- list(CK=c("W8_1","W8_2","W8_3"), S=c("W9_1","W9_2","W9_3"),
                 SC1=c("W10_1","W10_2","W10_3"), SC5=c("W14_1","W14_2","W14_3"))

ss_res  <- run_analysis("SS", ss_cols, nss_cols)
nss_res <- run_analysis("NSS", ss_cols, nss_cols)

cat(sprintf("\nFinal: SS=%d sig, NSS=%d sig\n", ss_res$n_sig, nss_res$n_sig))
cat("DONE!\n")
