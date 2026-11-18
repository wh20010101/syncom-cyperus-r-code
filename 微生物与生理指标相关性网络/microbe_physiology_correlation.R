# ================================================================
# 微生物-生理指标 Spearman 相关性分析
# 筛选：检出率≥50% + 丰度>0.01% → top100属
# SS/NSS 分开分析 + Raw P (不矫正) + 热图
# ================================================================

library(readxl)
library(igraph)
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

cat("============================================================\n")
cat("Microbe-Physiology Spearman Correlation\n")
cat("============================================================\n")

# ---------- 1. 读取数据 ----------
cat("Step 1: Reading data...\n")

# 微生物
df <- read.delim("../data/微生物OUT.xls", sep = "\t", header = TRUE, check.names = FALSE)
data_cols <- grep("^W[0-9]+_[0-9]+$", names(df), value = TRUE)
df$genus_clean <- gsub("^g__", "", df$genus)
genus_abund <- aggregate(df[, data_cols], by = list(genus = df$genus_clean), FUN = sum)
rownames(genus_abund) <- genus_abund$genus; genus_abund$genus <- NULL

# 生理指标
phys_raw <- read_excel("生理指标.xlsx")
phys_mat <- as.data.frame(phys_raw)
rownames(phys_mat) <- phys_mat$Sample
phys_mat$Sample <- NULL

# 分组
group_raw <- read_excel("分组信息.xlsx")
stopifnot(all(c("Sample", "Group") %in% names(group_raw)))

expected_groups <- c(
  "SS-CK", "SS-S", "SS-SC1", "SS-SC5",
  "NSS-CK", "NSS-S", "NSS-SC1", "NSS-SC5"
)
missing_groups <- setdiff(expected_groups, unique(as.character(group_raw$Group)))
if (length(missing_groups) > 0) {
  stop("分组信息中的 Group 缺少或拼写不一致: ",
       paste(missing_groups, collapse = ", "),
       "\n请使用 SS-CK/SS-S/SS-SC1/SS-SC5 与 NSS-... 的标准写法。")
}

# Sample → W column mapping
w_cols <- c(
  "W1_1","W1_2","W1_3",   # SS-CK     (1-3)
  "W2_1","W2_2","W2_3",   # SS-S      (4-6)
  "W3_1","W3_2","W3_3",   # SS-SC1    (7-9)
  "W7_1","W7_2","W7_3",   # SS-SC5    (10-12)
  "W8_1","W8_2","W8_3",   # NSS-CK    (13-15)
  "W9_1","W9_2","W9_3",   # NSS-S     (16-18)
  "W10_1","W10_2","W10_3",# NSS-SC1   (19-21)
  "W14_1","W14_2","W14_3" # NSS-SC5   (22-24)
)
sample_map <- setNames(w_cols, 1:24)
cat(sprintf("  Microbes: %d genera, Physiology: %d indicators, Samples: %d\n",
            nrow(genus_abund), ncol(phys_mat), length(sample_map)))

# ---------- 2. 属筛选 ----------
cat("\nStep 2: Filtering genera...\n")

# Remove unclassified
is_uncl <- grepl("unclassified", rownames(genus_abund), ignore.case = TRUE)
genus_clean <- genus_abund[!is_uncl, ]
cat(sprintf("  After removing unclassified: %d\n", nrow(genus_clean)))

# Reorder columns to match Sample 1-24
genus_ordered <- genus_clean[, w_cols]
colnames(genus_ordered) <- 1:24

# Relative abundance
rel_abund <- sweep(genus_ordered, 2, colSums(genus_ordered), "/") * 100

# Filter 1: prevalence >= 50%
detection_rate <- rowSums(genus_ordered > 0) / ncol(genus_ordered)
keep1 <- detection_rate >= 0.5
cat(sprintf("  Prevalence >= 50%%: %d\n", sum(keep1)))

# Filter 2: mean abundance > 0.01%
mean_abund <- rowMeans(rel_abund)
keep2 <- mean_abund > 0.01
cat(sprintf("  Mean abund > 0.01%%: %d\n", sum(keep2)))

# Combined
keep_both <- keep1 & keep2
genus_filtered <- genus_ordered[keep_both, ]
cat(sprintf("  Combined: %d genera\n", nrow(genus_filtered)))

# Top 100 by mean abundance
mean_filtered <- sort(rowMeans(genus_filtered), decreasing = TRUE)
top100 <- names(head(mean_filtered, 100))
genus_top <- genus_filtered[top100, ]
cat(sprintf("  Top 100 selected\n"))

# ---------- 3. SS / NSS 分开 ----------
cat("\nStep 3: Separating SS and NSS...\n")

ss_samples <- as.character(
  group_raw$Sample[group_raw$Group %in% c("SS-CK", "SS-S", "SS-SC1", "SS-SC5")]
)
nss_samples <- as.character(
  group_raw$Sample[group_raw$Group %in% c("NSS-CK", "NSS-S", "NSS-SC1", "NSS-SC5")]
)

stopifnot(
  all(ss_samples %in% rownames(phys_mat)),
  all(nss_samples %in% rownames(phys_mat)),
  all(ss_samples %in% colnames(genus_ordered)),
  all(nss_samples %in% colnames(genus_ordered))
)

# Clean genus names for display
clean_names <- function(x) {
  x <- gsub("^norank_f__|^norank_o__|^norank_c__|^norank_p__|^norank_", "", x)
  x <- gsub("^f__|^o__|^c__|^p__", "", x)
  x <- gsub("__.*$", "", x)
  gsub("_", " ", x)
}

run_correlation <- function(soil_label, sample_ids) {
  cat(sprintf("\n--- %s (n=%d) ---\n", soil_label, length(sample_ids)))

  # Microbe data for these samples
  microbe_mat <- genus_top[, sample_ids, drop = FALSE]
  microbe_mat <- microbe_mat[rowSums(microbe_mat) > 0, ]  # remove all-zero in subset
  cat(sprintf("  Genera present: %d\n", nrow(microbe_mat)))

  # Physiology data
  phys_sub <- phys_mat[sample_ids, , drop = FALSE]

  # Spearman correlation matrix: genera × physiology
  n_genus <- nrow(microbe_mat)
  n_phys <- ncol(phys_sub)
  r_mat <- matrix(NA, n_genus, n_phys)
  p_mat <- matrix(NA, n_genus, n_phys)

  for (i in 1:n_genus) {
    for (j in 1:n_phys) {
      test <- cor.test(as.numeric(microbe_mat[i, ]), phys_sub[[j]],
                       method = "spearman", exact = FALSE)
      r_mat[i, j] <- test$estimate
      p_mat[i, j] <- test$p.value
    }
  }
  rownames(r_mat) <- clean_names(rownames(microbe_mat))
  colnames(r_mat) <- colnames(phys_sub)
  rownames(p_mat) <- clean_names(rownames(microbe_mat))
  colnames(p_mat) <- colnames(phys_sub)

  # Raw P (no correction)
  p_raw_mat <- p_mat

  # Count significant edges (same criteria as the network: P + |r|)
  n_sig <- sum(p_raw_mat < 0.05 & abs(r_mat) >= 0.7, na.rm = TRUE)
  n_sig_01 <- sum(p_raw_mat < 0.01 & abs(r_mat) >= 0.7, na.rm = TRUE)
  cat(sprintf("  Significant edges (Raw P<0.05, |r|>=0.7): %d; (Raw P<0.01, |r|>=0.7): %d\n", n_sig, n_sig_01))

  list(r = r_mat, p = p_mat, p_raw = p_raw_mat, n_sig = n_sig)
}

ss_cor <- run_correlation("SS", ss_samples)
nss_cor <- run_correlation("NSS", nss_samples)

# ---------- 3.5. 相关性热图 (Raw P) ----------
cat("\nStep 3.5: Plotting correlation heatmaps (Raw P)...\n")
library(pheatmap)

plot_heatmap <- function(cor_obj, soil_label) {
  r_mat <- cor_obj$r
  p_raw_mat <- cor_obj$p_raw

  # 筛选至少有一个Raw P显著相关的菌属和指标
  has_sig_genus <- rowSums(p_raw_mat < 0.05, na.rm = TRUE) > 0
  has_sig_phys  <- colSums(p_raw_mat < 0.05, na.rm = TRUE) > 0

  if (sum(has_sig_genus) < 2 || sum(has_sig_phys) < 2) {
    cat(sprintf("  %s: not enough significant pairs for heatmap\n", soil_label))
    return()
  }

  r_plot <- r_mat[has_sig_genus, has_sig_phys, drop = FALSE]
  p_plot <- p_raw_mat[has_sig_genus, has_sig_phys, drop = FALSE]

  # 显著性标记 (based on Raw P)
  sig_mark <- matrix("", nrow(r_plot), ncol(r_plot))
  sig_mark[p_plot < 0.05]  <- "*"
  sig_mark[p_plot < 0.01]  <- "**"
  sig_mark[p_plot < 0.001] <- "***"

  nrow_h <- nrow(r_plot)
  ncol_w <- ncol(r_plot)
  fig_h <- max(100, nrow_h * 12)
  fig_w <- max(200, ncol_w * 18)

  svglite(sprintf("output/%s_correlation_heatmap.svg", soil_label),
          width = fig_w * MM, height = fig_h * MM)

  pheatmap(r_plot,
           cluster_rows = TRUE, cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = sig_mark,
           number_color = "black",
           fontsize_number = 8, fontsize = 8,
           fontsize_col = 9, fontsize_row = 7,
           main = sprintf("%s: Genus–Physiology Spearman (Raw P<0.05)\n%d genera × %d indicators",
                          soil_label, nrow(r_plot), ncol(r_plot)),
           angle_col = 45,
           border_color = NA,
           cellwidth = 16, cellheight = 10,
           legend_breaks = c(-1, -0.5, 0, 0.5, 1),
           legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
  dev.off()
  cat(sprintf("  -> %s_correlation_heatmap.svg (%d genera × %d indicators)\n",
              soil_label, nrow(r_plot), ncol(r_plot)))
}

plot_heatmap(ss_cor, "SS")
plot_heatmap(nss_cor, "NSS")

# ---------- 4. 相关性网络 ----------
cat("\nStep 4: Building correlation networks...\n")

# 边颜色
POS_COLOR <- "#D95F02"   # 正相关
NEG_COLOR <- "#4575B4"   # 负相关

# 节点颜色：菌=青色，生理指标=橙色
GENUS_COL  <- "#4DBBD5"
PHYS_COL   <- "#E64B35"

build_cor_network <- function(cor_obj, soil_label, p_cutoff = 0.05, r_cutoff = 0.7) {
  r_mat <- cor_obj$r
  p_raw_mat <- cor_obj$p_raw

  # 提取显著边（P + |r| 双筛选）
  edge_list <- data.frame()
  for (i in 1:nrow(r_mat)) {
    for (j in 1:ncol(r_mat)) {
      if (!is.na(p_raw_mat[i, j]) && p_raw_mat[i, j] < p_cutoff &&
          !is.na(r_mat[i, j]) && abs(r_mat[i, j]) >= r_cutoff) {
        edge_list <- rbind(edge_list, data.frame(
          Node1 = rownames(r_mat)[i],
          Node2 = colnames(r_mat)[j],
          Coefficient = round(r_mat[i, j], 3),
          P_value = signif(p_raw_mat[i, j], 3),
          Direction = ifelse(r_mat[i, j] > 0, "positive", "negative"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(edge_list) == 0) {
    cat(sprintf("  %s: no edges at P<%.2f & |r|≥%.1f\n", soil_label, p_cutoff, r_cutoff))
    return(NULL)
  }

  # 构建 igraph
  all_nodes <- unique(c(edge_list$Node1, edge_list$Node2))
  g <- make_empty_graph(n = length(all_nodes), directed = FALSE)
  V(g)$name <- all_nodes
  V(g)$type <- ifelse(all_nodes %in% rownames(r_mat), "Genus", "Physiology")
  V(g)$color <- ifelse(V(g)$type == "Genus", GENUS_COL, PHYS_COL)

  for (k in 1:nrow(edge_list)) {
    g <- add_edges(g, c(edge_list$Node1[k], edge_list$Node2[k]),
                   weight = abs(edge_list$Coefficient[k]),
                   coefficient = edge_list$Coefficient[k],
                   sign = edge_list$Direction[k])
  }

  cat(sprintf("  %s: %d nodes, %d edges (pos=%d, neg=%d)\n",
              soil_label, vcount(g), ecount(g),
              sum(E(g)$coefficient > 0), sum(E(g)$coefficient < 0)))

  list(graph = g, edge_list = edge_list, soil = soil_label)
}

ss_net <- build_cor_network(ss_cor, "SS")
nss_net <- build_cor_network(nss_cor, "NSS")

# ---------- 5. 绘图 ----------
cat("\nStep 5: Plotting networks...\n")

plot_network <- function(net, soil_label) {
  if (is.null(net)) return()
  g <- net$graph
  edge_list <- net$edge_list

  # 分离菌属和生理指标
  genuses   <- sort(unique(edge_list$Node1))
  phys_inds <- sort(unique(edge_list$Node2))
  n_g <- length(genuses)
  n_p <- length(phys_inds)

  # ---- 双圈坐标 ----
  r_outer <- 1.0
  r_inner <- 0.65

  angles_g <- seq(0, 2 * pi, length.out = n_g + 1)[1:n_g]
  coords_g <- data.frame(name = genuses, x = r_outer * cos(angles_g),
                         y = r_outer * sin(angles_g), type = "Genus",
                         stringsAsFactors = FALSE)

  angles_p <- seq(0, 2 * pi, length.out = n_p + 1)[1:n_p]
  coords_p <- data.frame(name = phys_inds, x = r_inner * cos(angles_p),
                         y = r_inner * sin(angles_p), type = "Physiology",
                         stringsAsFactors = FALSE)

  all_coords <- rbind(
    data.frame(name = coords_g$name, x = coords_g$x, y = coords_g$y, stringsAsFactors = FALSE),
    data.frame(name = coords_p$name, x = coords_p$x, y = coords_p$y, stringsAsFactors = FALSE)
  )
  rownames(all_coords) <- all_coords$name

  # ---- 节点度 → 大小 ----
  degs <- degree(g)
  deg_g <- degs[coords_g$name]
  deg_p <- degs[coords_p$name]
  deg_all <- c(deg_g, deg_p)
  deg_breaks <- quantile(deg_all, probs = c(0, 0.5, 0.85, 1))
  size_cat <- function(d) {
    if (d >= deg_breaks[3]) "Large"
    else if (d >= deg_breaks[2]) "Medium"
    else "Small"
  }
  size_map <- c("Small" = 1.2, "Medium" = 2.4, "Large" = 4.0)
  coords_g$cex <- size_map[sapply(deg_g, size_cat)]
  coords_p$cex <- size_map[sapply(deg_p, size_cat)]

  # ---- 绘图 ----
  n_pos <- sum(edge_list$Direction == "positive")
  n_neg <- sum(edge_list$Direction == "negative")
  svglite(sprintf("output/%s_correlation_network.svg", soil_label),
          width = 100 * MM, height = 100 * MM,
          pointsize = 8, system_fonts = list(serif = "Times New Roman"))

  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white", family = "serif")
  plot(0, 0, type = "n", xlim = c(-1.35, 1.35), ylim = c(-1.25, 1.25),
       xlab = "", ylab = "", axes = FALSE, asp = 1)

  # 边
  for (k in 1:nrow(edge_list)) {
    n1 <- edge_list$Node1[k]; n2 <- edge_list$Node2[k]
    if (!(n1 %in% rownames(all_coords) && n2 %in% rownames(all_coords))) next
    x1 <- all_coords[n1, "x"]; y1 <- all_coords[n1, "y"]
    x2 <- all_coords[n2, "x"]; y2 <- all_coords[n2, "y"]
    col <- ifelse(edge_list$Direction[k] == "positive", POS_COLOR, NEG_COLOR)
    alpha_val <- pmax(0.15, abs(edge_list$Coefficient[k]) * 0.5)
    rgb_vals <- col2rgb(col)[, 1] / 255
    edge_col <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha_val)
    lines(c(x1, x2), c(y1, y2), col = edge_col, lwd = abs(edge_list$Coefficient[k]) * 1.2)
  }

  # 外圈菌属
  points(coords_g$x, coords_g$y, pch = 21, bg = adjustcolor(GENUS_COL, 0.85),
         col = "grey40", cex = coords_g$cex)
  # 内圈指标
  points(coords_p$x, coords_p$y, pch = 21, bg = adjustcolor(PHYS_COL, 0.85),
         col = "grey40", cex = coords_p$cex)

  # 菌属标签（外圈圆心）
  for (i in 1:n_g) {
    ang <- angles_g[i]; tx <- r_outer * cos(ang); ty <- r_outer * sin(ang)
    text(tx, ty, labels = coords_g$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }
  # 指标标签（内圈圆心）
  for (i in 1:n_p) {
    ang <- angles_p[i]; tx <- r_inner * cos(ang); ty <- r_inner * sin(ang)
    text(tx, ty, labels = coords_p$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }

  # 图例
  n_sm <- sum(sapply(c(deg_g, deg_p), size_cat) == "Small")
  n_md <- sum(sapply(c(deg_g, deg_p), size_cat) == "Medium")
  n_lg <- sum(sapply(c(deg_g, deg_p), size_cat) == "Large")
  legend("bottomleft",
         legend = c(
           sprintf("Large (deg≥%.0f, n=%d)", deg_breaks[3], n_lg),
           sprintf("Medium (deg %.0f–%.0f, n=%d)", deg_breaks[2], deg_breaks[3]-1, n_md),
           sprintf("Small (deg<%.0f, n=%d)", deg_breaks[2], n_sm),
           "Genus", "Physiology",
           sprintf("Pos (n=%d)", n_pos), sprintf("Neg (n=%d)", n_neg)),
         pch = c(21, 21, 21, 21, 21, NA, NA),
         pt.bg = c("grey60", "grey60", "grey60", GENUS_COL, PHYS_COL, NA, NA),
         col = c("grey40", "grey40", "grey40", "grey40", "grey40", POS_COLOR, NEG_COLOR),
         pt.cex = c(4.0, 2.4, 1.2, 1.0, 1.0, NA, NA),
         lty = c(NA, NA, NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, NA, NA, 1.2, 1.2),
         cex = 0.55, bty = "n")

  dev.off()
  cat(sprintf("  -> %s_correlation_network.svg\n", soil_label))
}

plot_network(ss_net, "SS")
plot_network(nss_net, "NSS")

# ---------- 6. 导出 ----------
cat("\nStep 6: Exporting...\n")

for (net in list(ss_net, nss_net)) {
  if (is.null(net)) next
  soil_label <- net$soil
  if (nrow(net$edge_list) > 0) {
    write.csv(net$edge_list,
              sprintf("output/%s_correlation_edges.csv", soil_label),
              row.names = FALSE)
  }
  # Node attributes
  g <- net$graph
  node_df <- data.frame(
    Name = V(g)$name,
    Type = V(g)$type,
    Degree = degree(g),
    Betweenness = round(betweenness(g, normalized = TRUE), 4),
    stringsAsFactors = FALSE
  )
  node_df <- node_df[order(node_df$Degree, decreasing = TRUE), ]
  write.csv(node_df, sprintf("output/%s_correlation_nodes.csv", soil_label),
            row.names = FALSE)
  cat(sprintf("  %s: edges=%d, nodes=%d\n", soil_label, nrow(net$edge_list), nrow(node_df)))
}

# ---------- 7. 汇总 ----------
cat("\n============================================================\n")
cat(sprintf("Summary:\n"))
cat(sprintf("  SS:  %d genera x %d indicators, %d sig edges (Raw P<0.05, |r|≥0.7)\n",
            nrow(ss_cor$r), ncol(ss_cor$r), ss_cor$n_sig))
cat(sprintf("  NSS: %d genera x %d indicators, %d sig edges (Raw P<0.05, |r|≥0.7)\n",
            nrow(nss_cor$r), ncol(nss_cor$r), nss_cor$n_sig))
cat(sprintf("Output: %s/output/\n", getwd()))
cat("DONE!\n")
