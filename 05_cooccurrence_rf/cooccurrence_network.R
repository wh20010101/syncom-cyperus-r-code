# ================================================================
# 共现网络分析 — 4个处理组（CK, S, SC1, SC5），合并灭菌+非灭菌土
# Spearman |ρ| > 0.7, FDR校正后 P<0.05
# 前100属（排除unclassified），节点颜色按门，节点大小按丰度
# 输出 SVG 图 + 边列表 + 节点属性 + 拓扑对比表
# ================================================================

# ---------- 0. 加载/安装包 ----------
packages <- c("igraph", "Hmisc", "dplyr", "tidyr", "readxl", "svglite", "writexl")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org")
  library(pkg, character.only = TRUE)
}

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
dir.create("output_r", showWarnings = FALSE)
outdir <- "output_r"

cat("============================================================\n")
cat("Step 1: Reading data...\n")

# ---------- 1. 读取数据 ----------
df <- read.delim("../data/microbial_OTU_table.xls", sep = "\t", header = TRUE, check.names = FALSE)
data_cols <- grep("^W\\d+_\\d+$", names(df), value = TRUE)
cat(sprintf("  Samples: %d\n", length(data_cols)))

# ---------- 2. 属水平聚合 ----------
df$genus_clean <- gsub("^g__", "", df$genus)

genus_abund <- aggregate(df[, data_cols], by = list(genus = df$genus_clean), FUN = sum)
rownames(genus_abund) <- genus_abund$genus
genus_abund$genus <- NULL

# 属→门（众数）
phylum_df <- aggregate(list(phylum = df$phylum), by = list(genus = df$genus_clean),
                        FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
phylum_vec <- gsub("^p__", "", phylum_df$phylum)
names(phylum_vec) <- phylum_df$genus

cat(sprintf("  Genera: %d\n", nrow(genus_abund)))

# ---------- 3. 前100属：排除 unclassified，取整体平均丰度前100 ----------
cat("\nStep 3: Selecting top 100 genera (excluding unclassified)...\n")
is_uncl <- grepl("unclassified", rownames(genus_abund), ignore.case = TRUE)
cat(sprintf("  Removing %d unclassified genera\n", sum(is_uncl)))

genus_abund_clean <- genus_abund[!is_uncl, ]
overall_mean <- sort(rowMeans(genus_abund_clean), decreasing = TRUE)
N_TOP <- 100
topN <- names(overall_mean)[1:min(N_TOP, length(overall_mean))]
topN_phylum <- phylum_vec[topN]

cat(sprintf("  Selected %d genera from %d phyla\n", length(topN), length(unique(topN_phylum))))
print(table(topN_phylum))

# ---------- 4. 定义处理组 ----------
cat("\nStep 4: Building networks (|r| > 0.7, FDR P < 0.05)...\n")
CUTOFF <- 0.7

groups <- list(
  CK  = c("W1_1", "W1_2", "W1_3", "W8_1", "W8_2", "W8_3"),
  S   = c("W2_1", "W2_2", "W2_3", "W9_1", "W9_2", "W9_3"),
  SC1 = c("W3_1", "W3_2", "W3_3", "W10_1", "W10_2", "W10_3"),
  SC5 = c("W7_1", "W7_2", "W7_3", "W14_1", "W14_2", "W14_3")
)

# ---------- 网络构建（含 FDR 矫正）----------
build_network <- function(abund_mat, genera, group_name, cutoff) {
  dat <- t(abund_mat[genera, , drop = FALSE])

  cors <- rcorr(dat, type = "spearman")
  r_mat <- cors$r
  p_mat <- cors$P

  # FDR 矫正（仅下三角）
  p_vec <- p_mat[lower.tri(p_mat)]
  p_adj <- p.adjust(p_vec, method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat))
  p_adj_mat[lower.tri(p_adj_mat)] <- p_adj
  p_adj_mat[upper.tri(p_adj_mat)] <- t(p_adj_mat)[upper.tri(t(p_adj_mat))]

  g <- make_empty_graph(n = length(genera), directed = FALSE)
  V(g)$name <- genera
  V(g)$phylum <- topN_phylum[genera]
  V(g)$mean_abund <- rowMeans(abund_mat[genera, , drop = FALSE])

  n_genus <- length(genera)
  edge_list <- data.frame()
  for (i in 1:(n_genus - 1)) {
    for (j in (i + 1):n_genus) {
      r_val <- r_mat[i, j]
      p_fdr <- p_adj_mat[i, j]
      if (!is.na(r_val) && !is.na(p_fdr) &&
          abs(r_val) > cutoff && p_fdr < 0.05) {
        g <- add_edges(g, c(genera[i], genera[j]),
                       weight = abs(r_val),
                       coefficient = r_val,
                       pvalue = p_mat[i, j],
                       pvalue_fdr = p_fdr,
                       sign = ifelse(r_val > 0, "+", "-"))
        edge_list <- rbind(edge_list, data.frame(
          Node1 = genera[i],
          Node2 = genera[j],
          Coefficient = round(r_val, 4),
          Pvalue = round(p_mat[i, j], 6),
          Pvalue_FDR = round(p_fdr, 6),
          Direction = ifelse(r_val > 0, "positive", "negative"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  list(graph = g, edge_list = edge_list,
       n_nodes = vcount(g), n_edges = ecount(g),
       n_pos = if (ecount(g) > 0) sum(E(g)$coefficient > 0) else 0,
       n_neg = if (ecount(g) > 0) sum(E(g)$coefficient < 0) else 0)
}

# ---------- 建网 ----------
networks <- list()
for (grp_name in names(groups)) {
  cols <- groups[[grp_name]]
  cat(sprintf("\n  --- %s (n=%d) ---\n", grp_name, length(cols)))

  grp_abund <- genus_abund_clean[topN, cols, drop = FALSE]
  net <- build_network(grp_abund, topN, grp_name, cutoff = CUTOFF)

  deg <- degree(net$graph)
  net$avg_degree <- if (length(deg) > 0) mean(deg) else 0
  net$density <- edge_density(net$graph)
  comps <- components(net$graph)
  net$n_components <- comps$no
  net$largest_component <- max(comps$csize)
  networks[[grp_name]] <- net

  cat(sprintf("    Edges: %d (pos=%d, neg=%d)\n", net$n_edges, net$n_pos, net$n_neg))
  cat(sprintf("    Topology: avg_deg=%.2f, density=%.4f, components=%d\n",
              net$avg_degree, net$density, net$n_components))
}

# ================================================================
# ---------- 5. 绘图参数 ----------
cat("\nStep 5: Plotting networks...\n")

# 门 → 颜色（20色，确保18个门各不同）
phyla_list <- unique(unlist(lapply(networks, function(n) V(n$graph)$phylum)))
phylum_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",  # 1-5
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",  # 6-10
  "#FFA500", "#6A3D9A", "#8B4513", "#A0A0A0", "#228B22",  # 11-15
  "#CD5C5C", "#4682B4", "#DAA520", "#8FBC8F", "#B8860B"   # 16-20
)
phylum_colors <- setNames(phylum_palette[1:length(phyla_list)], phyla_list)

# 边颜色（新配色）
# Positive → 暖橙红  Negative → 冷钢蓝
POS_COLOR <- "#D95F02"   # 橙红，正相关
NEG_COLOR <- "#4575B4"   # 钢蓝，负相关

# 标题
title_map <- c(
  CK  = "CK (SS-CK + NSS-CK, n=6)",
  S   = "Salt Stress (SS-S + NSS-S, n=6)",
  SC1 = "SynCom1 (SS-SC1 + NSS-SC1, n=6)",
  SC5 = "SynCom5 (SS-SC5 + NSS-SC5, n=6)"
)

# 布局函数
get_layout <- function(g) {
  if (ecount(g) >= 2) {
    tryCatch(layout_with_kk(g), error = function(e) layout_with_fr(g, niter = 800))
  } else if (ecount(g) == 1) {
    layout_with_fr(g, niter = 500)
  } else {
    layout_in_circle(g)
  }
}

# ================================================================
# ---- 单图绘图 (100mm宽, Times New Roman 8pt) ----
# ================================================================
MM_TO_INCH <- 1 / 25.4
FIG_W <- 100 * MM_TO_INCH   # 100mm
FIG_H <- 95 * MM_TO_INCH

# 度分档 & 大小映射
size_bins <- list(
  "0-2" = list(lo=0, hi=2, vsize=4),
  "3-4" = list(lo=3, hi=4, vsize=8),
  "5+"  = list(lo=5, hi=999, vsize=14)
)

for (grp_name in names(networks)) {
  cat(sprintf("  Plotting %s...\n", grp_name))
  net <- networks[[grp_name]]
  g <- net$graph

  # ---- 节点颜色 ----
  node_cols <- phylum_colors[V(g)$phylum]
  node_cols[is.na(node_cols)] <- "#CCCCCC"

  # ---- 节点大小：3档 ----
  degs <- degree(g)
  node_size <- rep(4, length(degs))  # default small
  bin_labels <- character(length(degs))
  for (bin_name in names(size_bins)) {
    b <- size_bins[[bin_name]]
    idx <- degs >= b$lo & degs <= b$hi
    node_size[idx] <- b$vsize
    bin_labels[idx] <- bin_name
  }

  # 统计每档节点数
  bin_counts <- sapply(names(size_bins), function(bn) {
    b <- size_bins[[bn]]
    sum(degs >= b$lo & degs <= b$hi)
  })

  # ---- 边 ----
  if (ecount(g) > 0) {
    edge_coef <- E(g)$coefficient
    edge_base_color <- ifelse(edge_coef > 0, POS_COLOR, NEG_COLOR)
    edge_width <- abs(edge_coef) * 2.5
    edge_col_final <- character(length(edge_coef))
    for (k in seq_along(edge_coef)) {
      rgb_vals <- col2rgb(edge_base_color[k])[, 1] / 255
      edge_col_final[k] <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3],
                               alpha = pmax(0.3, abs(edge_coef[k]) * 0.8))
    }
  } else {
    edge_col_final <- character(0)
    edge_width <- numeric(0)
  }

  # ---- 布局 ----
  lay <- get_layout(g)

  # ---- 标签：仅标度最高的前8个属 ----
  deg_sorted <- names(sort(degs, decreasing = TRUE))
  top8_deg <- head(deg_sorted[degs[deg_sorted] > 0], 8)

  vertex_labels <- V(g)$name
  vertex_labels[!vertex_labels %in% top8_deg] <- ""
  # 简化属名：去 norank_ 前缀，只留属名本体
  vertex_labels <- gsub("^norank_f__|^norank_o__|^norank_c__|^norank_p__|^norank_", "", vertex_labels)
  vertex_labels <- gsub("^f__|^o__|^c__|^p__", "", vertex_labels)
  vertex_labels <- gsub("__.*$", "", vertex_labels)  # 去掉多余后缀
  vertex_labels <- gsub("_", " ", vertex_labels)

  # ---- SVG ----
  svglite(sprintf("%s/%s_network.svg", outdir, grp_name),
      width = FIG_W, height = FIG_H)

  # 布局：主图左 + 图例右
  layout_matrix <- rbind(c(1, 1, 2), c(1, 1, 3), c(1, 1, 4))
  layout(layout_matrix, widths = c(1, 1, 0.7))

  # --- 主图 ---
  par(mar = c(0.5, 0.5, 2, 0.5))
  plot(g,
       layout = lay,
       vertex.label = vertex_labels,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.size = node_size,
       vertex.color = adjustcolor(node_cols, alpha.f = 0.85),
       vertex.frame.color = "white",
       vertex.frame.width = 0.3,
       edge.color = edge_col_final,
       edge.width = edge_width,
       edge.curved = 0.1,
       margin = c(-0.05, -0.05, -0.05, -0.05))

  # 标题
  title(main = sprintf("%s  (|r| > %.1f, FDR P < 0.05)", title_map[grp_name], CUTOFF),
        cex.main = 0.75, line = 0.5)

  # --- 图例面板1: 节点大小（用 
  par(mar = c(0, 0.2, 0, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # 画3个参考圆
  ref_x <- 0.15
  ref_y <- c(0.85, 0.65, 0.45)
  ref_sizes <- c(4, 8, 14)
  for (i in 1:3) {
    points(ref_x, ref_y[i], cex = ref_sizes[i] / 3, pch = 21,
           bg = "gray75", col = "gray50")
  }

  # 文字标注
  for (i in 1:3) {
    bn <- names(size_bins)[i]
    txt <- sprintf("Deg %s (n=%d)", bn, bin_counts[bn])
    text(0.15, ref_y[i] - 0.08, txt, cex = 0.55,
         col = "gray30", adj = c(0.5, 1))
  }

  # 图例标题
  text(0.15, 0.95, "Node Size\n(Degree)", cex = 0.6,
       font = 2, col = "black", adj = c(0.5, 0))

  # --- 图例面板2: 门 ---
  par(mar = c(0, 0.2, 0, 0.2))
  plot.new()
  phyla_local <- names(sort(table(V(g)$phylum), decreasing = TRUE))
  cols_local <- phylum_colors[phyla_local]
  cols_local[is.na(cols_local)] <- "#CCCCCC"
  n_phyla <- length(phyla_local)
  y_pos <- seq(0.95, 0.05, length.out = n_phyla + 2)
  y_pos <- y_pos[2:(n_phyla + 1)]
  for (i in 1:n_phyla) {
    
    text(0.22, y_pos[i], sprintf("%s (%d)", phyla_local[i],
         sum(V(g)$phylum == phyla_local[i])),
         cex = 0.5, adj = c(0, 0.5))
  }
  text(0.1, 0.97, "Phylum", cex = 0.6, font = 2, adj = c(0.5, 0))

  # --- 图例面板3: 边 ---
  par(mar = c(0, 0.2, 0, 0.2))
  plot.new()
  segments(0.05, 0.8, 0.3, 0.8, col = POS_COLOR, lwd = 2)
  text(0.35, 0.8, sprintf("Positive (n=%d)", net$n_pos),
       cex = 0.55, adj = c(0, 0.5))
  segments(0.05, 0.55, 0.3, 0.55, col = NEG_COLOR, lwd = 2)
  text(0.35, 0.55, sprintf("Negative (n=%d)", net$n_neg),
       cex = 0.55, adj = c(0, 0.5))
  text(0.15, 0.95, "Edge", cex = 0.6, font = 2, adj = c(0.5, 0))

  dev.off()
  cat(sprintf("    -> %s/%s_network.svg\n", outdir, grp_name))
}

# ================================================================
# ========== 6. 四图并排 (180mm宽, Times New Roman 8pt) ==========
# ================================================================
cat("\nStep 6: Combined 4-panel SVG...\n")

svglite(sprintf("%s/combined_4_networks.svg", outdir),
    width = 180 * MM_TO_INCH, height = 180 * MM_TO_INCH)

layout_matrix <- rbind(c(1, 2), c(3, 4), c(5, 5))
layout(layout_matrix, heights = c(1, 1, 0.15))

size_bins_comb <- list(
  "0-2" = list(lo=0, hi=2, vsize=3.5),
  "3-4" = list(lo=3, hi=4, vsize=7),
  "5+"  = list(lo=5, hi=999, vsize=12)
)

for (grp_name in names(networks)) {
  net <- networks[[grp_name]]
  g <- net$graph

  node_cols <- phylum_colors[V(g)$phylum]
  node_cols[is.na(node_cols)] <- "#CCCCCC"

  degs <- degree(g)
  node_size <- rep(3.5, length(degs))
  for (bin_name in names(size_bins_comb)) {
    b <- size_bins_comb[[bin_name]]
    idx <- degs >= b$lo & degs <= b$hi
    node_size[idx] <- b$vsize
  }

  if (ecount(g) > 0) {
    edge_coef <- E(g)$coefficient
    edge_base_color <- ifelse(edge_coef > 0, POS_COLOR, NEG_COLOR)
    edge_width <- abs(edge_coef) * 1.8
    edge_col_final <- character(length(edge_coef))
    for (k in seq_along(edge_coef)) {
      rgb_vals <- col2rgb(edge_base_color[k])[, 1] / 255
      edge_col_final[k] <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3],
                               alpha = pmax(0.3, abs(edge_coef[k]) * 0.6))
    }
  } else {
    edge_col_final <- character(0)
    edge_width <- numeric(0)
  }

  lay <- get_layout(g)
  deg_sorted <- names(sort(degs, decreasing = TRUE))
  top5_deg <- head(deg_sorted[degs[deg_sorted] > 0], 5)

  vertex_labels <- V(g)$name
  vertex_labels[!vertex_labels %in% top5_deg] <- ""
  vertex_labels <- gsub("^norank_f__|^norank_o__|^norank_c__|^norank_p__|^norank_", "", vertex_labels)
  vertex_labels <- gsub("^f__|^o__|^c__|^p__", "", vertex_labels)
  vertex_labels <- gsub("__.*$", "", vertex_labels)

  par(mar = c(0.2, 0.2, 2, 0.2))
  plot(g, layout = lay,
       vertex.label = vertex_labels,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.size = node_size,
       vertex.color = adjustcolor(node_cols, alpha.f = 0.85),
       vertex.frame.color = "white",
       vertex.frame.width = 0.2,
       edge.color = edge_col_final,
       edge.width = edge_width,
       edge.curved = 0.08,
       margin = c(-0.05, -0.05, -0.05, -0.05))

  title(main = sprintf("%s  (Nodes=%d  Edges=%d  Pos=%d  Neg=%d)",
                       title_map[grp_name], net$n_nodes, net$n_edges,
                       net$n_pos, net$n_neg),
        cex.main = 0.6, line = 0.3)
}

# 底部共享图例
par(mar = c(0, 0.5, 0, 0.5))
plot.new()

phyla_all <- names(sort(table(unlist(lapply(networks, function(n) V(n$graph)$phylum))),
                         decreasing = TRUE))
phylum_cols_all <- phylum_colors[phyla_all]
phylum_cols_all[is.na(phylum_cols_all)] <- "#CCCCCC"

# 水平排列：门方块 + 大小参考 + 边
n_items <- length(phyla_all)
total_n <- n_items + 5  # phyla + spacer + 3 sizes + pos + neg
x_starts <- seq(0.02, 0.98, length.out = total_n)

# 门方块
for (i in 1:n_items) {
  xc <- x_starts[i]
  rect(xc-0.008, 0.6, xc+0.008, 0.8, col = phylum_cols_all[i], border = "white")
  text(xc, 0.52, phyla_all[i], cex = 0.35, srt = 45, adj = c(1, 0.5))
}

# 节点大小参考
sz_start <- x_starts[n_items + 1]
text(sz_start, 0.85, "Degree:", cex = 0.4, font = 2, adj = c(0.5, 0))
for (j in 1:3) {
  sx <- sz_start + (j-1)*0.025
  
  text(sx, 0.4, names(size_bins)[j], cex = 0.3, adj = c(0.5, 1))
}

# 边
edge_start <- x_starts[n_items + 4]
segments(edge_start-0.01, 0.75, edge_start+0.015, 0.75, col = POS_COLOR, lwd = 1.5)
text(edge_start+0.025, 0.75, "Pos", cex = 0.35, adj = c(0, 0.5))
segments(edge_start-0.01, 0.5, edge_start+0.015, 0.5, col = NEG_COLOR, lwd = 1.5)
text(edge_start+0.025, 0.5, "Neg", cex = 0.35, adj = c(0, 0.5))

dev.off()
cat(sprintf("  -> %s/combined_4_networks.svg\n", outdir))

# ================================================================
# ========== 7. 网络拓扑对比表 ==========
# ================================================================
cat("\nStep 7: Network topology comparison...\n")

safe_get <- function(lst, key) {
  sapply(lst, function(x) {
    v <- x[[key]]
    if (is.null(v) || length(v) == 0) return(0)
    as.numeric(v)
  })
}

topology <- data.frame(
  Group = names(networks),
  Nodes = safe_get(networks, "n_nodes"),
  Edges = safe_get(networks, "n_edges"),
  Positive = safe_get(networks, "n_pos"),
  Negative = safe_get(networks, "n_neg"),
  PosNegRatio = round(safe_get(networks, "n_pos") / pmax(1, safe_get(networks, "n_neg")), 2),
  AvgDegree = round(safe_get(networks, "avg_degree"), 2),
  Density = round(safe_get(networks, "density"), 4),
  Components = safe_get(networks, "n_components"),
  LargestComp = safe_get(networks, "largest_component"),
  row.names = NULL
)

print(topology)
write.csv(topology, sprintf("%s/network_topology_comparison.csv", outdir), row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(topology, sprintf("%s/network_topology_comparison.xlsx", outdir))
}

# ================================================================
# ========== 8. 导出边列表和节点属性 ==========
# ================================================================
cat("\nStep 8: Exporting edge lists and node attributes...\n")

for (grp_name in names(networks)) {
  net <- networks[[grp_name]]
  g <- net$graph

  if (nrow(net$edge_list) > 0) {
    write.csv(net$edge_list, sprintf("%s/%s_edge_list.csv", outdir, grp_name), row.names = FALSE)
  } else {
    write.csv(data.frame(Note = "No edges passed filter"),
              sprintf("%s/%s_edge_list.csv", outdir, grp_name), row.names = FALSE)
  }

  bet <- betweenness(g, normalized = TRUE)
  node_df <- data.frame(
    Genus = V(g)$name,
    Phylum = V(g)$phylum,
    Mean_Abundance = round(V(g)$mean_abund, 1),
    Degree = degree(g),
    Betweenness = round(bet, 4),
    stringsAsFactors = FALSE
  )
  node_df <- node_df[order(node_df$Degree, decreasing = TRUE), ]
  write.csv(node_df, sprintf("%s/%s_node_attributes.csv", outdir, grp_name), row.names = FALSE)

  cat(sprintf("  %s: %d edges, %d nodes\n", grp_name, nrow(net$edge_list), nrow(node_df)))
}

cat(sprintf("\nAll outputs saved to: %s\n", file.path(getwd(), outdir)))
cat("============================================================\n")
cat("DONE!\n")
