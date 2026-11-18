# ================================================================
# 模块 2：SS / NSS 代谢物—生理指标 Spearman 相关性网络分析
#
# 【输入】
#   metabolite_SS.csv、metabolite_NSS.csv
#       行 = 代谢物，列 = 注释信息 + 样本丰度（W1_1 … W14_6）
#   ../01_mantel/SS_physiology.xlsx、NSS_physiology.xlsx
#       行 = 样本（列名 Sample），列 = 各项生理指标
#
# 【输出】写入同目录 output/
#   {SS,NSS}_metabolite_network.svg   双圈网络图（外圈代谢物、内圈生理指标）
#   {SS,NSS}_metabolite_edges.csv     显著相关的边（节点对、相关系数、P 值）
#   {SS,NSS}_metabolite_nodes.csv     节点属性（类型、度、中介中心性）
#
# 【筛选标准】|Spearman rho| >= 0.7 且 Raw P < 0.05（未做 FDR 多重检验校正）
# 【说明】本分析不使用热图，仅输出网络图与数据表
# ================================================================

# ---- 依赖包 ----
library(readxl)    # 读取 xlsx 生理指标表
library(igraph)    # 构建与计算网络
library(svglite)   # 输出矢量图 SVG

# ================================================================
# 路径适配：自动把工作目录切换到脚本所在目录
#   - 用 Rscript 运行时，从 --file 参数解析路径
#   - 在 RStudio 中 Source 时，从当前编辑文档取路径
# 这样读者把仓库克隆到任意位置都能直接运行，无需改代码
# ================================================================
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

# ---- 全局参数 ----
dir.create("output", showWarnings = FALSE)  # 输出目录，已存在则不报错
MM <- 1/25.4                                 # 毫米转英寸，用于控制出图物理尺寸

# ---- 配色 ----
POS_COLOR  <- "#E64B35"  # 正相关边：红
NEG_COLOR  <- "#4575B4"  # 负相关边：蓝
METAB_COL  <- "#1B9E77"  # 代谢物节点：绿
PHYS_COL   <- "#D95F02"  # 生理指标节点：橙

# ================================================================
# 核心函数：对单一处理（SS 或 NSS）完成相关性计算与网络绘制
# 参数：
#   metab_file  代谢物表路径
#   phys_file   生理指标表路径
#   soil_label  处理标签（"SS" / "NSS"），用于输出文件命名
#   sample_map  分组 → 样本名向量的映射（决定用哪些样本参与计算）
# ================================================================
run_analysis <- function(metab_file, phys_file, soil_label, sample_map) {
  cat(sprintf("\n========== %s ==========\n", soil_label))

  # ---- 1. 读取代谢物表 ----
  # check.names = FALSE 保留列名原样（列名含空格、斜杠等特殊字符）
  metab_raw <- read.csv(metab_file, header = TRUE, check.names = FALSE)
  # 以代谢物名作为行名（原表第一列是 ID，第二列 Metabolite 才是名称）
  rownames(metab_raw) <- metab_raw$Metabolite
  # 只取 sample_map 中列出的样本列，丢弃注释列和其他样本
  metab_mat <- metab_raw[, unlist(sample_map), drop = FALSE]
  # 统一转为数值型，防止被读成字符
  metab_mat <- as.data.frame(lapply(metab_mat, as.numeric))
  rownames(metab_mat) <- rownames(metab_raw)
  cat(sprintf("  Metabolites: %d\n", nrow(metab_mat)))

  # ---- 2. 读取生理指标表 ----
  # 第一列 Sample 作为行名，其余列均为指标
  phys_raw <- read_excel(phys_file)
  phys_mat <- as.data.frame(phys_raw)
  rownames(phys_mat) <- phys_mat$Sample
  phys_mat$Sample <- NULL

  n_m <- nrow(metab_mat)   # 代谢物数
  n_p <- ncol(phys_mat)    # 生理指标数

  # ---- 3. 逐对计算 Spearman 相关 ----
  # exact = FALSE：样本量小且存在并列值时使用近似 P 值，避免计算失败
  r_mat <- matrix(NA, n_m, n_p)
  p_mat <- matrix(NA, n_m, n_p)
  for (i in 1:n_m) for (j in 1:n_p) {
    test <- cor.test(as.numeric(metab_mat[i, ]), phys_mat[[j]],
                     method = "spearman", exact = FALSE)
    r_mat[i, j] <- test$estimate   # 相关系数 rho
    p_mat[i, j] <- test$p.value    # 显著性 P 值
  }
  rownames(r_mat) <- rownames(metab_mat); colnames(r_mat) <- colnames(phys_mat)
  rownames(p_mat) <- rownames(metab_mat); colnames(p_mat) <- colnames(phys_mat)

  # ---- 4. P 值 ----
  # 使用原始 P（未做 FDR 校正）；如需 FDR，可对 p_mat 做 BH 校正。
  p_raw_mat <- p_mat

  n_sig <- sum(p_raw_mat < 0.05, na.rm = TRUE)
  cat(sprintf("  Tests: %d, Raw P<0.05: %d\n", n_m * n_p, n_sig))

  # ---- 5. 构建边表 ----
  # 同时满足两个条件才保留为一条边：
  #   ① P < 0.05  ② |rho| >= 0.7（强相关）
  edge_list <- data.frame()
  for (i in 1:nrow(r_mat)) for (j in 1:ncol(r_mat)) {
    if (!is.na(p_raw_mat[i, j]) && p_raw_mat[i, j] < 0.05 && abs(r_mat[i, j]) >= 0.7) {
      edge_list <- rbind(edge_list, data.frame(
        Node1 = rownames(r_mat)[i], Node2 = colnames(r_mat)[j],
        Coefficient = round(r_mat[i, j], 3), P_value = signif(p_raw_mat[i, j], 3),
        Direction = ifelse(r_mat[i, j] > 0, "positive", "negative"),
        stringsAsFactors = FALSE))
    }
  }

  # 该处理下没有强相关时，跳过绘图直接返回
  if (nrow(edge_list) == 0) {
    cat(sprintf("  %s: no significant edges\n", soil_label))
    return(list(n_sig = 0))
  }

  cat(sprintf("  Edges: %d (pos=%d, neg=%d)\n", nrow(edge_list),
              sum(edge_list$Direction == "positive"),
              sum(edge_list$Direction == "negative")))

  # ---- 6. 建立 igraph 网络对象 ----
  # 先建空图再加点、加边，避免节点名重复
  all_nodes <- unique(c(edge_list$Node1, edge_list$Node2))
  g <- make_empty_graph(n = length(all_nodes), directed = FALSE)
  V(g)$name  <- all_nodes
  # 节点类型：出现在 r_mat 行名中即代谢物，否则为生理指标
  V(g)$type  <- ifelse(all_nodes %in% rownames(r_mat), "Metabolite", "Physiology")
  V(g)$color <- ifelse(V(g)$type == "Metabolite", METAB_COL, PHYS_COL)

  for (k in 1:nrow(edge_list)) {
    g <- add_edges(g, c(edge_list$Node1[k], edge_list$Node2[k]),
                   weight      = abs(edge_list$Coefficient[k]),  # 边宽用相关强度
                   coefficient = edge_list$Coefficient[k],
                   sign        = edge_list$Direction[k])
  }

  # ---- 7. 双圈布局坐标 ----
  # 外圈半径 1.0 放代谢物，内圈半径 0.65 放生理指标，按等角度均匀分布
  metab_names <- sort(unique(edge_list$Node1))
  phys_names  <- sort(unique(edge_list$Node2))
  n_m <- length(metab_names); n_p <- length(phys_names)

  r_outer <- 1.0   # 外圈 = 代谢物
  r_inner <- 0.65  # 内圈 = 生理指标

  angles_m <- seq(0, 2 * pi, length.out = n_m + 1)[1:n_m]
  coords_m <- data.frame(name = metab_names, x = r_outer * cos(angles_m),
                         y = r_outer * sin(angles_m), stringsAsFactors = FALSE)
  angles_p <- seq(0, 2 * pi, length.out = n_p + 1)[1:n_p]
  coords_p <- data.frame(name = phys_names, x = r_inner * cos(angles_p),
                         y = r_inner * sin(angles_p), stringsAsFactors = FALSE)
  all_coords <- rbind(coords_m, coords_p)
  rownames(all_coords) <- c(coords_m$name, coords_p$name)

  # 节点大小按度分组：度 >= 85 分位数 → 大；>= 中位数 → 中；其余 → 小
  degs <- degree(g)
  deg_m <- degs[coords_m$name]; deg_p <- degs[coords_p$name]
  deg_all <- c(deg_m, deg_p)
  dbrk <- quantile(deg_all, probs = c(0, 0.5, 0.85, 1))
  sz_fn <- function(d) ifelse(d >= dbrk[3], 4.0, ifelse(d >= dbrk[2], 2.4, 1.2))
  coords_m$cex <- sz_fn(deg_m); coords_p$cex <- sz_fn(deg_p)

  n_pos <- sum(edge_list$Direction == "positive")
  n_neg <- sum(edge_list$Direction == "negative")

  # ---- 8. 绘图 ----
  # 100 mm × 100 mm 正方形画布，字体 Times New Roman（论文常用）
  svglite(sprintf("output/%s_metabolite_network.svg", soil_label),
          width = 100 * MM, height = 100 * MM,
          pointsize = 8, system_fonts = list(serif = "Times New Roman"))

  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white", family = "serif")
  plot(0, 0, type = "n", xlim = c(-1.35, 1.35), ylim = c(-1.25, 1.25),
       xlab = "", ylab = "", axes = FALSE, asp = 1)

  # 8.1 画边：颜色区分正负相关，透明度与线宽随相关系数绝对值变化
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

  # 8.2 画节点：实心圆点，灰色描边
  # 外圈代谢物
  points(coords_m$x, coords_m$y, pch = 21, bg = adjustcolor(METAB_COL, 0.85),
         col = "grey40", cex = coords_m$cex)
  # 内圈生理指标
  points(coords_p$x, coords_p$y, pch = 21, bg = adjustcolor(PHYS_COL, 0.85),
         col = "grey40", cex = coords_p$cex)

  # 8.3 画节点标签
  # 代谢物标签（外圈）
  for (i in 1:n_m) {
    ang <- angles_m[i]; tx <- r_outer * cos(ang); ty <- r_outer * sin(ang)
    text(tx, ty, labels = coords_m$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }
  # 指标标签（内圈）
  for (i in 1:n_p) {
    ang <- angles_p[i]; tx <- r_inner * cos(ang); ty <- r_inner * sin(ang)
    text(tx, ty, labels = coords_p$name[i], srt = 0, cex = 0.4, col = "black", font = 2)
  }

  # 8.4 图例：节点大小分级 + 节点类型 + 边的正负，并标注各类数量
  n_sm <- sum(c(deg_m, deg_p) < dbrk[2])
  n_md <- sum(c(deg_m, deg_p) >= dbrk[2] & c(deg_m, deg_p) < dbrk[3])
  n_lg <- sum(c(deg_m, deg_p) >= dbrk[3])
  legend("bottomleft",
         legend = c(
           sprintf("Large (deg≥%.0f, n=%d)", dbrk[3], n_lg),
           sprintf("Medium (deg %.0f–%.0f, n=%d)", dbrk[2], dbrk[3]-1, n_md),
           sprintf("Small (deg<%.0f, n=%d)", dbrk[2], n_sm),
           "Metabolite", "Physiology",
           sprintf("Pos (n=%d)", n_pos), sprintf("Neg (n=%d)", n_neg)),
         pch = c(21, 21, 21, 21, 21, NA, NA),
         pt.bg = c("grey60", "grey60", "grey60", METAB_COL, PHYS_COL, NA, NA),
         col = c("grey40", "grey40", "grey40", "grey40", "grey40", POS_COLOR, NEG_COLOR),
         pt.cex = c(4.0, 2.4, 1.2, 1.0, 1.0, NA, NA),
         lty = c(NA, NA, NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, NA, NA, 1.2, 1.2),
         cex = 0.55, bty = "n")

  dev.off()
  cat(sprintf("  -> %s_metabolite_network.svg\n", soil_label))

  # ---- 9. 导出数据表 ----
  # 边表按 P 值升序排列，最显著的在前
  edge_list <- edge_list[order(edge_list$P_value), ]
  write.csv(edge_list, sprintf("output/%s_metabolite_edges.csv", soil_label), row.names = FALSE)

  # 节点表：类型、度（连接数）、中介中心性（节点在网络中的桥接作用）
  node_df <- data.frame(Name = V(g)$name, Type = V(g)$type, Degree = degree(g),
                         Betweenness = round(betweenness(g, normalized = TRUE), 4))
  node_df <- node_df[order(node_df$Degree, decreasing = TRUE), ]
  write.csv(node_df, sprintf("output/%s_metabolite_nodes.csv", soil_label), row.names = FALSE)

  # 返回结果供主流程汇总
  list(n_sig = n_sig, edge_list = edge_list, graph = g)
}

# ================================================================
# 主流程：分别对 SS（灭菌）和 NSS（非灭菌）两组执行分析
# ================================================================
cat("============================================================\n")
cat("SS & NSS Metabolite-Physiology Correlation\n")
cat("============================================================\n")

# 样本映射：CK / S / SC1 / SC5 四个处理，每个处理取前 3 个重复
# 注意：原始表每个处理有 6 个重复（W*_1 … W*_6）；本分析使用前 3 个（详见论文 Methods）
ss_map  <- list(CK  = c("W1_1", "W1_2", "W1_3"),  S   = c("W2_1", "W2_2", "W2_3"),
                SC1 = c("W3_1", "W3_2", "W3_3"),  SC5 = c("W7_1", "W7_2", "W7_3"))
nss_map <- list(CK  = c("W8_1", "W8_2", "W8_3"),  S   = c("W9_1", "W9_2", "W9_3"),
                SC1 = c("W10_1", "W10_2", "W10_3"), SC5 = c("W14_1", "W14_2", "W14_3"))

# 生理指标表来自模块 1 的输出目录，使用相对路径跨目录读取
phys_ss  <- "../01_mantel/SS_physiology.xlsx"
phys_nss <- "../01_mantel/NSS_physiology.xlsx"

ss_res  <- run_analysis("metabolite_SS.csv",  phys_ss,  "SS",  ss_map)
nss_res <- run_analysis("metabolite_NSS.csv", phys_nss, "NSS", nss_map)

cat(sprintf("\nFinal: SS=%d sig, NSS=%d sig\n", ss_res$n_sig, nss_res$n_sig))
cat("DONE!\n")
