# ================================================================
# 随机森林 — 合并 SS+NSS，4组（n=6/组）识别关键指示属
# ================================================================

library(randomForest)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
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
dir.create("output_r", showWarnings = FALSE)
outdir <- "output_r"

cat("============================================================\n")
cat("Random Forest (SS+NSS combined, n=6/group)\n")
cat("============================================================\n")

# ---------- 1. 预处理 ----------
cat("Step 1: Reading data...\n")
df <- read.delim("../data/microbial_OTU_table.xls", sep = "\t", header = TRUE, check.names = FALSE)
data_cols <- grep("^W[0-9]+_[0-9]+$", names(df), value = TRUE)
df$genus_clean <- gsub("^g__", "", df$genus)
genus_abund <- aggregate(df[, data_cols], by = list(genus = df$genus_clean), FUN = sum)
rownames(genus_abund) <- genus_abund$genus; genus_abund$genus <- NULL
is_uncl <- grepl("unclassified", rownames(genus_abund), ignore.case = TRUE)
genus_abund_clean <- genus_abund[!is_uncl, ]
overall_mean <- sort(rowMeans(genus_abund_clean), decreasing = TRUE)

N_FEAT <- 30
topN <- names(overall_mean)[1:min(N_FEAT, length(overall_mean))]
cat(sprintf("  Features: top %d genera\n", length(topN)))

# ---------- 2. 分组 ----------
groups <- list(
  CK  = c("W1_1","W1_2","W1_3","W8_1","W8_2","W8_3"),
  S   = c("W2_1","W2_2","W2_3","W9_1","W9_2","W9_3"),
  SC1 = c("W3_1","W3_2","W3_3","W10_1","W10_2","W10_3"),
  SC5 = c("W7_1","W7_2","W7_3","W14_1","W14_2","W14_3")
)

all_samples <- unlist(groups)
feat_mat <- t(genus_abund_clean[topN, all_samples, drop = FALSE])
y <- factor(rep(names(groups), each = 6), levels = names(groups))
cat(sprintf("  Samples=%d, Features=%d, Classes=%d\n", nrow(feat_mat), ncol(feat_mat), length(levels(y))))

# ---------- 3. 训练 ----------
cat("Step 3: Training RF...\n")
set.seed(42)
rf <- randomForest(x = feat_mat, y = y, ntree = 2000, mtry = 6, importance = TRUE)
cat(sprintf("  OOB error: %.1f%%\n", rf$err.rate[2000,1]*100))
print(rf$confusion)

imp <- importance(rf, type = 1)
imp_sorted <- sort(imp[,1], decreasing = TRUE)

# ---------- 4. 绘图 ----------
cat("Step 4: Plotting...\n")
MM <- 1/25.4

clean_names <- function(x) {
  x <- gsub("^norank_f__|^norank_o__|^norank_c__|^norank_p__|^norank_", "", x)
  x <- gsub("^f__|^o__|^c__|^p__", "", x)
  x <- gsub("__.*$", "", x); x
}

# 4a. Importance (positive only, orange→blue gradient)
imp_pos <- imp_sorted[imp_sorted > 0]                         # 只保留正值
n_pos <- length(imp_pos)
imp_df <- data.frame(
  Genus = factor(clean_names(names(imp_pos)), levels = clean_names(names(imp_pos))),
  Importance = as.numeric(imp_pos)
)
# 从上到下渐变：暗橙→暗蓝（无白过渡）
grad_cols <- colorRampPalette(c("#A63603", "#D95F02", "#2C6A9E", "#08519C"))(n_pos)
bar_cols <- rev(grad_cols)  # 最大值(顶部)橙色，底部蓝色

svglite(sprintf("%s/RF_importance.svg", outdir), width = 150*MM, height = max(80, n_pos*7)*MM)
par(mar = c(4, 12, 1, 2), mgp = c(2.8, 0.6, 0))
bp <- barplot(rev(imp_df$Importance), names.arg = rev(imp_df$Genus),
              horiz = TRUE, las = 1, col = bar_cols, border = NA,
              xlab = "Mean Decrease Accuracy (%)",
              ylab = "", xlim = c(0, 10),
              cex.names = 0.65, cex.axis = 0.65, cex.lab = 0.75)
# 标签在条形外侧
text(rev(imp_df$Importance) + max(imp_df$Importance)*0.04, bp,
     labels = sprintf("%.1f", rev(imp_df$Importance)), cex = 0.55, adj = c(0, 0.5),
     xpd = TRUE)
# 只保留左轴和下轴
axis(1, labels = NA, tck = -0.01)
axis(2, labels = NA, tck = -0.01)
box(bty = "L", col = "gray50")
dev.off(); cat(sprintf("  -> RF_importance.svg (%d positive genera)\n", n_pos))

# 4b. Confusion matrix
cm_pct <- sweep(rf$confusion[, 1:4], 1, rowSums(rf$confusion[, 1:4]), "/") * 100
svglite(sprintf("%s/RF_confusion.svg", outdir), width = 85*MM, height = 75*MM)
pheatmap(cm_pct, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.0f", number_color = "white",
         fontsize_number = 9, fontsize = 8, fontfamily = "sans",
         color = colorRampPalette(c("#F7FBFF", "#3182BD", "#08519C"))(50),
         main = "Confusion Matrix (%)", angle_col = 0,
         border_color = "white", cellwidth = 32, cellheight = 32)
dev.off(); cat("  -> RF_confusion.svg\n")

# 4c. Heatmap (top 20 by importance, all values)
top20g <- names(head(imp_sorted, 20))
mean_mat <- matrix(NA, nrow = 20, ncol = 4)
rownames(mean_mat) <- clean_names(top20g)
colnames(mean_mat) <- names(groups)
for (i in 1:4) {
  grp_samples <- groups[[i]]
  grp_data <- feat_mat[grp_samples, top20g, drop = FALSE]
  mean_mat[, i] <- colMeans(grp_data)
}
	svglite(sprintf("%s/RF_heatmap.svg", outdir), width = 140*MM, height = 110*MM)
	pheatmap(log10(mean_mat + 1), cluster_rows = TRUE, cluster_cols = FALSE,
	         scale = "row", fontsize = 7, fontfamily = "sans", fontsize_col = 9,
	         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
	         main = "Top 20 Biomarker Genera", angle_col = 0,
	         border_color = NA, cellwidth = 26, cellheight = 11,
	         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
	         legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
dev.off(); cat("  -> RF_heatmap.svg\n")

# 4d. Convergence
svglite(sprintf("%s/RF_convergence.svg", outdir), width = 100*MM, height = 70*MM)
par(mar = c(4, 4, 2, 1))
plot(rf$err.rate[,1], type = "l", col = "#4575B4", lwd = 1.5,
     xlab = "Trees", ylab = "OOB Error",
     main = sprintf("RF Convergence (final OOB=%.1f%%)", rf$err.rate[2000,1]*100),
     cex.lab = 0.7, cex.axis = 0.6, cex.main = 0.75)
grid(lty = "dotted", col = "gray80")
dev.off(); cat("  -> RF_convergence.svg\n")

# ---------- 5. Export ----------
cat("Step 5: Exporting...\n")
imp_out <- data.frame(Genus = names(imp_sorted),
                       MeanDecreaseAccuracy = round(as.numeric(imp_sorted), 2))
write.csv(imp_out, sprintf("%s/RF_importance.csv", outdir), row.names = FALSE)
write.csv(as.data.frame(rf$confusion), sprintf("%s/RF_confusion.csv", outdir))
cat(sprintf("  All saved to %s\n", file.path(getwd(), outdir)))
cat("DONE!\n")
