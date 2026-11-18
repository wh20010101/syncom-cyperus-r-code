library(readxl)
library(linkET)
library(tidyverse)
library(openxlsx)

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

# ============================================
# 通用 Mantel 分析函数
# ============================================
run_mantel <- function(micro_file, metab_file, phys_file, label) {
  cat(sprintf("\n========== %s ==========\n", label))

  env   <- read_excel(phys_file) %>% column_to_rownames("Sample")
  micro <- read_excel(micro_file) %>% column_to_rownames("Sample") %>% drop_na()
  metab <- read_excel(metab_file) %>% column_to_rownames("Sample")

  cat(sprintf("  Microbes: %d x %d\n", nrow(micro), ncol(micro)))
  cat(sprintf("  Metabolites: %d x %d\n", nrow(metab), ncol(metab)))
  cat(sprintf("  Physiology: %d x %d\n", nrow(env), ncol(env)))

  # 合并微生物+代谢物作为 spec
  spec <- cbind(micro, metab)
  n_micro <- ncol(micro)
  n_metab <- ncol(metab)

  spec_group <- list(
    Microbiome = 1:n_micro,
    Metabolome = (n_micro + 1):(n_micro + n_metab)
  )

  # Mantel test
  mantel <- mantel_test(spec, env,
                        spec_select = spec_group,
                        spec_dist   = "bray",
                        env_dist    = "euclidean",
                        seed        = 123) %>%
    mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

  # Spearman between physiology indicators
  corr <- correlate(env, method = "spearman")

  # Plot
  p <- qcorrplot(corr, type = "upper") +
    geom_square() +
    geom_text(
      data = {
        r_mat <- corr$r
        p_mat <- corr$p
        vnames <- colnames(r_mat)
        expand.grid(row = vnames, col = vnames) %>%
          filter(as.character(row) != as.character(col)) %>%
          rowwise() %>%
          mutate(
            p_value = p_mat[as.character(row), as.character(col)],
            sig = case_when(
              p_value < 0.001 ~ "***",
              p_value < 0.01  ~ "**",
              p_value < 0.05  ~ "*",
              TRUE            ~ ""
            )
          ) %>%
          ungroup() %>%
          filter(match(row, vnames) < match(col, vnames))
      },
      aes(x = col, y = row, label = sig),
      size = 4, color = "black", fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = 0.15) +
    scale_fill_gradient2(low  = "#2166AC", mid  = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman's r") +
    scale_size_manual(values = c(0.5, 2, 4)) +
    scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
    guides(size   = guide_legend(title = "Mantel's r",
                                 override.aes = list(colour = "grey35"), order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3), order = 1),
           fill   = guide_colorbar(title = "Spearman's r", order = 3)) +
    labs(title = label)

  ggsave(sprintf("Mantel_plot_%s.svg", label), plot = p,
         width = 14, height = 12, dpi = 300)

  # Export tables
  corr_r <- as.data.frame(corr$r)
  corr_r$env1 <- rownames(corr_r)
  corr_r <- corr_r %>% pivot_longer(-env1, names_to = "env2", values_to = "Spearman_r")

  corr_p <- as.data.frame(corr$p)
  corr_p$env1 <- rownames(corr_p)
  corr_p <- corr_p %>% pivot_longer(-env1, names_to = "env2", values_to = "p_value")

  spearman_table <- left_join(corr_r, corr_p, by = c("env1", "env2")) %>%
    filter(match(env1, colnames(corr$r)) < match(env2, colnames(corr$r))) %>%
    mutate(Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))

  mantel_table <- mantel %>% select(spec, env, r, p, rd, pd)

  write.xlsx(list(Spearman = spearman_table, Mantel = mantel_table),
             sprintf("Mantel_results_%s.xlsx", label))

  cat(sprintf("  -> Mantel_plot_%s.svg + Mantel_results_%s.xlsx\n\n", label, label))
  return(list(mantel = mantel, corr = corr))
}

# ============================================
# 灭菌土 (SS)
# ============================================
run_mantel("灭菌微生物.xlsx", "灭菌代谢物.xlsx", "灭菌指标.xlsx", "SS")

# ============================================
# 非灭菌土 (NSS)
# ============================================
run_mantel("非灭菌微生物.xlsx", "非灭菌代谢物.xlsx", "非灭菌指标.xlsx", "NSS")

cat("DONE!\n")
