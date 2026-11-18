###############################################################################
##  06_factorial_anova.R
##  2 × 7 双因素方差分析（38 个指标）
##  ---------------------------------------------------------------------------
##  设计：土壤背景 SS/NSS × 处理 CK,S,SC1-SC5 × 3 重复 = 42 个观测
##  全模型：y ~ soil * trt，误差 df = 28
##  计算：
##    - 三效应（soil、trt、soil:trt）的 F、P
##    - 偏 eta^2 = SS_effect / (SS_effect + SS_error)
##    - 偏 eta^2 的 95% CI：非中心 F 反解
##    - Levene（Brown-Forsythe，median）与 Shapiro-Wilk
##    - 每个土壤背景内 7 个处理的 Tukey HSD
##      HSD = q(0.95, k=7, df=28) * sqrt(MSE / n)，n = 3
##    - Piepho 插入-吸收法生成显著性字母
##  输入：优先 ../data/所有数据.xlsx；否则当前目录 所有数据.xlsx
##  输出：output/1_三效应_F_P_偏eta2_CI.xlsx
##        output/2_各处理均值与字母.xlsx
##  依赖：readxl, writexl；car 可选（平衡设计下 Type I/III 等价）
###############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
})

## ---- 自动定位到脚本所在目录 ----
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

## ---- 路径 ----
DATA_XLSX <- if (file.exists("../data/所有数据.xlsx")) {
  "../data/所有数据.xlsx"
} else {
  "所有数据.xlsx"
}
OUTDIR <- "output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DATA_XLSX)) {
  stop("Input data not found: ", DATA_XLSX,
       "\n请把 所有数据.xlsx 放到 ../data/ 或脚本所在目录。")
}

## ---- 读入与清洗 ----
d <- as.data.frame(read_excel(DATA_XLSX, sheet = "Sheet1", .name_repair = "minimal"))
stopifnot(nrow(d) == 42, ncol(d) == 39, "Sample" %in% names(d))

# 修正第 7–9 行原标 W2、实为 W3 的样本（与 Python 脚本一致）
# 仅在原表确实标为 W2 时触发，避免数据顺序变化时误改
if (grepl("^W2", as.character(d$Sample[7]))) {
  d[7:9, "Sample"] <- c("W3_1", "W3_2", "W3_3")
}

d$W <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", as.character(d$Sample)))
TRT <- c("CK", "S", "SC1", "SC2", "SC3", "SC4", "SC5")
d$trt <- factor(TRT[(d$W - 1) %% 7 + 1], levels = TRT)
d$soil <- factor(ifelse(d$W <= 7, "SS", "NSS"), levels = c("SS", "NSS"))

# 列名修正（与论文表述一致）
rename_map <- c(PH = "土壤pH", 中洗 = "中性洗涤纤维", 酸洗 = "酸性洗涤纤维")
for (nm in names(rename_map)) {
  if (nm %in% names(d)) names(d)[names(d) == nm] <- rename_map[[nm]]
}

traits <- setdiff(names(d), c("Sample", "W", "trt", "soil"))
stopifnot(length(traits) == 38, sum(is.na(d[traits])) == 0)
for (tr in traits) d[[tr]] <- as.numeric(d[[tr]])
stopifnot(all(table(d$soil, d$trt) == 3))

cat("Data:", nrow(d), "observations x", length(traits), "traits\n")
cat("Design: 2 soil backgrounds x 7 treatments x 3 replicates; error df = 28\n")

## ---- 偏 eta^2 的 95% CI ----
eta2_ci <- function(F, df1, df2, alpha = 0.05) {
  sf_fun <- function(lambda) {
    suppressWarnings(pf(F, df1, df2, ncp = lambda, lower.tail = FALSE))
  }
  solve_lambda <- function(target) {
    f <- function(l) sf_fun(l) - target
    if (!is.finite(f(1e-9))) return(NA_real_)
    if (f(1e-9) >= 0) return(0)
    if (f(1e9) <= 0) return(Inf)
    tryCatch(
      uniroot(f, c(1e-9, 1e9), tol = 1e-10)$root,
      error = function(e) NA_real_
    )
  }
  lo <- solve_lambda(alpha / 2)
  hi <- solve_lambda(1 - alpha / 2)
  conv <- function(l) {
    if (is.na(l)) return(NA_real_)
    if (is.infinite(l)) return(1)
    l / (l + df2)
  }
  c(lower = conv(lo), upper = conv(hi))
}

## ---- Levene（Brown-Forsythe，median）----
levene_p <- function(y, group) {
  med <- ave(y, group, FUN = median)
  z <- abs(y - med)
  fitz <- aov(z ~ group)
  summary(fitz)[[1]][["Pr(>F)"]][1]
}

## ---- Piepho 插入-吸收 CLD ----
cld_piepho <- function(nms, vals, sig) {
  ord <- order(-vals)
  cols <- list()
  for (i in ord) {
    placed <- FALSE
    for (ci in seq_along(cols)) {
      if (all(!sig[i, cols[[ci]]])) {
        cols[[ci]] <- c(cols[[ci]], i)
        placed <- TRUE
        break
      }
    }
    if (!placed) cols[[length(cols) + 1L]] <- i
  }

  for (i in ord) {
    for (j in ord) {
      if (i != j && !sig[i, j]) {
        already <- any(vapply(cols, function(cc) i %in% cc && j %in% cc, logical(1)))
        if (!already) {
          for (ci in seq_along(cols)) {
            if (i %in% cols[[ci]] && all(!sig[cols[[ci]], j])) {
              cols[[ci]] <- c(cols[[ci]], j)
              break
            }
          }
        }
      }
    }
  }

  cols <- cols[order(vapply(cols, function(cc) min(match(cc, ord)), numeric(1)))]
  out <- setNames(rep("", length(nms)), nms)
  for (ci in seq_along(cols)) {
    L <- letters[ci]
    for (i in cols[[ci]]) out[nms[i]] <- paste0(out[nms[i]], L)
  }
  vapply(out, function(x) paste(sort(strsplit(x, "")[[1]]), collapse = ""), character(1))
}

## ---- Tukey HSD（全模型 MSE，df = 28）----
tukey_full <- function(means, mse, dfe, n, k) {
  q <- qtukey(0.95, nmeans = k, df = dfe)
  hsd <- q * sqrt(mse / n)
  vals <- as.numeric(means)
  nms <- names(means)
  sig <- outer(seq_along(vals), seq_along(vals), function(i, j) abs(vals[i] - vals[j]) > hsd)
  list(hsd = hsd, letters = cld_piepho(nms, vals, sig), sig = sig)
}

## ---- 逐指标计算 ----
rows <- list()
letters_rows <- list()
k <- 0L
kl <- 0L

for (tr in traits) {
  d$y <- d[[tr]]
  fit <- aov(y ~ soil * trt, data = d)

  # 平衡设计下 Type I/II/III 等价；若 car 可用，用 type III 与 Methods 对齐
  if (requireNamespace("car", quietly = TRUE)) {
    tab <- car::Anova(fit, type = 3)
  } else {
    tab <- anova(fit)
  }

  ss_e <- tab["Residuals", "Sum Sq"]
  df_e <- tab["Residuals", "Df"]
  mse <- ss_e / df_e
  ss_cell <- sum(tab[c("soil", "trt", "soil:trt"), "Sum Sq"])
  r2 <- ss_cell / (ss_cell + ss_e)

  p_lev <- levene_p(d$y, interaction(d$soil, d$trt, drop = TRUE))
  p_sha <- shapiro.test(resid(fit))$p.value
  p_inter <- tab["soil:trt", "Pr(>F)"]

  get_effect <- function(term) {
    ss <- tab[term, "Sum Sq"]
    df1 <- tab[term, "Df"]
    Fv <- tab[term, "F value"]
    Pv <- tab[term, "Pr(>F)"]
    eta <- ss / (ss + ss_e)
    ci <- if (is.finite(Fv) && is.finite(df1) && df1 > 0) {
      eta2_ci(Fv, df1, df_e)
    } else {
      c(lower = NA_real_, upper = NA_real_)
    }
    list(SS = ss, df = df1, F = Fv, P = Pv, eta = eta, lo = ci["lower"], hi = ci["upper"])
  }

  a_soil <- get_effect("soil")
  a_trt  <- get_effect("trt")
  a_int  <- get_effect("soil:trt")

  fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  fmt_ci <- function(lo, hi) sprintf("%.3f-%.3f", lo, hi)

  k <- k + 1L
  rows[[k]] <- data.frame(
    指标 = tr,
    土壤_F = round(a_soil$F, 2), 土壤_P = fmt_p(a_soil$P),
    土壤_偏eta2 = round(a_soil$eta, 3), 土壤_eta2_95CI = fmt_ci(a_soil$lo, a_soil$hi),
    处理_F = round(a_trt$F, 2), 处理_P = fmt_p(a_trt$P),
    处理_偏eta2 = round(a_trt$eta, 3), 处理_eta2_95CI = fmt_ci(a_trt$lo, a_trt$hi),
    交互_F = round(a_int$F, 2), 交互_P = fmt_p(a_int$P),
    交互_偏eta2 = round(a_int$eta, 3), 交互_eta2_95CI = fmt_ci(a_int$lo, a_int$hi),
    MSE = round(mse, 4), 误差df = df_e, 模型R2 = round(r2, 3),
    Levene_P = round(p_lev, 3), 残差Shapiro_P = round(p_sha, 3),
    stringsAsFactors = FALSE
  )

  # 每个土壤背景内 7 个处理的 Tukey；全模型 MSE、df=28、n=3
  for (s in levels(d$soil)) {
    sub <- d[d$soil == s, , drop = FALSE]
    mu <- tapply(sub$y, sub$trt, mean)[TRT]
    sdv <- tapply(sub$y, sub$trt, sd)[TRT]
    nn <- tapply(sub$y, sub$trt, length)[TRT]
    tuk <- tukey_full(mu, mse, df_e, n = 3, k = 7)

    for (tt in TRT) {
      kl <- kl + 1L
      letters_rows[[kl]] <- data.frame(
        指标 = tr, 土壤 = s, 处理 = tt,
        均值 = round(as.numeric(mu[tt]), 3),
        SD = round(as.numeric(sdv[tt]), 3),
        n = as.integer(nn[tt]),
        字母 = tuk$letters[tt],
        HSD = round(tuk$hsd, 3),
        交互_P = fmt_p(p_inter),
        stringsAsFactors = FALSE
      )
    }
  }
}

df_res <- do.call(rbind, rows)
df_let <- do.call(rbind, letters_rows)

f1 <- file.path(OUTDIR, "1_三效应_F_P_偏eta2_CI.xlsx")
f2 <- file.path(OUTDIR, "2_各处理均值与字母.xlsx")
write_xlsx(list("ANOVA" = df_res), f1)
write_xlsx(list("means_letters" = df_let), f2)

cat("\nOutput:\n  ", normalizePath(f1), "\n  ", normalizePath(f2), "\n", sep = "")
p_int_num <- ifelse(df_res$交互_P == "<0.001", 0, suppressWarnings(as.numeric(df_res$交互_P)))
cat("\nInteraction significant (P<0.05): ",
    sum(p_int_num < 0.05, na.rm = TRUE),
    " / ", nrow(df_res), "\n", sep = "")

cat("\nSession info:\n")
print(sessionInfo())
