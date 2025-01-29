###############################################################
# 
#    ライブラリ・データの読み込みと作業ディレクトリ設定
# 
###############################################################

# 必要なライブラリの読み込み
library(drc)  # 用量反応曲線用
library(tidyverse)
library(rstatix)  # 統計検定用
library(ggpubr)   # プロット用
library(dplyr)
library(stats)

# 作業ディレクトリを設定
setwd("/media/sato-lab/ボリューム/developments/ERG/ERG_R")

# データの読み込み
data <- read.csv("./ERG.csv", header = TRUE, fileEncoding = "UTF-8")

###############################################################
# 
#    対数光強度の計算と、logisticのフィッティング
# 
###############################################################

# 光の対数強度の計算
data$log_intensity <- log10(data$LightIntensity)

# 因子型への変換
data$Infection_Status <- as.factor(data$Infection_Status)
data$Polarized <- as.factor(data$Polarized)
data$Year <- as.factor(data$Year)

# contrastsの設定を変更
contrasts(data$Infection_Status) <- NULL
contrasts(data$Polarized) <- NULL
contrasts(data$Year) <- NULL

# 各IDごとにフィッティングを行う関数
fit_logistic_by_id <- function(data) {
  results <- data.frame()

  for (id in unique(data$ID)) {
    # 各IDのInfection_Status取得
    inf_status <- unique(data$Infection_Status[data$ID == id])

    for (pol in levels(data$Polarized)) {
      subset_data <- data[data$ID == id & data$Polarized == pol, ]

      if (nrow(subset_data) < 4) next

      tryCatch({
        model <- drm(ERG ~ log_intensity,
                     data = subset_data,
                     fct = L.4(),
                     control = drmc(maxIt = 2000,
                                    method = "BFGS",
                                    relTol = 1e-7))

        params <- coef(model)

        results <- rbind(results, data.frame(
          ID = id,
          Infection_Status = inf_status,
          Polarized = pol,
          Lower_Asymptote = params[1],
          Upper_Asymptote = params[2],
          Midpoint = params[3],
          Slope = params[4],
          R_squared = 1 - sum((residuals(model))^2) /
            sum((subset_data$ERG - mean(subset_data$ERG))^2)
        ))
      }, error = function(e) {
        warning(paste("ID", id, "Polarized", pol, "でフィッティングに失敗:", e$message))
      })
    }
  }
  return(results)
}

# フィッティング実行
fitting_results <- fit_logistic_by_id(data)

# 結果の表示
print(fitting_results)


###############################################################
# 
#    フィッティング結果の図示
# 
###############################################################


# IDごとのフィッティングプロット関数
plot_fits_by_id <- function(data, fitting_results) {
  for (id in unique(fitting_results$ID)) {
    id_data <- data[data$ID == id,]
    id_fits <- fitting_results[fitting_results$ID == id,]
    
    p <- ggplot(id_data, aes(x = log_intensity, y = ERG, color = Polarized)) +
      geom_point() +
      geom_smooth(method = "drm", 
                  method.args = list(fct = L.4()),
                  se = FALSE) +
      facet_wrap(~Polarized) +
      labs(title = paste("ID:", id),
           subtitle = paste("Infection Status:", unique(id_data$Infection_Status)),
           x = "Log Intensity",
           y = "ERG Response") +
      theme_minimal()
    
    ggsave(paste0("fit_ID_", id, ".png"), p, width = 10, height = 6, dpi = 300)
  }
}

# 実行
plot_fits_by_id(data, fitting_results)

# パラメータごとの箱ひげ図を作成
plot_params <- function(results) {
  results_long <- tidyr::pivot_longer(results, 
                                      cols = c(Lower_Asymptote, Upper_Asymptote, Midpoint, Slope),
                                      names_to = "Parameter",
                                      values_to = "Value")
  
  ggplot(results_long, aes(x = Polarized, y = Value, fill = Infection_Status)) +
    geom_boxplot() +
    facet_wrap(~Parameter, scales = "free_y") +
    theme_minimal() +
    labs(x = "Polarization", y = "Parameter Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# プロット作成とpng保存
p <- plot_params(fitting_results)
ggsave("ERG_parameters_boxplot.png", p, width = 10, height = 8, dpi = 300)

# CSVにも保存
write.csv(fitting_results, "ERG_fitting_results.csv", row.names = FALSE)

analyze_and_plot_midpoint <- function(results) {
  midpoint_data <- results %>% 
    dplyr::select("ID", "Infection_Status", "Polarized", "Midpoint")
  
  # Mann-Whitney U検定
  stat.test <- data.frame()
  for(pol in unique(midpoint_data$Polarized)) {
    subset <- midpoint_data[midpoint_data$Polarized == pol,]
    w.test <- wilcox.test(Midpoint ~ Infection_Status, data = subset)
    stat.test <- rbind(stat.test, 
                       data.frame(Polarized = pol,
                                  p = w.test$p.value))
  }
  stat.test$p.adj <- p.adjust(stat.test$p, method = "bonferroni")
  
  # Kruskal-Wallis検定
  kw.test <- kruskal.test(Midpoint ~ interaction(Infection_Status, Polarized), 
                          data = midpoint_data)
  
  # プロット作成
  p <- ggplot(midpoint_data, aes(x = Polarized, y = Midpoint, fill = Infection_Status)) +
    geom_boxplot() +
    labs(title = "Midpoint comparison across conditions",
         subtitle = paste("Kruskal-Wallis p-value:", round(kw.test$p.value, 4))) +
    theme_minimal()
  
  # 結果表示用の関数
  print("Kruskal-Wallis test results:")
  print(kw.test)
  
  print("\nMann-Whitney U test results (by polarization):")
  print(stat.test)
  
  list(plot = p, 
       kruskal = kw.test,
       mann_whitney = stat.test)
}

# 分析実行
results <- analyze_and_plot_midpoint(fitting_results)

# プロット保存
ggsave("ERG_midpoint_comparison.png", results$plot, width = 10, height = 6, dpi = 300)


###############################################################
# 
#    ART ANOVA で、Midpoint をランクに変換（ノンパラメトリックに処理）し、
#    交互作用を考慮して分析
# 
#    帰無仮説
#       主効果: 感染の有無は変曲点に影響しない
#       主効果: 偏光の種類は変曲点に影響しない
#       交互作用: 感染状態と偏光種の組み合わせによる変曲点への影響はない
# 
###############################################################


# 必要なライブラリ
library(ARTool)
library(ggplot2)

# データの準備
fitting_results$Polarized <- factor(fitting_results$Polarized)
fitting_results$Infection_Status <- factor(fitting_results$Infection_Status)

# ARTによる分析
m <- art(Midpoint ~ Infection_Status * Polarized + (1|ID), data = fitting_results)
art_result <- anova(m)

# 順位データを抽出してデータフレームに追加
fitting_results$art_ranks <- m$aligned.ranks

# 結果をテキストファイルに保存
sink("art_anova_results.txt")
print(art_result)
sink()


###############################################################
# 
#    ART ANOVA の分析に対応した可視化
# 
###############################################################

library(ggplot2)
library(reshape2)

# データの前処理：行名をインデックスに変換
fitting_results <- data.frame(lapply(fitting_results, as.vector))
rownames(fitting_results) <- NULL

# オリジナルのMidpointを使用したプロット
p0 <- ggplot(fitting_results, aes(x=Infection_Status, y=Midpoint, fill=Polarized)) +
  geom_violin() +
  geom_point(position=position_jitterdodge()) +
  theme_minimal() +
  labs(title="Raw Midpoint Values by Condition",
       y="Midpoint")
ggsave("midpoint_violin_plot.png", p0)

# Midpointの単純なランク付け
fitting_results$simple_rank <- rank(fitting_results$Midpoint)

# Violin plot
p1 <- ggplot(fitting_results, aes(x=Infection_Status, y=simple_rank, fill=Polarized)) +
  geom_violin() +
  geom_point(position=position_jitterdodge()) +
  theme_minimal() +
  labs(title="Aligned Ranks by Infection Status",
       y="Aligned Rank Response")
ggsave("ranked_midpoint_violin_plot.png", p1)

# Interaction plot
p2 <- ggplot(fitting_results, aes(x=Infection_Status, y=simple_rank, color=Polarized, group=Polarized)) +
  stat_summary(fun=mean, geom="point", size=5) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  theme_minimal() +
  theme(text = element_text(size = 24)) +
  labs(title="Interaction Effect on Ranks",
       y="Aligned Rank Response")
ggsave("ranked_midpoint_interaction_plot.png", p2)

# Heatmap
rank_matrix <- dcast(fitting_results, Infection_Status ~ Polarized, 
                     value.var="simple_rank", fun.aggregate=mean)
p3 <- ggplot(melt(rank_matrix), aes(x=Infection_Status, y=variable, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title="Rank Distribution Heatmap", 
       y="Polarized",
       fill="Mean Aligned Rank")
ggsave("ranked_midpoint_heatmap.png", p3)


# art_ranks.Infection_Status
# 感染状態（Infection_Status）の主効果を評価するためのアラインド（調整済み）ランク
# 他の要因（PolarizedとID）の影響を除去して、感染状態の効果のみを反映
# 
# art_ranks.Polarized
# 偏光（Polarized）の主効果を評価するためのアラインド・ランク
# 他の要因（Infection_StatusとID）の影響を除去して、偏光の効果のみを反映
# 
# art_ranks.Infection_Status:Polarized
# 交互作用効果を評価するためのアラインド・ランク
# 主効果を除去して、感染状態と偏光の組み合わせによる特異的な効果を反映


# Violin plot
p1 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Infection_Status, fill=Polarized)) +
  geom_violin() +
  geom_point(position=position_jitterdodge()) +
  theme_minimal() +
  labs(title="Aligned Ranks by Infection Status",
       y="Aligned Rank Response")
ggsave("infection_Status_violin_plot.png", p1)

# Interaction plot
p2 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Infection_Status, color=Polarized, group=Polarized)) +
  stat_summary(fun=mean, geom="point", size=5) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  theme_minimal() +
  theme(text = element_text(size = 24)) +
  labs(title="Interaction Effect on Ranks",
       y="Aligned Rank Response")
ggsave("infection_Status_interaction_plot.png", p2)

# Heatmap
rank_matrix <- dcast(fitting_results, Infection_Status ~ Polarized, 
                     value.var="art_ranks.Infection_Status", fun.aggregate=mean)
p3 <- ggplot(melt(rank_matrix), aes(x=Infection_Status, y=variable, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title="Rank Distribution Heatmap", 
       y="Polarized",
       fill="Mean Aligned Rank")
ggsave("infection_Status_heatmap.png", p3)



# Violin plot
p1 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Polarized, fill=Polarized)) +
  geom_violin() +
  geom_point(position=position_jitterdodge()) +
  theme_minimal() +
  labs(title="Aligned Ranks by Infection Status",
       y="Aligned Rank Response")
ggsave("polarized_violin_plot.png", p1)

# Interaction plot
p2 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Polarized, color=Polarized, group=Polarized)) +
  stat_summary(fun=mean, geom="point", size=3) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  theme_minimal() +
  labs(title="Interaction Effect on Ranks",
       y="Aligned Rank Response")
ggsave("polarized_interaction_plot.png", p2)

# Heatmap
rank_matrix <- dcast(fitting_results, Infection_Status ~ Polarized, 
                     value.var="art_ranks.Polarized", fun.aggregate=mean)
p3 <- ggplot(melt(rank_matrix), aes(x=Infection_Status, y=variable, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title="Rank Distribution Heatmap", 
       y="Polarized",
       fill="Mean Aligned Rank")
ggsave("polarized_heatmap.png", p3)



# Violin plot
p1 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Infection_Status.Polarized, fill=Polarized)) +
  geom_violin() +
  geom_point(position=position_jitterdodge()) +
  theme_minimal() +
  labs(title="Aligned Ranks by Infection Status",
       y="Aligned Rank Response")
ggsave("infection_status-polarized_violin_plot.png", p1)

# Interaction plot
p2 <- ggplot(fitting_results, aes(x=Infection_Status, y=art_ranks.Infection_Status.Polarized, color=Polarized, group=Polarized)) +
  stat_summary(fun=mean, geom="point", size=3) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  theme_minimal() +
  labs(title="Interaction Effect on Ranks",
       y="Aligned Rank Response")
ggsave("infection_status-polarized_interaction_plot.png", p2)

# Heatmap
rank_matrix <- dcast(fitting_results, Infection_Status ~ Polarized, 
                     value.var="art_ranks.Infection_Status.Polarized", fun.aggregate=mean)
p3 <- ggplot(melt(rank_matrix), aes(x=Infection_Status, y=variable, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title="Rank Distribution Heatmap", 
       y="Polarized",
       fill="Mean Aligned Rank")
ggsave("infection_status-polarized_heatmap.png", p3)



# 有意差が見つかった場合は事後分析（post-hoc）として対比較を行う
# 1. データをランク変換
# 2. 主効果や交互作用の影響を除去（align）
# 3. 変換されたデータに対して対応のあるt検定を実行 (半パラメトリックな検定)
# ただし、事後比較では (1|ID) が適切に扱われない可能性
# 同一個体から得られた複数のデータは独立でない
# ARTの事後比較は、データの独立性を前提としているので、今回は不適
# art_post_hoc.Infection_Status <- art.con(m, "Infection_Status")
# 結果をテキストファイルに保存
# sink("art_anova_post-hoc_Infection_Status.txt")
# print(art_post_hoc.Infection_Status)
# sink()

# 各個体のPolarizedデータの平均を計算
avg_by_id <- aggregate(Midpoint ~ ID + Infection_Status, 
                       data = fitting_results, 
                       FUN = mean)
# 平均値でMann-Whitney検定
result <- wilcox.test(Midpoint ~ Infection_Status, data = avg_by_id)

# 結果の保存
sink("infection_Mann-Whitney-U-test.txt")
print(result)
sink()

