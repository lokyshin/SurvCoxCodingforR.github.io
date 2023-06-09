window.R_strata_or_group = `# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

#安装必要的包
suppressPackageStartupMessages({
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(ggpubr)) install.packages("ggpubr")
  if (!require(survival)) install.packages("survival")
  if (!require(survminer)) install.packages("survminer")
})

#加载必要的包
library("ggplot2")
library("ggpubr")
library("survival")
library("survminer")

options(warn = -1)

cat("全人群与基于Strata_factor的生存分析（含单因素Cox回归分析）R统计程序\\nby Lokyshin\\n\\n")

#定义数据集，命名为df
df = data.frame(
    Patient_id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    Strata_factor = c(1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0),
    Time_to_event = c(2.6,9.3,8.2,16.8,7.2,17.4,14.6,18.9,9.2,16.4,2.8,4.8,13.4,12.3,19.7,14.4,8.7,13.6,0.2,7.8),
    Event = c(1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,0)
)
cat("研究数据集\\n")
print(df)

#数据整理
df$surv = with(df,Surv(Time_to_event, Event)) #添加标记
surv = df$surv

cat("\\n全人群分析\\n")
fit = survfit (surv ~ 1, df)
print(fit)
#print(summary(fit))

#取得总体结果，取小数点后3位
n = paste(fit$n, sep = " ")
median_total = round(surv_median(fit)[2]$median,3)
lcl95_total = round(surv_median(fit)[3]$lower,3)
ucl95_total = round(surv_median(fit)[4]$upper,3)
cl95_total = paste(lcl95_total, ucl95_total, sep = " - ")
#计算样本结果，取小数点后3位
median_sample = round(median(df$Time_to_event),3)
#计算样本中Time_to_event的95%可信区间

#构建结果
#median_ci = data.frame( Sample_Size = n, Median_total = median_total, conf_interval = cl95_total, Median_sample = median_sample)
#kable(median_ci, align = "l", column.width = c(10,10,10,10), col.names = c("研究的样本量","总体的中位数","总体可信区间","样本的中位数"))
cat(
    "\\n",
    "统计结果","\\n",
    "----------------------","\\n",
    "样本数量        ", as.character(n), "\\n",
    "----------------------","\\n",
    "总体","\\n",
    "中位数值        ", as.character(median_total), "\\n",
    "95%可信区间","\\n",
    "可信下限        ", as.character(lcl95_total), "\\n",
    "可信上限        ", as.character(ucl95_total), "\\n",
    "----------------------","\\n",
    "样本","\\n",
    "中位数值        ", as.character(median_sample), "\\n",
    "----------------------","\\n"
    )

km_curve = ggsurvplot(fit, df,
         pval = TRUE, #添加分层对比的p值
           pval.method = TRUE, #添加p值计算方法
           conf.int = TRUE, #添加可信区间
           risk.table = TRUE, #添加风险表
           tables.height = 0.2,#定义风险表
           surv.median.line = "hv", #添加中位值参考线
           ggtheme = theme_bw(), #更改ggplot2的主题
           palette = c("#9E006E"), #定义颜色
           title="Kaplan-Meier Curve", #定义图标题
           legend.title = "Patients", #定义图例名称
           ylab="Survival(%)",#定义纵坐标名称
           xlab = " Time (months)" #定义横坐标名称
           )       
km_curve

cat("\\n分层分析\\n")
fit_Strata_factor = survfit (surv ~ Strata_factor, df)
cat("分层的中位数\\n")
print(fit_Strata_factor)
cat("分层中位值的统计检验\\n")
survd = survdiff(surv ~ Strata_factor, df) #分层后的中位值统计检验，必须使用该格式，不能替换成fit_Strata_factor
survd
survp = surv_pvalue(fit_Strata_factor)
Strata_factor_kmp = survp$pval
n = survd$n[[dimnames(survd$n)[[1]][1]]] + survd$n[[dimnames(survd$n)[[1]][2]]]
Strata_factor0 = dimnames(survd$n)[[1]][1]
Strata_factor1 = dimnames(survd$n)[[1]][2]
samplesize_Strata_factor0 = survd$n[[dimnames(survd$n)[[1]][1]]]
samplesize_Strata_factor1 = survd$n[[dimnames(survd$n)[[1]][2]]]
obevents_Strata_factor0 = survd$obs[1]
obevents_Strata_factor1 = survd$obs[2]
expevents_Strata_factor0 = survd$exp[1]
expevents_Strata_factor1 = survd$exp[2]
table_fit_Strata_factor = summary(fit_Strata_factor)$table
median_Strata_factor0 = table_fit_Strata_factor[, "median"][1]
median_Strata_factor1= table_fit_Strata_factor[, "median"][2]
lcl_Strata_factor0 = table_fit_Strata_factor[, "0.95LCL"][1]
ucl_Strata_factor0 = table_fit_Strata_factor[, "0.95UCL"][1]
lcl_Strata_factor1 = table_fit_Strata_factor[, "0.95LCL"][2]
ucl_Strata_factor1 = table_fit_Strata_factor[, "0.95UCL"][1]
chisq_Strata_factor = survd$chisq
pval_Strata_factor = survp$pval
if (is.na(Strata_factor_kmp) || is.na(median_Strata_factor0) || is.na(median_Strata_factor1)) {
  km_curves_differ = "请谨慎看待组间p值，因为事件结局数量太少。"
} else {
  if (Strata_factor_kmp < 0.05) {
    km_curves_differ = "经统计学推断，分层后两组中位时间【有】统计学差异！"
  } else {
    km_curves_differ = "经统计学推断，分层后两组中位时间【无】统计学差异。"
  }
}

cat(
    "\\n",
    "生存函数相等性检验结果","\\n",
    "---------------------------", "\\n",
    "总样本量        ", as.character(n), "\\n",
    "---------------------------", "\\n",
    "分层因素        ", Strata_factor0,"\\n",
    "样本数量        ", samplesize_Strata_factor0,"\\n",
    "观测事件        ", as.character(round(obevents_Strata_factor0,3)), "\\n",
    "预期事件        ", as.character(round(expevents_Strata_factor0,3)), "\\n",
    "中位数值        ", as.character(round(median_Strata_factor0,3)), "\\n",
    "可信下限        ", as.character(round(lcl_Strata_factor0,3)),"*", "\\n",
    "可信上限        ", as.character(round(ucl_Strata_factor0,3)),"*", "\\n",
    "---------------------------", "\\n",
    "分层因素        ", Strata_factor1,"\\n",
    "样本数量        ", samplesize_Strata_factor1,"\\n",
    "观测事件        ", as.character(round(obevents_Strata_factor1,3)), "\\n",
    "预期事件        ", as.character(round(expevents_Strata_factor1,3)), "\\n",
    "中位数值        ", as.character(round(median_Strata_factor1,3)), "\\n",
    "可信下限        ", as.character(round(lcl_Strata_factor1,3)),"*", "\\n",
    "可信上限        ", as.character(round(ucl_Strata_factor1,3)),"*", "\\n",
    "---------------------------","\\n",
    "统计量","\\n",
    "卡方数值        ", as.character(round(chisq_Strata_factor,3)), "\\n",
    "检验概率        ", as.character(round(pval_Strata_factor,3)),"#", "\\n",
    "---------------------------","\\n",
    "*95%可信区间", "\\n",
    "#Log-rank p value", "\\n\\n",
    km_curves_differ, "\\n"
    )

km_curve_Strata_factor = ggsurvplot(fit_Strata_factor, df,
           pval = TRUE, #添加分层对比的p值
           pval.method = TRUE, #添加p值计算方法
           conf.int = FALSE, #可信区间是否显示
           risk.table = TRUE, #添加风险表
           risk.table.col = "Strata_factor", #根据分层更改风险表颜色
           linetype = "Strata_factor", # 根据分层更改线型
           surv.median.line = "hv", #同时显示垂直和水平参考线
           ggtheme = theme_bw(), #更改ggplot2的主题
           palette = c("#9E006E", "#2E9FDF"),#定义颜色
           title="Kaplan-Meier Curve by Strata_factor", #定义图标题
           ylab="Survival(%)",#定义纵坐标名称
           xlab = " Time (months)", #定义横坐标名称
           )       
km_curve_Strata_factor


cat("\\nCox回归的单因素分析\\n")
Strata_factor_Cox = coxph(surv ~ Strata_factor, df) #计算Strata_factor的Cox
Strata_factor_Sum = summary(Strata_factor_Cox) #总结一下
cat("分层单因素分析结果\\n")
Strata_factor_Cox

nevent = Strata_factor_Cox$nevent
lrchi2 = Strata_factor_Sum$logtest[1]
lrpval = Strata_factor_Sum$logtest[3]
hr_cox = Strata_factor_Sum$coefficients[2]
z_cox = Strata_factor_Sum$coefficients[4]
p_cox = Strata_factor_Sum$coefficients[5]
lcl_cox = Strata_factor_Sum$conf.int[3]
ucl_cox = Strata_factor_Sum$conf.int[4]
if (is.na(p_cox)) {
  unicox = "无法进行单因素Cox回归分析，可能是事件结局数量太少。"
} else {
  if (p_cox < 0.05) {
    unicox = "该因素【是】事件结局的影响因素，应进入多因素分析是否独立影响因素！"
  } else {
    unicox = "该因素【并非】事件结局的影响因素。"
  }
}

cat(
    "\\n",
    "Cox回归检验结果","\\n",
    "---------------------------", "\\n",
    "样本数量        ", as.character(Strata_factor_Cox$n), "\\n",
    "事件数量        ", as.character(nevent), "\\n",
    "---------------------------", "\\n",
    "似然卡方        ", as.character(round(lrchi2,3)),"*1", "\\n",
    "检验概率        ", as.character(round(lrpval,3)),"*2", "\\n",
    "---------------------------", "\\n",
    "风险比值        ", as.character(round(hr_cox,3)),"\\n",
    "检验变量        ", as.character(round(z_cox,3)),"*3", "\\n",
    "检验概率        ", as.character(round(p_cox,3)),"*4", "\\n",
    "可信下限        ", as.character(round(lcl_cox,3)), "\\n",
    "可信上限        ", as.character(round(ucl_cox,3)), "\\n",
    "---------------------------","\\n",
    "*1 最大似然估计卡方值", "\\n",
    "*2 最大似然估计p值", "\\n",
    "*3 z值", "\\n",
    "*4 单因素Cox回归的p值", "\\n\\n",
    unicox, "\\n"
    )

cat("\\n分析结论\\n")
cat(
    "汇总报告","\\n",
    "---------------------------", "\\n",
    "样本数量        ", as.character(n), "\\n",
    "事件数量        ", as.character(nevent), "\\n",
    "---------------------------", "\\n",
    "---------------------------", "\\n",
    "总体情况","\\n",
    "中位数值        ", as.character(median_total), "\\n",
    "可信下限        ", as.character(lcl95_total), "*1", "\\n",
    "可信上限        ", as.character(ucl95_total), "*1", "\\n",
    "---------------------------", "\\n",
    "---------------------------", "\\n",
    "分层分析","\\n",
    "分层因素        ", Strata_factor0,"\\n",
    "样本数量        ", samplesize_Strata_factor0,"\\n",
    "中位数值        ", as.character(round(median_Strata_factor0,3)), "\\n",
    "可信下限        ", as.character(round(lcl_Strata_factor0,3)),"*1", "\\n",
    "可信上限        ", as.character(round(ucl_Strata_factor0,3)),"*1", "\\n",
    "分层因素        ", Strata_factor1,"\\n",
    "样本数量        ", samplesize_Strata_factor1,"\\n",
    "中位数值        ", as.character(round(median_Strata_factor1,3)), "\\n",
    "可信下限        ", as.character(round(lcl_Strata_factor1,3)),"*1", "\\n",
    "可信上限        ", as.character(round(ucl_Strata_factor1,3)),"*1", "\\n",
    "统计量","\\n",
    "卡方数值        ", as.character(round(chisq_Strata_factor,3)), "\\n",
    "检验概率        ", as.character(round(pval_Strata_factor,3)),"*2", "\\n",
    "---------------------------","\\n",
    "---------------------------","\\n",
    "单因素Cox回归","\\n",
    "风险比值        ", as.character(round(hr_cox,3)),"\\n",
    "检验概率        ", as.character(round(p_cox,3)),"*3", "\\n",
    "可信下限        ", as.character(round(lcl_cox,3)), "\\n",
    "可信上限        ", as.character(round(ucl_cox,3)), "\\n",
    "---------------------------","\\n",
    "*1 95%可信区间", "\\n",
    "*2 Log-rank p value", "\\n",
    "*3 单因素Cox回归的p值", "\\n\\n"
    )

cat(
    "总结报告","\\n",
    "1 如全人群推断总体的Kaplan-Meier曲线所示，本研究总体预后时间（Time_to_event）为：", as.character(median_total),"[", as.character(lcl95_total), " - ", as.character(ucl95_total), "] 月。","\\n",
    "2 按照本研究Strata分层因素，两组预后时间（Time_to_event）分别为：", "\\n",
    "    ", Strata_factor0, ": ", as.character(round(median_Strata_factor0,3)), "[", as.character(round(lcl_Strata_factor0,3)), " - ", as.character(round(ucl_Strata_factor0,3)), "] 月;","\\n",
    "    ", Strata_factor1, ": ", as.character(round(median_Strata_factor1,3)), "[", as.character(round(lcl_Strata_factor1,3)), " - ", as.character(round(ucl_Strata_factor1,3)), "] 月;","\\n",
    "    ", km_curves_differ, "\\n",
    "3 对Strata分层因素进行单因素Cox回归时发现，", unicox, "\\n"
    )

cat("\\n", "感谢您使用本程序进行统计分析，再见。", "\\n")
`;

window.R_cox = `# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

#安装必要的包
suppressPackageStartupMessages({
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(ggpubr)) install.packages("ggpubr")
  if (!require(survival)) install.packages("survival")
  if (!require(survminer)) install.packages("survminer")
})

#加载必要的包
library("ggplot2")
library("ggpubr")
library("survival")
library("survminer")

options(warn = -1)

cat("Cox回归分析R统计程序\\nby Lokyshin\\n\\n")

#定义数据集，命名为df
df = data.frame(
    Patient_id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    Strata_factor = c(1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0),
    Time_to_event = c(2.6,9.3,8.2,16.8,7.2,17.4,14.6,18.9,9.2,16.4,2.8,4.8,13.4,12.3,19.7,14.4,8.7,13.6,0.2,7.8),
    Event = c(1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,0)
)
cat("研究数据集\\n")
print(df)

#数据整理
df$surv = with(df,Surv(Time_to_event, Event)) #添加标记
surv = df$surv

cat("\\n全人群分析\\n")
fit = survfit (surv ~ 1, df)
print(fit)
#print(summary(fit))

#取得总体结果，取小数点后3位
n = paste(fit$n, sep = " ")
median_total = round(surv_median(fit)[2]$median,3)
lcl95_total = round(surv_median(fit)[3]$lower,3)
ucl95_total = round(surv_median(fit)[4]$upper,3)
cl95_total = paste(lcl95_total, ucl95_total, sep = " - ")
#计算样本结果，取小数点后3位
median_sample = round(median(df$Time_to_event),3)
#计算样本中Time_to_event的95%可信区间
##计算可信区间，由于在线R没有此包，代码注释pr掉
#lcl95_sample = round(CI(df$Time_to_event)[3],3)由于在线R没有此包，代码注释掉
#ucl95_sample = round(CI(df$Time_to_event)[1],3)由于在线R没有此包，代码注释掉
#构建结果
#median_ci = data.frame( Sample_Size = n, Median_total = median_total, conf_interval = cl95_total, Median_sample = median_sample)
#kable(median_ci, align = "l", column.width = c(10,10,10,10), col.names = c("研究的样本量","总体的中位数","总体可信区间","样本的中位数"))
cat(
    "\\n",
    "统计结果","\\n",
    "----------------------","\\n",
    "样本数量        ", as.character(n), "\\n",
    "----------------------","\\n",
    "总体","\\n",
    "中位数值        ", as.character(median_total), "\\n",
    "95%可信区间","\\n",
    "可信下限        ", as.character(lcl95_total), "\\n",
    "可信上限        ", as.character(ucl95_total), "\\n",
    "----------------------","\\n",
    "样本","\\n",
    "中位数值        ", as.character(median_sample), "\\n",
    "----------------------","\\n"
    )


cat("\\nCox回归分析结果\\n")
# 进行coxph分析并获取结果摘要
Strata_Cox = coxph(surv ~ Strata_factor, df) #计算Strata_factor的Cox
Strata_Sum = summary(Strata_Cox) #总结一下
print(Strata_Sum)

# 提取摘要中的系数表
coef_table = Strata_Sum$coefficients

# 查找p值小于0.05的行
significant_results = coef_table[coef_table[, "Pr(>|z|)"] < 0.05, ]

# 如果找到了p值小于0.05的变量，则打印它们的名称和p值；否则，打印“没有找到”
if (nrow(significant_results) > 0) {
  for (i in row.names(significant_results)) {
    cat(paste("变量", i, "的p值为", round(significant_results[i, "Pr(>|z|)"], 3), "是结局事件的风险因素，建议放入多因素分析Cox回归模型\\n"))
  }
} else {
  cat("没有找到p值小于0.05的变量（结局事件的风险因素）。\\n")
}

cat("\\n感谢您使用本程序进行统计分析，再见。", "\\n")
`;

window.R_des_continuous= `# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

# 安装所需的包
suppressPackageStartupMessages({
  if(!require(tidyverse)) install.packages("tidyverse")
  if(!require(dplyr)) install.packages("dplyr")
  if(!require(broom)) install.packages("broom")
  if (!require(purrr)) install.packages("purrr")
})

# 加载所需的包
library(tidyverse)
library(dplyr)
library(broom)
library(purrr)

options(warn = -1)

cat("连续变量统计分析的R程序\\nby Lokyshin\\n")

#定义数据集，命名为df
df = data.frame(
    RDCvariate = c(3.4,2.5,6.7,8,5.7,3.2,3.1,9.9,10,11,12,4.5,6.7,8.3,6.9,4.7,6.9,2.1,12.3,7.4),
    RDCgroup = c(1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,1,0)
)
cat("\\n研究数据集\\n")
print(df)

# 检查RDCgroup中的唯一值的数量
unique_groups = length(unique(df$RDCgroup))

summary_stats = df %>%
  group_by(RDCgroup) %>%
  summarize(
    Mean = round(mean(RDCvariate),3),
    SD = round(sd(RDCvariate),3),
    SE = round(sd(RDCvariate) / sqrt(n()),3),
    Q1 = round(quantile(RDCvariate, 0.25),3),
    Median = round(median(RDCvariate),3),
    Q3 = round(quantile(RDCvariate, 0.75),3),
    IQR = round(IQR(RDCvariate),3),
    \`95%CI\` = broom::tidy(t.test(RDCvariate)) %>% 
                summarise(CI = paste(round(conf.low, 3), "to", round(conf.high, 3))) %>%
                pull(CI),
    .groups = 'drop'
  )

# 转换为数据框
summary_stats = as.data.frame(summary_stats)

cat("\\nRDCvariate按RDCgroup分组统计描述", "\\n")
print(summary_stats)

cat("\\n下面将进行两组或多组的均数检验。这通常用于比较两组或多组的基线是否平齐", "\\n")
cat("以下是两组或多组的均数比较流程：\\n",
    "1. 独立性假设\\n",
    "2. 正态分布假设\\n",
    "\t如果不符合正态分布，采用非参数检验，本程序采用Kruskal-Wallis检验\\n",
    "\t如果符合正态分布，进入下一步\\n",
    "3. 方差齐性假设\\n",
    "\t如果方差齐，且为两组，进行t检验；如果为三组或更多，进行F检验（单因素方差分析）\\n",
    "\t如果方差不齐，且为两组，进行Welch's t检验；如果为三组或更多，本程序采用Kruskal-Wallis检验\\n",
    sep = ""
)

# 进行Shapiro-Wilk检验(正态性检验)
normality_test_results = tapply(df$RDCvariate, df$RDCgroup, shapiro.test)

# 创建一个数据框来保存结果
normality_test_results_df = data.frame(
  RDCgroup = names(normality_test_results),
  W = sapply(normality_test_results, function(x) x$statistic),
  p.value = sapply(normality_test_results, function(x) x$p.value)
)
cat("\\nRDCvariate按RDCgroup分组的正态性检验结果", "\\n\\n")
print(normality_test_results_df)

# 进行方差齐性检验
bartlett_result = bartlett.test(RDCvariate ~ RDCgroup, data = df)
bartlett_p = bartlett_result$p.value

#Kruskal-Wallis检验
kruskal_result = kruskal.test(RDCvariate ~ RDCgroup, data = df)
kruskal_p = kruskal_result$p.value

#Welch's t检验
t_test_result = t.test(RDCvariate ~ RDCgroup, data = df, var.equal = FALSE)
t_test_result_p = t_test_result$p.value  

#F检验（单因素方差分析）
aov_result = aov(RDCvariate ~ as.factor(RDCgroup), data = df)
aovm = summary(aov_result)
aovm_p = aovm[[1]]$\`Pr(>F)\`[1]

#成组t检验
t_test_result = t.test(RDCvariate ~ RDCgroup, data = df, var.equal = TRUE)
t_test_result_p = t_test_result$p.value

#如果不符合正态分布，采用非参数检验，本程序采用Kruskal-Wallis检验
if(any(normality_test_results_df$p.value <= 0.05)) {
  cat("\\nKruskal-Wallis检验结果", "\\n")
  print(kruskal_result)
  if(kruskal_p<0.05){
    cat("
    组间基线情况：\\n分组存在不符合正态分布的情况。\\n当前基线不平齐。基线描述建议采用中位数和四分位间距。")}else{cat("基线情况：\\n分组存在不符合正态分布的情况。\\n当前基线平齐。基线描述建议采用中位数和四分位间距。
    ")
  }
}

#如果符合正态分布,方差不齐，为三组或更多，本程序采用Kruskal-Wallis检验
if((any(normality_test_results_df$p.value > 0.05)) && (bartlett_p <= 0.05) && (unique_groups != 2)) {
  cat("\\nRDCvariate按RDCgroup分组的方差齐性检验结果", "\\n")
  print(bartlett_result)
  cat("Kruskal-Wallis检验结果", "\\n")
  print(kruskal_result)
  if(kruskal_p<0.05){
    cat("
    组间基线情况：\\n符合正态分布，但方差不齐。\\n当前基线不平齐。基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。")}else{cat("基线情况：\\n符合正态分布，但方差不齐。\\n当前基线平齐。基线描述基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。
    ")
  }
}

#如果符合正态分布,方差不齐，且为两组，进行Welch's t检验
if((any(normality_test_results_df$p.value > 0.05)) && (bartlett_p <= 0.05) && (unique_groups == 2)) {
  cat("\\nRDCvariate按RDCgroup分组的方差齐性检验结果", "\\n")
  print(bartlett_result)
  cat("Welch's t检验结果", "\\n")
  print(t_test_result)
  if(t_test_result_p<0.05){
    cat("
    组间基线情况：\\n符合正态分布，但方差不齐。\\n当前基线不平齐。基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。")}else{cat("基线情况：\\n符合正态分布，但方差不齐。\\n当前基线平齐。基线描述基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。
    ")
  }
}

#如果符合正态分布,方差齐，为三组或更多，进行F检验（单因素方差分析）
if((any(normality_test_results_df$p.value > 0.05)) && (bartlett_p > 0.05) && (unique_groups != 2)) {
  cat("\\nRDCvariate按RDCgroup分组的方差齐性检验结果", "\\n")
  print(bartlett_result)
  cat("F检验（单因素方差分析）结果", "\\n")
  print(t_test_result)
  if(aovm_p<0.05){cat("
    组间基线情况：\\n符合正态分布，方差齐。\\n当前基线不平齐。基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。")}else{cat("基线情况：\\n符合正态分布，方差齐。\\n当前基线平齐。基线描述基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。
    ")
  }
}

#如果符合正态分布,方差齐，且为两组，进行t检验
if((any(normality_test_results_df$p.value > 0.05)) && (bartlett_p > 0.05) && (unique_groups == 2)) {
  cat("\\nRDCvariate按RDCgroup分组的方差齐性检验结果", "\\n")
  print(bartlett_result)
  cat("t检验结果", "\\n")
  print(t_test_result)
  if(t_test_result_p<0.05){
    cat("
    组间基线情况：\\n符合正态分布，方差齐。\\n当前基线不平齐。基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。")}else{cat("基线情况：\\n符合正态分布，方差齐。\\n当前基线平齐。基线描述基线描述可采用均数和95%可信区间（生存分析用）或均数与标准差（一般）。
    ")
  }
}

# 绘图箱图
# 将RDCgroup列转换为因子类型
df$RDCgroup <- as.factor(df$RDCgroup)
# 创建箱图
ggplot(df, aes(x=RDCgroup, y=RDCvariate, fill=RDCgroup)) +
    geom_boxplot() +
    labs(title="Boxplot of RDCvariate by Group", 
         x="Group", 
         y="RDCvariate",
         fill="Group") +
    theme_minimal()
cat("\\n\n您可在下方查看分组数据的箱图。", "\\n")    

cat("\\n\\n感谢您使用本程序进行统计分析，再见。", "\\n")
`;

window.R_des_count= `# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

options(warn = -1)

cat("计数资料统计描述的R程序\\nby Lokyshin\\n")

#定义数据集，命名为df
df = data.frame(
    RD_characteristic = c(1,1,0,0,1,0,1,1,1,0,1,0,1,0,1,0),
    RD_classfication = c(0,1,1,0,0,1,0,1,1,0,1,0,1,0,1,0)
)

# 创建行列表
cross_table = xtabs(~ RD_characteristic + RD_classfication, data = df)
rownames(cross_table) <- c("RD_characteristic 0", "RD_characteristic 1")  # 调整行名称
colnames(cross_table) = c("RD_classfication 0", "RD_classfication 1")  # 调整列名称

#下面开计算格子的理论频数等信息，用于判断对应的统计方法
#首先，我们定义计算理论频数的函数
calc_expected_frequency = function(TableData){
  row_totals = rowSums(TableData)
  col_totals = colSums(TableData)
  grand_total = sum(TableData)
  
  expected_frequency = outer(row_totals, col_totals) / grand_total
  return(expected_frequency)
}

# 定义函数判断是否为四格表
is_2by2 = function(TableData) {
  return(dim(TableData)[1] == 2 & dim(TableData)[2] == 2)
}

# 定义函数进行相应的检验
test_table <- function(TableData){
  expected_frequency <- calc_expected_frequency(TableData)
  min_frequency <- min(expected_frequency)
  total_sample <- sum(TableData)
  
  if(is_2by2(TableData)){
    if(min_frequency >= 5 && total_sample >= 40){
      result <- chisq.test(TableData, correct = FALSE)
      cat("\\n表格是四格表，所有的理论频数≥5，并且总样本量≥40，因此使用Pearson卡方进行检验。\\n")
      cat("\\nThe table is a 2x2 table, all theoretical frequencies T≥5 and total sample size n≥40, so Pearson's chi-square test is used.\\n")
    } else if (min_frequency >= 1 && total_sample >= 40) {
      result <- chisq.test(TableData, correct = TRUE)
      cat("\\n表格是四格表，1≤理论频数<5，并且总样本量因此使用连续性校正的卡方进行检验\\n")
      cat("\\nThe table is a 2x2 table, the theoretical frequency T<5 but T≥1, and the total sample size n≥40, so the chi-square test with continuity correction is used.\\n")
    } else {
      result <- fisher.test(TableData)
      cat("\\n表格是四格表，但有理论频数<<1 或总样本量<40，因此使用 Fisher’s 检验。\\n")
      cat("\\nThe table is a 2x2 table, but there are theoretical numbers T<1 or total sample size n<40, so Fisher's test is used.\\n")
    }
  } else {
    if(sum(expected_frequency < 5) <= 0.2 * length(TableData) && min_frequency >= 1){
      result <- chisq.test(TableData, correct = FALSE)
      cat("\\n表格不是四格表，但 R×C 表中论频数<5的格子不超过1/5并且任何格子没有小于1的理论频数，因此使用 Pearson 卡方进行检验。\\n")
      cat("\\nThe table is not a 2x2 table, but in the R×C table, the cells with theoretical frequencies less than 5 do not exceed 1/5, and no cells have theoretical frequencies less than 1, so Pearson's chi-square test is used.\\n")
    } else {
      result <- fisher.test(TableData)
      cat("\\n表格不是四格表，并且不满足使用 Pearson 卡方进行检验的条件，因此使用 Fisher’s 检验。\\n")
      cat("\\nThe table is not a 2x2 table and does not meet the conditions for using Pearson's chi-square test, so Fisher's test is used.\\n")
    }
  }
  
  # 提取p值并进行判断
  p_value <- result$p.value
  if(p_value < 0.05){
    cat("\\np =", round(p_value,3), "差异有统计学意义。请查看下方检验结果：\\n")
  } else {
    cat("\\np =", round(p_value,3), "，没有统计学差异。请查看下方检验结果：\\n")
  }
  
  return(result)
}

# 打印行列表
cat("\\nRD_characteristic X RD_classfication\\n行列表\\n")
cross_table

# 打印判断标准
cat("1. 计算期望频率和最小频率，以及总样本量\\n")
cat("\\n2. 检查表格类型\\n")
cat("2.1 如果是2x2表格，执行以下操作：\\n")
cat("2.1.1 如果所有的理论频数≥5，并且总样本量≥40，进行pearson卡方检验\\n")
cat("2.1.2 如果1≤理论频数<5，并且总样本量≥40，进行连续性校正的卡方检验\\n")
cat("2.1.3 如果理论频数<1 或总样本量<40，进行Fisher’s 检验\\n")
cat("\\n3. 如果不是2x2表格，执行以下操作：\\n")
cat("3.1 如果R×C表中理论频数<5的格子不超过总格子数的1/5，并且所有格子的理论频数≥1，进行pearson卡方检验\\n")
cat("3.2 如果不满足上述条件，进行Fisher’s 检验\\n")

# 对cross_table进行检验
test_table(cross_table)

cat("\\n感谢您使用本程序进行统计分析，再见。", "\\n")
`;


window.R_des_ranked= `# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

options(warn = -1)

cat("等级资料统计描述的R程序\\nby Lokyshin\\n")

#定义数据集，命名为df
df = data.frame(
    RD_grade_value = c(1,1,0,0,1,0,1,1,1,0,1,0,1,0,1,0),
    RD_grade_group = c(0,1,1,0,0,1,0,1,1,0,1,0,1,0,1,0)
)

# 创建行列表
cross_table = xtabs(~ RD_grade_value + RD_grade_group, data = df)
rownames(cross_table) <- c("RD_grade_value 0", "RD_grade_value 1")  # 调整行名称
colnames(cross_table) = c("RD_grade_group 0", "RD_grade_group 1")  # 调整列名称

#Kruskal-Wallis检验
kruskal_result = kruskal.test(RD_grade_value ~ RD_grade_group, data = df)
# 提取p值
kruskal_p = kruskal_result$p.value
  
# 打印行列表
cat("\\nRD_characteristic X RD_classfication\\n行列表\\n")
cross_table

# 打印检验结果
kruskal_result

# 进行结果判断
if(kruskal_p < 0.05){
  cat("p值<0.05，差异有统计学意义。\\n")
} else {
  cat("p =", round(kruskal_p,3), "，没有统计学差异。\\n")
}

cat("\\n感谢您使用本程序进行统计分析，再见。", "\\n")
`;

window.R_wfplot= `
# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html


# 安装所需的包
suppressPackageStartupMessages({
  if(!require(ggplot2)) install.packages("ggplot2")
})

# 加载所需的包
library(ggplot2)

options(warn = -1)

cat("绘制肿瘤瀑布图的R程序\\nby Lokyshin\\n")

# 定义数据集，命名为df
df = data.frame(
    Patient_id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),
    Time_to_event = c(2.6,9.3,8.2,12.8,7.2,17.4,11.6,17.9,9.2,10.4,2.8,4.8,16.1,12.3,11.7,9.8,12.3,13.4,16.7,18.9,22),
    Event = c(1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
    Groups = c(1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0),
    Factor = c(0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,1),
    Response_Rate = c(26,22,20,19,15,6,0,0,-8,-15,-20,-24,-28,-30,-37,-42,-28,-30,-37,-42,-50),
    Response_evaluation = c('NE','PD','SD','SD','SD','PD','SD','SD','SD','SD','SD','SD','SD','PR','PR','PR','PR','PR','PD','PR','PR')
)

# 按照Response从大到小排序
df = df[order(-df$Response_Rate),]

# 定义颜色映射
color_mapping = c('CR' = "#CCEAFC", #淡蓝，#淡粉#FCD9FC
                   'PR' = "#D2E8E3", #淡绿，或者#DDFCCC
                   'SD' = "#FFE8B5", #淡黄，或者#F2E3D5 或者#FCE8CC
                   'PD' = "#F4E2DE", #淡红
                   'NE' = "#5D6A73", #淡褐，#深褐#594438 
                   'NA' = "#C2C3BE") #淡灰
                  

# 计算y轴最大值
y_max = ifelse(max(df$Response_Rate) < 30, 40, 
                ifelse(max(df$Response_Rate) + 20 < 100, max(df$Response_Rate) + 20, 
                       max(df$Response_Rate)))

# 计算y轴最小值
if (min(df$Response_Rate) > -30) {
  y_min = -40
} else if (-40 <= min(df$Response_Rate) && min(df$Response_Rate) < -30) {
  y_min = -50
} else if (-50 <= min(df$Response_Rate) && min(df$Response_Rate) < -40) {
  y_min = -60
} else if (-60 <= min(df$Response_Rate) && min(df$Response_Rate) < -50) {
  y_min = -80
} else if (min(df$Response_Rate) < -60) {
  y_min = -100
}

# 使用颜色映射对Response_evaluation进行映射
colors = color_mapping[df$Response_evaluation]

# 创建透明度向量，group为1的设为0.6，其余设为1
alpha = ifelse(df$Groups == 1, 0.6, 1)

# 预设分组图标
group_shapes = c(1, 2, 0, 6, 5)
# 0 空心方框
# 1 空心圆圈
# 2 空心三角
# 3 十字符号
# 4 交叉符号
# 5 旋转45度的空心正方形
# 6 空心倒三角

# 增加图形下方空间放置x轴标题
#par(mar=c(5, 0, 0, 0))

# 绘制条形图
bp = barplot(df$Response_Rate, names.arg=df$Patient_id, col=colors, ylim=c(y_min, y_max),
             xaxt="n", main="Waterfall plot of best overall response", ylab="Response Rate (%)",
             cex.main=1.3, cex.sub=1.3, cex.lab=1, cex.names=0.8, cex.axis=0.8, border='#D9CDBF')


# 添加x轴
axis(1, at=bp, labels=1:length(df$Response_Rate))

# 添加水平线和注释
abline(h=20, lty=2, col = rgb(0, 0, 0, alpha = 0.2))#设置透明度80%
abline(h=-30, lty=2, col = rgb(0, 0, 0, alpha = 0.2))#设置透明度80%

# 绘制图例
legend_labels = names(color_mapping)[names(color_mapping) %in% unique(df$Response_evaluation)]
legend("bottomleft", inset=c(0.1,0.1), legend=legend_labels, 
       fill=color_mapping[legend_labels], 
       title="Response", bty="n", cex=0.8, border='#D9CDBF')

# 初始化绘制组别和组别图例开关
addgroup <- 1
# 绘制组别和组别图例开关（待替换）
addgroup = 1


if (addgroup == 1) {
  # 在条形图上添加分组图形
  points(bp, ifelse(df$Response_Rate > 0, df$Response_Rate + 3, 3), pch = group_shapes[df$Groups+1], cex = 1.2)

  # 添加关于组别的图例
  legend("bottomleft", inset=c(0.25,0.1), legend=paste("Groups", 0:(length(unique(df$Groups))-1)), 
         pch = group_shapes[1:length(unique(df$Groups))], 
         title="Groups", bty="n", cex=0.8,border='#D9CDBF')
}

title(xlab="Patients id")

# 出图
bp

cat("\\n感谢您使用本程序进行统计分析，再见。", "\\n")
`;

window.R_swplot= `

# Code by Lokyshin (Lokyshin.net@202306)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html


# 安装所需的包
suppressPackageStartupMessages({
  if(!require(ggplot2)) install.packages("ggplot2")
})

# 加载所需的包
library(ggplot2)

options(warn = -1)

cat("绘制肿瘤游泳图的R程序\\nby Lokyshin\\n")

# 定义数据集，命名为df
df = data.frame(
    Patient_id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),
    Time_to_event = c(2.6,9.3,8.2,12.8,7.2,17.4,11.6,17.9,9.2,10.4,2.8,4.8,16.1,12.3,11.7,9.8,12.3,13.4,16.7,18.9,22),
    Event = c(1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
    Groups = c(1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0),
    Factor = c(0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,1),
    Response_Rate = c(26,22,20,19,15,6,0,0,-8,-15,-20,-24,-28,-30,-37,-42,-28,-30,-37,-42,-50),
    Response_evaluation = c('NE','PD','SD','SD','SD','PD','SD','SD','SD','SD','SD','SD','SD','PR','PR','PR','PR','PR','PD','PR','PR')
)

# 按照Response从大到小排序
df = df[order(df$Time_to_event),]

# 定义颜色映射
color_mapping = c('CR' = "#D9E7F4", 
                   'PR' = "#9FCCC2", 
                   'SD' = "#F4E0BD", 
                   'PD' = "#D8E0F0", 
                   'NE' = "#C2C3BE", 
                   'NA' = "#5D6A73") 
                  
# 使用颜色映射对Response_evaluation进行映射
colors = color_mapping[df$Response_evaluation]

# 创建透明度向量，group为1的设为0.6，其余设为1
alpha = ifelse(df$Groups == 1, 0.6, 1)

# 计算x轴最大值
x_max = max(df$Time_to_event)*1.1

# 计算x轴最小值
x_min = 0

# 预设分组图标
group_shapes = c(1, 2, 0, 6, 5)
# 0 空心方框
# 1 空心圆圈
# 2 空心三角
# 3 十字符号
# 4 交叉符号
# 5 旋转45度的空心正方形
# 6 空心倒三角

# 预设结局事件图标
Ev_shapes = c(0, 15)

# 为事件点创建颜色映射
event_colors = c("#594438", "#594438")  # 这里根据图标进行区分，因此都设置了相同的颜色。

# 0 空心方框
# 15 实心矩形

# 绘制条形图
par(las=2, mar=c(7, 5, 2, 2))

sw = barplot(df$Time_to_event, col=colors, horiz=TRUE, xlim=c(x_min, x_max),
             xaxt="n", yaxt="n",
             cex.main=1.3, cex.sub=1.3, cex.lab=1, cex.names=0.8, cex.axis=0.8, border='#F3F3F3')

title(main="Swimmer plot of all patients", xlab="Duration(months)", ylab="Patients id", cex.main=1.1)
axis(1, at=seq(x_min, x_max), las=1, tck=-0.01, cex.axis=0.8)
axis(2, at=sw, labels=df$Patient_id, tick=FALSE, cex.axis=0.8)

df$Event = as.factor(df$Event)
# 在条形图上绘制Event
for (i in 1:length(sw)) {
  # 寻找对应的Event图形
  Ev_shape = Ev_shapes[which(levels(df$Event) == df$Event[i])]
  
  # 寻找对应的颜色
  Ev_color = event_colors[which(levels(df$Event) == df$Event[i])]

  # 绘制Event图形
  points(df$Time_to_event[i], sw[i], pch=Ev_shape, col=Ev_color)
}

# 绘制图例
legend_labels = names(color_mapping)[names(color_mapping) %in% unique(df$Response_evaluation)]
legend("bottomright", inset=c(0.1,0.1), legend=legend_labels, 
       fill=color_mapping[legend_labels], 
       title="Response", bty="n", cex=0.8, border='#F3F3F3')

# 添加Event图例
legend("bottomright", inset=c(0.25,0.1), legend=paste("Event", 0:(length(unique(df$Event))-1)), 
       pch = Ev_shapes[1:length(unique(df$Event))], 
       col = event_colors[1:length(unique(df$Event))], 
       title="Event", bty="n", cex=0.8,border='#594438')

# 绘制全人群游泳图
sw

# 初始化绘制组别和组别图例开关
addgroup <- 1
# 绘制组别和组别图例开关（待替换）
addgroup = 1

if (addgroup == 1) {

  # 分组绘制魅族的游泳图
  # 将数据拆分成不同的 df
  grouped_data = split(df, df$Groups)

  # 循环中绘制每个 df 的条形图
  for(i in seq_along(grouped_data)) {

    # 绘制条形图
    par(las=2, mar=c(7, 5, 2, 2))
    
    sw = barplot(grouped_data[[i]]$Time_to_event, col=colors, horiz=TRUE, xlim=c(x_min, x_max),
                 xaxt="n", yaxt="n",
                 cex.main=1.3, cex.sub=1.3, cex.lab=1, cex.names=0.8, cex.axis=0.8, border='#F3F3F3')

    title(main=paste("Swimmer plot of Group", i), xlab="Duration(months)", ylab="Patients id", cex.main=1.1) 
    axis(1, at=seq(x_min, x_max), las=1, tck=-0.01, cex.axis=0.8)
    axis(2, at=sw, labels=grouped_data[[i]]$Patient_id, tick=FALSE, cex.axis=0.8)

    grouped_data[[i]]$Event = as.factor(grouped_data[[i]]$Event)
    
    # 在条形图上绘制Event图形
    for (j in 1:length(sw)) {
      # 寻找对应的Event图形
      Ev_shape = Ev_shapes[which(levels(grouped_data[[i]]$Event) == grouped_data[[i]]$Event[j])]
      
      # 绘制Event图形
      points(grouped_data[[i]]$Time_to_event[j], sw[j], pch=Ev_shape, col=Ev_color)
    }

    # 绘制图例
    legend_labels = names(color_mapping)[names(color_mapping) %in% unique(grouped_data[[i]]$Response_evaluation)]
    legend("bottomright", inset=c(0.1,0.1), legend=legend_labels, 
           fill=color_mapping[legend_labels], 
           title="Response", bty="n", cex=0.8, border='#F3F3F3')

    # 添加Event图例
    legend("bottomright", inset=c(0.25,0.1), legend=paste("Event", 0:(length(unique(df$Event))-1)), 
           pch = Ev_shapes[1:length(unique(df$Event))], 
           col = event_colors[1:length(unique(df$Event))], 
           title="Event", bty="n", cex=0.8,border='#594438')

  }
}

cat("\\n感谢您使用本程序进行统计分析，再见。", "\\n")

`;
