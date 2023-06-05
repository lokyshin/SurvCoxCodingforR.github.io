window.Rsurv_cox = `# Code by Lokyshin (Lokyshin.net@202306051232)
# Reference
# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer
# https://cran.r-project.org/web/packages/survminer/index.html

library("ggplot2")
library("ggpubr")
library("survival")
library("survminer")
#library("Rmic") #用于计算正态分布的可信区间。在线R无此包，因此注释掉
options(warn = -1)

#定义数据集，命名为df
df = data.frame(
    Patient_id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    Strata = c(1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0),
    Time_to_event_month = c(2.6,9.3,8.2,16.8,7.2,17.4,14.6,18.9,9.2,16.4,2.8,4.8,13.4,12.3,19.7,14.4,8.7,13.6,0.2,7.8),
    Event = c(1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,0)
)
cat("研究数据集\n")
print(df)

#数据整理
df$surv = with(df,Surv(Time_to_event_month, Event)) #添加标记
surv = df$surv

cat("\n全人群分析\n")
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
median_sample = round(median(df$Time_to_event_month),3)
#计算样本中Time_to_event_month的95%可信区间
##计算可信区间，由于在线R没有此包，代码注释pr掉
#lcl95_sample = round(CI(df$Time_to_event_month)[3],3)由于在线R没有此包，代码注释掉
#ucl95_sample = round(CI(df$Time_to_event_month)[1],3)由于在线R没有此包，代码注释掉
#构建结果
#median_ci = data.frame( Sample_Size = n, Median_total = median_total, conf_interval = cl95_total, Median_sample = median_sample)
#kable(median_ci, align = "l", column.width = c(10,10,10,10), col.names = c("研究的样本量","总体的中位数","总体可信区间","样本的中位数"))
cat(
    "\n",
    "统计结果","\n",
    "----------------------","\n",
    "样本数量        ", as.character(n), "\n",
    "----------------------","\n",
    "总体","\n",
    "中位数值        ", as.character(median_total), "\n",
    "95%可信区间","\n",
    "可信下限        ", as.character(lcl95_total), "\n",
    "可信上限        ", as.character(ucl95_total), "\n",
    "----------------------","\n",
    "样本","\n",
    "中位数值        ", as.character(median_sample), "\n",
    "----------------------","\n"
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

cat("\n分层分析\n")
fit_Strata = survfit (surv ~ Strata, df)
cat("分层的中位数\n")
print(fit_Strata)
cat("分层中位值的统计检验\n")
survd = survdiff(surv ~ Strata, df) #分层后的中位值统计检验，必须使用该格式，不能替换成fit_Strata
survd
survp = surv_pvalue(fit_Strata)
Stratakmp = survp$pval
n = survd$n[[dimnames(survd$n)[[1]][1]]] + survd$n[[dimnames(survd$n)[[1]][2]]]
Strata0 = dimnames(survd$n)[[1]][1]
Strata1 = dimnames(survd$n)[[1]][2]
samplesize_Strata0 = survd$n[[dimnames(survd$n)[[1]][1]]]
samplesize_Strata1 = survd$n[[dimnames(survd$n)[[1]][2]]]
obevents_strata0 = survd$obs[1]
obevents_strata1 = survd$obs[2]
expevents_strata0 = survd$exp[1]
expevents_strata1 = survd$exp[2]
table_fit_strata = summary(fit_Strata)$table
median_strata0 = table_fit_strata[, "median"][1]
median_strata1= table_fit_strata[, "median"][2]
lcl_strata0 = table_fit_strata[, "0.95LCL"][1]
ucl_strata0 = table_fit_strata[, "0.95UCL"][1]
lcl_strata1 = table_fit_strata[, "0.95LCL"][2]
ucl_strata1 = table_fit_strata[, "0.95UCL"][1]
chisq_strata = survd$chisq
pval_strata = survp$pval
if (is.na(Stratakmp) || is.na(median_strata0) || is.na(median_strata1)) {
  km_curves_differ = "无法进行统计推断，可能是事件结局数量太少。"
} else {
  if (Stratakmp < 0.05) {
    km_curves_differ = "经统计学推断，分层后两组中位时间【有】统计学差异！"
  } else {
    km_curves_differ = "经统计学推断，分层后两组中位时间【无】统计学差异。"
  }
}

cat(
    "\n",
    "生存函数相等性检验结果","\n",
    "---------------------------", "\n",
    "总样本量        ", as.character(n), "\n",
    "---------------------------", "\n",
    "分层因素        ", Strata0,"\n",
    "样本数量        ", samplesize_Strata0,"\n",
    "观测事件        ", as.character(round(obevents_strata0,3)), "\n",
    "预期事件        ", as.character(round(expevents_strata0,3)), "\n",
    "中位数值        ", as.character(round(median_strata0,3)), "\n",
    "可信下限        ", as.character(round(lcl_strata0,3)),"*", "\n",
    "可信上限        ", as.character(round(ucl_strata0,3)),"*", "\n",
    "---------------------------", "\n",
    "分层因素        ", Strata1,"\n",
    "样本数量        ", samplesize_Strata1,"\n",
    "观测事件        ", as.character(round(obevents_strata1,3)), "\n",
    "预期事件        ", as.character(round(expevents_strata1,3)), "\n",
    "中位数值        ", as.character(round(median_strata1,3)), "\n",
    "可信下限        ", as.character(round(lcl_strata1,3)),"*", "\n",
    "可信上限        ", as.character(round(ucl_strata1,3)),"*", "\n",
    "---------------------------","\n",
    "统计量","\n",
    "卡方数值        ", as.character(round(chisq_strata,3)), "\n",
    "检验概率        ", as.character(round(pval_strata,3)),"#", "\n",
    "---------------------------","\n",
    "*95%可信区间", "\n",
    "#Log-rank p value", "\n\n",
    km_curves_differ, "\n"
    )

km_curve_Strata = ggsurvplot(fit_Strata, df,
           pval = TRUE, #添加分层对比的p值
           pval.method = TRUE, #添加p值计算方法
           conf.int = FALSE, #可信区间是否显示
           risk.table = TRUE, #添加风险表
           risk.table.col = "strata", #根据分层更改风险表颜色
           linetype = "strata", # 根据分层更改线型
           surv.median.line = "hv", #同时显示垂直和水平参考线
           ggtheme = theme_bw(), #更改ggplot2的主题
           palette = c("#9E006E", "#2E9FDF"),#定义颜色
           title="Kaplan-Meier Curve by Strata", #定义图标题
           ylab="Survival(%)",#定义纵坐标名称
           xlab = " Time (months)", #定义横坐标名称
           )       
km_curve_Strata


cat("\nCox回归的单因素分析\n")
StrataCox = coxph(surv ~ Strata, df) #计算sex的Cox
StrataSum = summary(StrataCox) #总结一下
cat("分层单因素分析结果\n")
StrataCox

nevent = StrataCox$nevent
lrchi2 = StrataSum$logtest[1]
lrpval = StrataSum$logtest[3]
hr_cox = StrataSum$coefficients[2]
z_cox = StrataSum$coefficients[4]
p_cox = StrataSum$coefficients[5]
lcl_cox = StrataSum$conf.int[3]
ucl_cox = StrataSum$conf.int[4]
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
    "\n",
    "Cox回归检验结果","\n",
    "---------------------------", "\n",
    "样本数量        ", as.character(StrataCox$n), "\n",
    "事件数量        ", as.character(nevent), "\n",
    "---------------------------", "\n",
    "似然卡方        ", as.character(round(lrchi2,3)),"*1", "\n",
    "检验概率        ", as.character(round(lrpval,3)),"*2", "\n",
    "---------------------------", "\n",
    "风险比值        ", as.character(round(hr_cox,3)),"\n",
    "检验变量        ", as.character(round(z_cox,3)),"*3", "\n",
    "检验概率        ", as.character(round(p_cox,3)),"*4", "\n",
    "可信下限        ", as.character(round(lcl_cox,3)), "\n",
    "可信上限        ", as.character(round(ucl_cox,3)), "\n",
    "---------------------------","\n",
    "*1 最大似然估计卡方值", "\n",
    "*2 最大似然估计p值", "\n",
    "*3 z值", "\n",
    "*4 单因素Cox回归的p值", "\n\n",
    unicox, "\n"
    )

#多因素回归
#coxph(surv ~ factor1 + factor2 + ... + factorn, df) #单因素分析中p<0.05的df数据集中的变量，可用加号链接放在左侧

cat("\n分析结论\n")
cat(
    "汇总报告","\n",
    "---------------------------", "\n",
    "样本数量        ", as.character(n), "\n",
    "事件数量        ", as.character(nevent), "\n",
    "---------------------------", "\n",
    "---------------------------", "\n",
    "总体情况","\n",
    "中位数值        ", as.character(median_total), "\n",
    "可信下限        ", as.character(lcl95_total), "*1", "\n",
    "可信上限        ", as.character(ucl95_total), "*1", "\n",
    "---------------------------", "\n",
    "---------------------------", "\n",
    "分层分析","\n",
    "分层因素        ", Strata0,"\n",
    "样本数量        ", samplesize_Strata0,"\n",
    "中位数值        ", as.character(round(median_strata0,3)), "\n",
    "可信下限        ", as.character(round(lcl_strata0,3)),"*1", "\n",
    "可信上限        ", as.character(round(ucl_strata0,3)),"*1", "\n",
    "分层因素        ", Strata1,"\n",
    "样本数量        ", samplesize_Strata1,"\n",
    "中位数值        ", as.character(round(median_strata1,3)), "\n",
    "可信下限        ", as.character(round(lcl_strata1,3)),"*1", "\n",
    "可信上限        ", as.character(round(ucl_strata1,3)),"*1", "\n",
    "统计量","\n",
    "卡方数值        ", as.character(round(chisq_strata,3)), "\n",
    "检验概率        ", as.character(round(pval_strata,3)),"*2", "\n",
    "---------------------------","\n",
    "---------------------------","\n",
    "单因素Cox回归","\n",
    "风险比值        ", as.character(round(hr_cox,3)),"\n",
    "检验概率        ", as.character(round(p_cox,3)),"*3", "\n",
    "可信下限        ", as.character(round(lcl_cox,3)), "\n",
    "可信上限        ", as.character(round(ucl_cox,3)), "\n",
    "---------------------------","\n",
    "*1 95%可信区间", "\n",
    "*2 Log-rank p value", "\n",
    "*3 单因素Cox回归的p值", "\n\n"
    )

cat(
    "总结报告","\n",
    "1 如全人群推断总体的Kaplan-Meier曲线所示，本研究总体预后时间（Time_to_event_month）为：", as.character(median_total),"[", as.character(lcl95_total), " - ", as.character(ucl95_total), "] 月。","\n",
    "2 按照本研究Strata分层因素，两组预后时间（Time_to_event_month）分别为：", "\n",
    "    ", Strata0, ": ", as.character(round(median_strata0,3)), "[", as.character(round(lcl_strata0,3)), " - ", as.character(round(ucl_strata0,3)), "] 月;","\n",
    "    ", Strata1, ": ", as.character(round(median_strata1,3)), "[", as.character(round(lcl_strata1,3)), " - ", as.character(round(ucl_strata1,3)), "] 月;","\n",
    "    ", km_curves_differ, "\n",
    "3 对Strata分层因素进行单因素Cox回归时发现，", unicox, "\n"
    )

cat("\n", "感谢您使用本程序进行统计分析，再见。", "\n")
`;
