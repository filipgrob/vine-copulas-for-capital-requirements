#libraries
rm(list = ls())
set.seed(129)
library(VineCopula)
library(rugarch)
library(mvtsplot)
library(portvine)
library(stats)
library(rvinecopulib)
library(skewt)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggraph)















#RISK MEASURES COMPUTATION AND BACKTESTING
#import datasets
db <- read_excel("db.xlsx", sheet = "db2")
db=db[nrow(db):1,]
db <- column_to_rownames(db, var = "Date")

db_p <- read_excel("db.xlsx", sheet = "db2 price")
db_p=db_p[nrow(db_p):1,]
db_p <- column_to_rownames(db_p, var = "Date")
db_p_2 <- db_p %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  gather(key = "Constituent", value = "log_return", -date)

#plot univariate price series
ggplot(data = db_p_2, aes(x = date, y = log_return)) +
  geom_line() +
  labs(x = "Date", y = "Log Return") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 3) +
  ggtitle("Univariate Price Series of Portfolio Assets") +
  labs(subtitle = "GCJ4 and CCK4 are the conditional market indexes")

#compute correlation matrices
threshold_test <- as.Date("2007-09-4")
dates=as.Date(row.names(db_p))
realized_train=as.numeric(rowSums(subset(db_p, row.names(db_p)<= threshold_test)[,1:4])/4)
realized_test=as.numeric(rowSums(subset(db_p, row.names(db_p) > threshold_test)[,1:4])/4)
gold_train=subset(db_p, row.names(db_p)<= threshold_test)[,5]
gold_test=subset(db_p, row.names(db_p)> threshold_test)[,5]
cocoa_train=subset(db_p, row.names(db_p)<= threshold_test)[,6]
cocoa_test=subset(db_p, row.names(db_p)> threshold_test)[,6]

cor_spearman <- matrix(NA, nrow = 2, ncol = 2)
cor_kendall <- matrix(NA, nrow = 2, ncol = 2)
cor_spearman[1, 1] <- cor.test(realized_train, gold_train, method = "spearman",exact = FALSE)$estimate
cor_kendall[1, 1] <- cor.test(realized_train, gold_train, method = "kendall",exact = FALSE)$estimate
cor_spearman[1, 2] <- cor.test(realized_train, cocoa_train, method = "spearman",exact = FALSE)$estimate
cor_kendall[1, 2] <- cor.test(realized_train, cocoa_train, method = "kendall",exact = FALSE)$estimate
cor_spearman[2, 1] <- cor.test(realized_test, gold_test, method = "spearman",exact = FALSE)$estimate
cor_kendall[2, 1] <- cor.test(realized_test, gold_test, method = "kendall",exact = FALSE)$estimate
cor_spearman[2, 2] <- cor.test(realized_test, cocoa_test, method = "spearman",exact = FALSE)$estimate
cor_kendall[2, 2] <- cor.test(realized_test, cocoa_test, method = "kendall",exact = FALSE)$estimate
rownames(cor_spearman) <- c("realized_train", "realized_test")
colnames(cor_spearman) <- c("gold_train", "cocoa_train")
rownames(cor_kendall) <- c("realized_train", "realized_test")
colnames(cor_kendall) <- c("gold_train", "cocoa_train")
cor_spearman
#db2 price
#               gold_train cocoa_train
#realized_train  0.8845225   0.5328657
#realized_test  -0.1222775  -0.3032364
cor_kendall
#db2 price
#               gold_train cocoa_train
#realized_train  0.6932501   0.3852222
#realized_test  -0.1098750  -0.1762041

















#train/test set construction
train <- db %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>% 
  filter(date <= as.Date("2007-09-4"))
test <- db %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  filter(date > as.Date("2007-09-4"))
train_test=rbind(train,test)
cat("Dimensione train set:", nrow(train)) 
cat("Dimensione test set:", nrow(test))
cat("Dimensione dataset:", nrow(train_test)) 

#unconditional var,es estimation using a rolling window
#marginal settings
marg_settings <- marginal_settings(
  train_size =nrow(train), 
  refit_size = 50, #length of forecasting window of marginal models
  #individual_spec = spec_list #uncomment this line if ARMA-GARCH gridsearch is done
  default_spec = default_garch_spec()
)

#vine settings
uncond_vine_settings <- vine_settings(
  train_size = nrow(train),
  refit_size = 50, #how many times use the same copula
)

col_index <- setdiff(2:ncol(train_test), 8)
col_sampled <- c(2,3,4,5,6,7)
asset_names=names(train_test)[col_sampled]
weights_portaf=setNames(c(rep(1/(length(asset_names)-2), length(asset_names) - 2), 0, 0), asset_names)
realized=rowSums(test[,col_sampled[1:4]]/4)

#rolling estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
uncond_risk_roll <- estimate_risk_roll(
  data = train_test[,col_sampled[1:4]],
  weights=weights_portaf[1:4],
  marginal_settings = marg_settings,
  vine_settings = uncond_vine_settings,
  alpha = c(0.05),
  risk_measures = c("VaR", "ES_mean"),
  n_samples = 500,
  trace = TRUE
)
df_risk=risk_estimates(uncond_risk_roll,exceeded = TRUE)

#ljung-box test for serial autocorrelation at different lags (H0: no autocorrelation)
marginals <- fitted_marginals(uncond_risk_roll)
vines=fitted_vines(uncond_risk_roll)
lag=c(1,5,10,15,20)
ljung_mat_std <- array(NA, dim=c(4,5,4)) #asset x lag x rolling_window
ljung_mat_sqr_std <- array(NA, dim=c(4,5,4)) #asset x lag x rolling_window
for (i in 1:4){
  for (j in 1:4){ #portfolio with 4 asset
    std_resid <- roll_residuals(marginals[[names(train_test)[col_sampled[j]]]], roll_num = i)
    sqr_std_resid <- roll_residuals(marginals[[names(train_test)[col_sampled[j]]]], roll_num = i)**2
    for (t in 1:5){
      ljung_mat_std[j,t,i]=Box.test(std_resid, lag = lag[t], type = "Lju")$p.value
      ljung_mat_sqr_std[j,t,i]=Box.test(sqr_std_resid, lag = lag[t], type = "Lju")$p.value
    }
  }
}
cat("Ljung-Box test p-values on residuals\n")
ljung_mat_std
cat("Null hypothesis rejected",sum(ljung_mat_std < 0.05),"times")
#db2: Null hypothesis rejected 11 times

cat("Ljung-Box test pvalues on squared residuals\n")
ljung_mat_sqr_std
cat("Null hypothesis rejected",sum(ljung_mat_sqr_std < 0.05),"times")
#db2: Null hypothesis rejected 7 times

#plot of results
x_geom_point=df_risk$row_num[df_risk$exceeded]
y_geom_point=df_risk$realized[df_risk$exceeded]
df_geom_point=data.frame(df_risk$realized[df_risk$exceeded])
df_risk %>% 
  ggplot() +
  geom_line(aes(x = row_num, y = realized), col = "grey") +
  geom_line(aes(x = row_num, y = risk_est, col = factor(risk_measure))) +
  scale_fill_manual() +
  geom_point(aes(x = x_geom_point, y = y_geom_point), 
             data = df_geom_point, 
             col = "#db4f59")+
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Risk measure",
       title = "Unconditional risk measures, conf. level 5%",
       subtitle = "Exceedances in red, portfolio log return in grey")

#conditional var,es estimation using a rolling window
#vine settings
cond_vine_settings <- vine_settings(
  train_size = nrow(train),
  refit_size = 50, #how many times use the same copula
  family_set = c("parametric"),
  vine_type = "dvine")

#confidence level of the estimated quantile from the marginal distribution
pelcov=0.9
pelcov_str=as.character(pelcov)
#rolling estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
cond_risk_roll <- estimate_risk_roll(
  data = train_test[col_sampled[1:6]],
  weights =weights_portaf ,
  marginal_settings = marg_settings,
  vine_settings = cond_vine_settings,
  alpha = c(0.05),
  risk_measures = c("VaR", "ES_mean"),
  n_samples = 500,
  cond_vars = c(names(train_test)[col_sampled[5]],names(train_test)[col_sampled[6]]),
  cond_u=pelcov, #value inside pelcov is associated to both conditioning variables 
  ##QUANTILE STRATEGY=marginal market index on copula scale
  #(marginals for copulas are unif-->quantile is the confidence level itself),
  prior_resid_strategy = TRUE,
  #RESIDUAL STRATEGY=conditioning values for the forecast at time t
  #are the PIT of index I residual at time t-1
  trace = TRUE
)

#QUANTILE STRATEGY PLOT
df_cond=risk_estimates(cond_risk_roll,exceeded = TRUE)
ggplot(df_cond) +
  geom_line(aes(x = row_num, y = realized), col = "grey") +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("VaR_",pelcov_str))) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_str))) +
  geom_line(data = subset(df_cond, cond_u == pelcov_str), 
            aes(x = row_num, y =!!as.name(names(train_test)[col_sampled[5]]),
                col=paste0(names(train_test)[col_sampled[5]],"_",pelcov))) +
  geom_line(data = subset(df_cond, cond_u == pelcov_str), 
            aes(x = row_num, y =!!as.name(names(train_test)[col_sampled[6]]),
                col=paste0(names(train_test)[col_sampled[6]],"_",pelcov))) +
  geom_point(data=subset(df_cond,exceeded== TRUE & cond_u == pelcov_str),
             aes(x = row_num, y=realized),col="#db4f59")+
  scale_fill_manual(
    name = "Risk Measure",
    labels = c(paste0(names(train_test)[col_sampled[6]],"_",pelcov),
               paste0(names(train_test)[col_sampled[5]],"_",pelcov),
               paste0("ES_mean_",pelcov_str),paste0("VaR_",pelcov_str)))+
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Risk measure",
       title = "Quantile based conditional risk measures, conf. level 5%",
       subtitle = "Exceedances in red, portfolio realized log return in grey")

#PRIOR RESIDUAL COPULA SCALE STRATEGY PLOT
#It should be noted that the conditional series based on the fitted 
#residuals of the time unit before will most likely exaggerate 
#sudden high volatility situations
ggplot(df_cond, aes(x = row_num)) +
  geom_line(aes(y = realized), col = "grey") +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == "prior_resid"), aes(y = risk_est,col="VaR_prior_resid")) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == "prior_resid"), aes(y = risk_est,col="ES_mean_prior_resid")) +
  geom_point(data=subset(df_cond,exceeded== TRUE & cond_u == "prior_resid"),aes(y=realized),col="#db4f59")+
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Risk measure",
       title = "Prior Residual conditional risk measures, conf. level 5%",
       subtitle = "Exceedances in red, portfolio realized log return in grey")+
  scale_fill_manual(
    name = "Risk Measure",
    labels = c("ES_mean_prior_resid","VaR_prior_resid")
  )

#backtesting
uncond_var=df_risk[df_risk$risk_measure=="VaR",2]
cond_var_quantile=df_cond[df_cond$risk_measure=="VaR" & df_cond$cond_u==pelcov_str ,2]
cond_var_residual=df_cond[df_cond$risk_measure=="VaR" & df_cond$cond_u=="prior_resid",2]
bckt_var_uncond=VaRTest(alpha = 0.05, actual=realized, 
                        VaR=uncond_var, conf.level = 0.95) #0.05)#
bckt_var_quantile=VaRTest(alpha = 0.05, actual=realized, 
                          VaR=cond_var_quantile, conf.level = 0.95)
bckt_var_residual=VaRTest(alpha = 0.05, actual=realized, 
                          VaR=cond_var_residual, conf.level = 0.95)
cat("unconditional strategy for VaR estimation\n", str(bckt_var_uncond))
# db2:
# $ expected.exceed: num 16
# $ actual.exceed  : num 15
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 0.0958
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0.757
# $ uc.Decision    : chr "Fail to Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 0.225
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0.893
# $ cc.Decision    : chr "Fail to Reject H0"
cat("conditional quantile strategy for VaR estimation\n", str(bckt_var_quantile))
# db2 0.1: 
# $ expected.exceed: num 16
# $ actual.exceed  : num 11
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 1.97
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0.16
# $ uc.Decision    : chr "Fail to Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 2.75
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0.253
# $ cc.Decision    : chr "Fail to Reject H0"

# db2 0.9: 
# $ expected.exceed: num 16
# $ actual.exceed  : num 24
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 3.47
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0.0627
# $ uc.Decision    : chr "Fail to Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 3.5
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0.174
# $ cc.Decision    : chr "Fail to Reject H0"
cat("conditional prior residual strategy for VaR estimation\n", str(bckt_var_residual))

uncond_es=df_risk[df_risk$risk_measure=="ES_mean",2]
cond_es_quantile=df_cond[df_cond$risk_measure=="ES_mean" & df_cond$cond_u==pelcov_str ,2]
cond_es_residual=df_cond[df_cond$risk_measure=="ES_mean" & df_cond$cond_u=="prior_resid",2]
bckt_es_uncond=ESTest(alpha = 0.05, actual=realized, ES=uncond_es, 
                      VaR=uncond_var, conf.level = 0.95,boot = TRUE, n.boot = 1000)
bckt_es_quantile=ESTest(alpha = 0.05, realized, cond_es_quantile, 
                        cond_var_quantile, conf.level = 0.95, boot = TRUE, n.boot = 1000)
bckt_es_residual=ESTest(alpha = 0.05, realized,cond_es_residual, 
                        cond_var_residual, conf.level = 0.95, boot = TRUE, n.boot = 1000)
cat("unconditional strategy for ES estimation\n", str(bckt_es_uncond))
# db2: 
# $ expected.exceed: num 16
# $ actual.exceed  : int 15
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.309
# $ p.value        : num 0.227
# $ Decision       : chr "Fail to Reject H0"

cat("conditional quantile strategy for ES estimation\n", str(bckt_es_quantile))
# db2 0.1: 
# $ expected.exceed: num 16
# $ actual.exceed  : int 11
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.275
# $ p.value        : num 0.182
# $ Decision       : chr "Fail to Reject H0"

# db2 0.9: 
# $ expected.exceed: num 16
# $ actual.exceed  : int 24
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.144
# $ p.value        : num 0.0825
# $ Decision       : chr "Fail to Reject H0"
cat("conditional prior residual strategy for ES estimation\n", str(bckt_es_residual))



#plot vine structures
uncond_vines=fitted_vines(uncond_risk_roll)
cond_vines=fitted_vines(cond_risk_roll)
plot(uncond_vines[[1]], tree="ALL", var_names="legend", edge_labels="tau")
plot(cond_vines[[1]], tree=1:5, var_names="legend", edge_labels="tau")



#PROBABILITY EQUIVALENT LEVEL ANALYSIS
#definition of risk levels at which pelcov should be detected
risk_levels_v=c(0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07)
final_pel_uv=array(NA,dim=c(4,length(risk_levels_v)))
pel_plots <- vector("list", length(risk_levels_v)) #vector containing all plots to see graphically PELs

for (v in 1:length(risk_levels_v)){
  #rolling estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
  uncond_risk_roll <- estimate_risk_roll(
    data = train_test[,col_sampled[1:4]],
    weights=weights_portaf[1:4],
    marginal_settings = marg_settings,
    vine_settings = uncond_vine_settings,
    alpha = c(risk_levels_v[v]),
    risk_measures = c("VaR", "ES_mean"),
    n_samples = 500,
    trace = TRUE
  )
  df_risk=risk_estimates(uncond_risk_roll,exceeded = TRUE)
  uncond_var=df_risk[df_risk$risk_measure=="VaR",2]
  uncond_es=df_risk[df_risk$risk_measure=="ES_mean",2]
  intersections=array(0,dim=c(4,9))
  
  
  
  #1D CASE
  pelcov_1d= c(0.1,0.9)
  pelcov_1d_str=as.character(pelcov_1d)
  df_cond_list_1d=vector("list", length = length(pelcov_1d))
  var_tests_1d=vector("list", length = length(pelcov_1d))
  es_tests_1d=vector("list", length = length(pelcov_1d))
  
  #rolling estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
  for (i in 1:length(pelcov_1d)){
    cond_risk_roll <- estimate_risk_roll(
      data = train_test[col_sampled[1:5]],
      weights =weights_portaf[1:5],
      marginal_settings = marg_settings,
      vine_settings = cond_vine_settings,
      alpha = c(risk_levels_v[v]),
      risk_measures = c("VaR", "ES_mean"),
      n_samples = 500,
      cond_vars = c(names(train_test)[col_sampled[5]]),
      cond_u=pelcov_1d[i], 
      prior_resid_strategy = TRUE,
      trace = TRUE
    )
    df_cond_list_1d[[i]]=risk_estimates(cond_risk_roll,exceeded = TRUE)
    cond_var_quantile=df_cond_list_1d[[i]][df_cond_list_1d[[i]]$risk_measure=="VaR" & df_cond_list_1d[[i]]$cond_u==pelcov_1d_str[i] ,2]
    var_tests_1d[[i]]=VaRTest(alpha = risk_levels_v[v], actual=realized, 
                              VaR=cond_var_quantile, conf.level = 1-risk_levels_v[v])
    cond_es_quantile=df_cond_list_1d[[i]][df_cond_list_1d[[i]]$risk_measure=="ES_mean" & df_cond_list_1d[[i]]$cond_u==pelcov_1d_str[i] ,2]
    es_tests_1d[[i]]=ESTest(alpha = risk_levels_v[v], realized, cond_es_quantile, 
                            cond_var_quantile, conf.level = 1-risk_levels_v[v], boot = TRUE, n.boot = 1000)
  }
  
  #pelcov
  dfs_var_1d<- list()
  for (i in 1:2) {
    dfs_var_1d[[i]] <- subset(df_cond_list_1d[[i]], risk_measure == "VaR" & cond_u == pelcov_1d_str[i])
    diff_series <- dfs_var_1d[[i]][,2] - uncond_var
    for (j in 2:length(diff_series)) {
      if (diff_series[j - 1] * diff_series[j] < 0) {
        intersections[1,i] <- intersections[1,i] + 1
      }
    }
  }
  
  first_elements_var_1d <- sapply(dfs_var_1d, function(df) df[1, "risk_est"])
  order_index_var_1d <- order(first_elements_var_1d)
  colors <- c("#33a02c", "#fb9a99")
  legend_data_var_1d <- data.frame(
    labels = c(paste0("Var_", pelcov_1d_str[order_index_var_1d])),
    colors = colors[order_index_var_1d]
  )
  
  #stores graphical results
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "VaR"), 
              aes(x = row_num, y = risk_est), color="black")+
    geom_line(data = subset(df_cond_list_1d[[1]], risk_measure == "VaR" & cond_u == pelcov_1d_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[1])))+
    geom_line(data = subset(df_cond_list_1d[[2]], risk_measure == "VaR" & cond_u == pelcov_1d_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[2])))+
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcov research graphically, 1 conditional asset",
         subtitle = paste0("Unconditional VaR in black, VaR conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_var_1d$colors, labels = legend_data_var_1d$labels)
  pel_plots[[v]][[1]] <- plot
  
  #pelcoes
  dfs_es_1d <- list()
  for (i in 1:2) {
    dfs_es_1d[[i]] <- subset(df_cond_list_1d[[i]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[i])
    diff_series <- dfs_es_1d[[i]][,2] - uncond_es
    for (j in 2:length(diff_series)) {
      if (diff_series[j - 1] * diff_series[j] < 0) {
        intersections[2,i] <- intersections[2,i] + 1
      }
    }
  }
  
  first_elements_es_1d <- sapply(dfs_es_1d, function(df) df[1, "risk_est"])
  order_index_es_1d <- order(first_elements_es_1d)
  legend_data_es_1d <- data.frame(
    labels = c(paste0("ES_mean_", pelcov_1d_str[order_index_es_1d])),
    colors = colors[order_index_es_1d]
  )
  
  #stores graphical results
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "ES_mean"), 
              aes(x = row_num, y = risk_est), color="black")+
    geom_line(data = subset(df_cond_list_1d[[1]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[1])))+
    geom_line(data = subset(df_cond_list_1d[[2]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[2])))+
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcoes research graphically, 1 conditional asset",
         subtitle = paste0("Unconditional ES in black, ES conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_es_1d$colors, labels = legend_data_es_1d$labels)
  pel_plots[[v]][[2]] <- plot
  
  
  
  #2D CASE
  pelcov_2d_a= c(0.1,0.9)
  pelcov_2d_a_str=as.character(pelcov_2d_a)
  df_cond_list_2d_a=vector("list", length = length(pelcov_2d_a))
  var_tests_2d_a=vector("list", length = length(pelcov_2d_a))
  es_tests_2d_a=vector("list", length = length(pelcov_2d_a))
  
  #rolling estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
  for (i in 1:length(pelcov_2d_a)){
    cond_risk_roll <- estimate_risk_roll(
      data = train_test[col_sampled],
      weights =weights_portaf,
      marginal_settings = marg_settings,
      vine_settings = cond_vine_settings,
      alpha = c(risk_levels_v[v]),
      risk_measures = c("VaR", "ES_mean"),
      n_samples = 500,
      cond_vars = c(names(train_test)[col_sampled[5]],names(train_test)[col_sampled[6]]),
      cond_u=pelcov_2d_a[i],
      prior_resid_strategy = TRUE,
      trace = TRUE
    )
    df_cond_list_2d_a[[i]]=risk_estimates(cond_risk_roll,exceeded = TRUE)
    cond_var_quantile=df_cond_list_2d_a[[i]][df_cond_list_2d_a[[i]]$risk_measure=="VaR" & df_cond_list_2d_a[[i]]$cond_u==pelcov_2d_a_str[i] ,2]
    var_tests_2d_a[[i]]=VaRTest(alpha = risk_levels_v[v], actual=realized, 
                                VaR=cond_var_quantile, conf.level = 1-risk_levels_v[v])
    cond_es_quantile=df_cond_list_2d_a[[i]][df_cond_list_2d_a[[i]]$risk_measure=="ES_mean" & df_cond_list_2d_a[[i]]$cond_u==pelcov_2d_a_str[i] ,2]
    es_tests_2d_a[[i]]=ESTest(alpha = risk_levels_v[v], realized, cond_es_quantile, 
                              cond_var_quantile, conf.level = 1-risk_levels_v[v], boot = TRUE, n.boot = 1000)
  }
  
  #pelcov
  dfs_var_2d_a<- list()
  for (i in 1:2) {
    dfs_var_2d_a[[i]] <- subset(df_cond_list_2d_a[[i]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[i])
    diff_series <- dfs_var_2d_a[[i]][,2] - uncond_var
    for (j in 2:length(diff_series)) {
      if (diff_series[j - 1] * diff_series[j] < 0) {
        intersections[3,i] <- intersections[3,i] + 1
      }
    }
  }
  
  first_elements_var_2d_a <- sapply(dfs_var_2d_a, function(df) df[1, "risk_est"])
  order_index_var_2d_a <- order(first_elements_var_2d_a)
  legend_data_var_2d_a <- data.frame(
    labels = c(paste0("Var_", pelcov_2d_a_str[order_index_var_2d_a])),
    colors = colors[order_index_var_2d_a]
  )
  
  #stores graphical results
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "VaR"), 
              aes(x = row_num, y = risk_est), color="black")+
    geom_line(data = subset(df_cond_list_2d_a[[1]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[1])))+
    geom_line(data = subset(df_cond_list_2d_a[[2]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[2])))+
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcov research graphically, 2 conditional assets with same value",
         subtitle = paste0("Unconditional VaR in black, VaR conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_var_2d_a$colors, labels = legend_data_var_2d_a$labels)
  pel_plots[[v]][[3]] <- plot
  
  #pelcoes
  dfs_es_2d_a <- list()
  for (i in 1:2) {
    dfs_es_2d_a[[i]] <- subset(df_cond_list_2d_a[[i]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[i])
    diff_series <- dfs_es_2d_a[[i]][,2] - uncond_es
    for (j in 2:length(diff_series)) {
      if (diff_series[j - 1] * diff_series[j] < 0) {
        intersections[4,i] <- intersections[4,i] + 1
      }
    }
  }
  
  first_elements_es_2d_a <- sapply(dfs_es_2d_a, function(df) df[1, "risk_est"])
  order_index_es_2d_a <- order(first_elements_es_2d_a)
  legend_data_es_2d_a <- data.frame(
    labels = c(paste0("ES_mean_", pelcov_2d_a_str[order_index_es_2d_a])),
    colors = colors[order_index_es_2d_a]
  )
  
  #stores graphical results
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "ES_mean"), 
              aes(x = row_num, y = risk_est), color="black")+
    geom_line(data = subset(df_cond_list_2d_a[[1]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[1])))+
    geom_line(data = subset(df_cond_list_2d_a[[2]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[2])))+
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcoes research graphically, 2 conditional assets with same value",
         subtitle = paste0("Unconditional ES in black, ES conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_es_2d_a$colors, labels = legend_data_es_2d_a$labels)
  pel_plots[[v]][[4]] <- plot
  
  
  
  #pel computation at confidence level v
  for (u in 1:4){
    final_pel_uv[u,v]=sum(pelcov_1d*intersections[u,])/sum(intersections[u,])
  }
}

#save plots in pdfs
for (i in 1:length(risk_levels_v)) {
  pdf_name <- paste("risk_level_", risk_levels_v[i], ".pdf", sep = "")
  pdf(pdf_name)
  for (j in 1:length(pel_plots[[i]])) {
    print(pel_plots[[i]][[j]])
  }
  dev.off()
}
