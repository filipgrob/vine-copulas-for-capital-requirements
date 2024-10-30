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
data("dji30ret")















#RISK MEASURES COMPUTATION AND BACKTESTING
#plot djia constituents and index
data_inizio_crisi <- as.Date("2007-01-01")
dates=as.Date(row.names(dji30ret))
dji30ret_07_09 <- dji30ret[dates > data_inizio_crisi, ]

djia <- dji30ret_07_09 %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  gather(key = "Constituent", value = "log_return", -date)
ggplot(data = djia, aes(x = date, y = log_return, color = Constituent)) +
  geom_line() +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "DJIA Constituents Log Returns 2007-2009",
       subtitle = "LOESS smoothing line in black to check weak stationarity",
       x = "Year",
       y = "Log Return")

index <- data.frame(data = rowSums(dji30ret_07_09)) %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) 
ggplot(data = index, aes(x = date, y = data)) +
  geom_line() +
  labs(title = "DJIA Index Log returns 2007-2009",
       x = "Year",
       y = "Log Return")



#train/test set construction
train <- dji30ret %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>% 
  filter(date <= as.Date("2008-04-18")) %>% #2008-07-01***
  tail(1000) #1200***
#train=train[1:1000,1:length(train)] #***
test <- dji30ret %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  tail(200) #200***
train_test=rbind(train,test)
cat("Dimensione train set:", nrow(train))
cat("Dimensione test set:", nrow(test))
cat("Dimensione dataset:", nrow(train_test))


# #gridsearch for arma(p,q)-garch(x,y) orders selection (commented due to high computational cost)
# max_order <- 3 #possible values of p, q, x, y (from 1 to 3)
# combinaz <- expand.grid(p=1:max_order,q=1:max_order,x=1:max_order,y=1:max_order)
# aic_mat <- array(NA, dim=c(dim(dji30ret)[2],dim(combinaz)[1]))
# bic_mat <- array(NA, dim=c(dim(dji30ret)[2],dim(combinaz)[1]))
# best_aic=array(NA,dim=dim(dji30ret)[2])
# best_bic=array(NA,dim=dim(dji30ret)[2])
# best_combinaz=array(NA,dim=dim(dji30ret)[2])
# for (i in 2:dim(train)[2]){
#   print(i)
#   for (j in 1:dim(combinaz)[1]) {
#     model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(combinaz[j,3], combinaz[j,4])),
#                         mean.model = list(armaOrder = c(combinaz[j,1], combinaz[j,2])),
#                         distribution.model="sstd")
#     fit_model <- ugarchfit(spec = model, data = train[,i], solver = "hybrid")
#     aic_mat[i,j] <- infocriteria(fit_model)[1]
#     bic_mat[i,j] <- infocriteria(fit_model)[2]
#   }
# }
# 
# print(sum(is.na(aic_mat)))
# print(sum(is.na(bic_mat)))
# for (i in 1:dim(dji30ret)[2]){
#   best_aic[i]=which.min(aic_mat[i,])
#   best_bic[i]=which.min(bic_mat[i,])
#   if (best_aic[i]<best_bic[i]){
#     best_combinaz[i]=best_aic[i]
#   }
#   else {
#     best_combinaz[i]=best_bic[i]
#   }
# }
# print(best_aic)
# print(best_bic)
# print(best_combinaz)
# #ARMA-GARCH specification given the best orders
# spec_list <- lapply(1:30, function(i) {
#   aic_values <- combinaz[best_bic[i],]
#   spec <- default_garch_spec(ar = aic_values$p, ma = aic_values$q, arch = aic_values$x, garch = aic_values$y)
#   return(spec)
#   })
# names(spec_list)=names(dji30ret)


















#unconditional var,es estimation using a rolling window
#marginal settings
marg_settings <- marginal_settings(
  train_size = 1000,
  refit_size = 50, #length of forecasting window of marginal models
  #individual_spec = spec_list #uncomment this line if ARMA-GARCH gridsearch is done
  default_spec = default_garch_spec()
)

#vine settings
uncond_vine_settings <- vine_settings(
  train_size = 1000,
  refit_size = 50, #how many times use the same copula
)

col_index <- setdiff(2:ncol(train_test), 8)
col_sampled <- c(2,3,4,5,6,7)
#col_sampled <- c(15,16,17,18,19,20) #***
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
#roll window 1
#0.6380067 0.5463938 0.6029006 0.6344385 0.6718554
#0.8014358 0.9053650 0.1593253 0.3227831 0.5922897
#0.4676530 0.9880981 0.9515401 0.7557235 0.9220219
#0.5154959 0.9828692 0.4915308 0.7299531 0.8291133
#roll window 2
#0.7574107 0.4250166 0.5197523 0.5792903 0.5595016
#0.9649823 0.7189497 0.1273333 0.3255091 0.6412178
#0.3584866 0.9509725 0.7782875 0.6425933 0.8159217
#0.5889626 0.9750401 0.2606740 0.5787592 0.7132500
#roll window 3
#0.7303668 0.4845496 0.5032916 0.3530659 0.5161353
#0.9758730 0.6804913 0.2132273 0.3709773 0.6907675
#0.3565942 0.9610094 0.7203465 0.6127089 0.8014826
#0.6870309 0.9555553 0.4242265 0.6489335 0.8104965
#roll window 4
#0.4785026 0.8102419 0.6801353 0.6230811 0.5774870
#0.7419728 0.7221892 0.5306022 0.7773300 0.9391385
#0.4223057 0.8533465 0.9064178 0.7290525 0.7063390
#0.7907871 0.9558479 0.1554015 0.3116655 0.5337789
#Null hypothesis rejected 0 times

cat("Ljung-Box test pvalues on squared residuals\n")
ljung_mat_sqr_std
cat("Null hypothesis rejected",sum(ljung_mat_sqr_std < 0.05),"times")
#roll window 1
#0.3720645 0.3588128 0.5520004 0.7587870 0.8277244
#0.1763380 0.5868463 0.3012886 0.4141761 0.6076664
#0.3100602 0.7714612 0.9323075 0.9715605 0.9816345
#0.9376997 0.9005088 0.6384118 0.5807796 0.7340830
#roll window 2
#0.2901060 0.4123180 0.6495486 0.8579638 0.9231744
#0.1825294 0.6219297 0.3051888 0.4206254 0.6002521
#0.2898876 0.8078642 0.9686650 0.9793846 0.9912506
#0.8328646 0.8739790 0.7051593 0.6253267 0.7566121
#roll window 3
#0.2781621 0.3542505 0.4614083 0.7665455 0.9266283
#0.1737278 0.6032993 0.2955229 0.3834031 0.5551994
#0.3166856 0.7722314 0.8726943 0.9458312 0.9755234
#0.8350609 0.9535959 0.9202004 0.8726341 0.9193350
#roll window 4
#0.3840283 0.5761943 0.27725871 0.52219124 0.7443184
#0.1481619 0.6964143 0.06170925 0.09706642 0.1774252
#0.5043815 0.9213580 0.77727810 0.75624577 0.6348034
#0.8154017 0.7972667 0.81759214 0.87994837 0.8821706
#Null hypothesis rejected 0 times

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
  train_size = 1000,
  refit_size = 50, #how many times use the same copula
  family_set = c("parametric"),
  vine_type = "dvine")

#confidence level of the estimated quantile from the marginal distribution
pelcov=0.1
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

#PLOT OF ALL STRATEGIES
df_cond=risk_estimates(cond_risk_roll,exceeded = TRUE)
ggplot(df_cond) +
  geom_line(data=df_cond, aes(x = row_num, y = realized), col = "grey") +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == "prior_resid"), 
            aes(x = row_num, y = risk_est,col="Var_prior_resid")) +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("VaR_",pelcov_str))) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == "prior_resid"), 
            aes(x = row_num, y = risk_est,col="ES_mean_prior_resid")) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_str))) +
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Risk measure",
       title = "All conditional risk measures, conf. level 5%",
       subtitle = "Exceedances in red, portfolio realized log return in grey")+
  scale_fill_manual(
    name = "Risk Measure",
    labels = c(paste0("ES_mean_",pelcov_str),"ES_mean_prior_resid",paste0("VaR_",pelcov_str),"Var_prior_resid")
  )

#QUANTILE STRATEGY PLOT
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

#var, covar comparison
lowest_serie=df_cond$risk_est[df_cond$risk_measure=="VaR" & df_cond$cond_u==pelcov_str]
middle_serie=df_risk$risk_est[df_risk$risk_measure=="VaR"]
upper_serie=df_cond$risk_est[df_cond$risk_measure=="VaR" & df_cond$cond_u=="prior_resid"]
ggplot() +
  geom_line(data = subset(df_risk, risk_measure == "VaR"), 
            aes(x = row_num, y = risk_est,col="Var")) +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == "prior_resid"), 
            aes(x = row_num, y = risk_est,col="Var_prior_resid")) +
  geom_line(data = subset(df_cond, risk_measure == "VaR" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("VaR_",pelcov_str))) +
  labs(x = "trading day",
       y = "portfolio VaR",
       col = "Risk measure",
       title = "Comparison VaR & CoVaR, conf. level 5%")+
  scale_fill_manual(
    name = "Risk Measure",
    labels = c("VaR",paste0("VaR_",pelcov_str),"Var_prior_resid")
  )

#es, coes comparison
ggplot() +
  #geom_line(data=df_risk, aes(x = row_num, y = realized), col = "black") +
  #geom_line(data=df_cond, aes(x = row_num, y = realized), col = "yellow") +
  #geom_line(data=test, aes(x = (501:650), y =rowSums(test[,col_sampled[1:4]]/4)), col = "red") +
  geom_line(data = subset(df_risk, risk_measure == "ES_mean"), 
            aes(x = row_num, y = risk_est,col="ES_mean")) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == "prior_resid"), 
            aes(x = row_num, y = risk_est,col="ES_mean_prior_resid")) +
  geom_line(data = subset(df_cond, risk_measure == "ES_mean" & cond_u == pelcov_str), 
            aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_str))) +
  #geom_point(data=subset(df_cond, exceeded== TRUE),
  #           aes(x = row_num,y=realized),col="#db4f59")+
  #geom_point(data=subset(df_cond,exceeded== TRUE & risk_measure == "VaR" & cond_u == "prior_resid" & realized>middle_serie),
  #           aes(x = row_num,y=realized),col="green3")+
  #geom_point(data=subset(df_risk,exceeded== TRUE & risk_measure == "VaR" & realized>=lowest_serie),
  #           aes(x = row_num,y=realized),col="yellow")+
  labs(x = "trading day",
       y = "portfolio VaR",
       col = "Risk measure",
       title = "Comparison ES & CoES, conf. level 5%")+
  scale_fill_manual(
    name = "Risk Measure",
    labels = c("ES_mean",paste0("ES_mean_",pelcov_str),"ES_mean_prior_resid")
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
# $ expected.exceed: num 10
# $ actual.exceed  : num 24
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 15.1
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0.000103
# $ uc.Decision    : chr "Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 17.1
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0.000194
# $ cc.Decision    : chr "Reject H0"
cat("conditional quantile strategy for VaR estimation\n", str(bckt_var_quantile))
# $ expected.exceed: num 10
# $ actual.exceed  : num 7
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 1.05
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0.305
# $ uc.Decision    : chr "Fail to Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 1.56
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0.457
# $ cc.Decision    : chr "Fail to Reject H0"
cat("conditional prior residual strategy for VaR estimation\n", str(bckt_var_residual))
# $ expected.exceed: num 10
# $ actual.exceed  : num 50
# $ uc.H0          : chr "Correct Exceedances"
# $ uc.LRstat      : num 90
# $ uc.critical    : num 3.84
# $ uc.LRp         : num 0
# $ uc.Decision    : chr "Reject H0"
# $ cc.H0          : chr "Correct Exceedances & Independent"
# $ cc.LRstat      : num 110
# $ cc.critical    : num 5.99
# $ cc.LRp         : num 0
# $ cc.Decision    : chr "Reject H0"

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
# $ expected.exceed: num 10
# $ actual.exceed  : int 24
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.49
# $ p.value        : num 0.452
# $ Decision       : chr "Fail to Reject H0"
cat("conditional quantile strategy for ES estimation\n", str(bckt_es_quantile))
# $ expected.exceed: num 10
# $ actual.exceed  : int 7
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.156
# $ p.value        : num 0.0762
# $ Decision       : chr "Fail to Reject H0"
cat("conditional prior residual strategy for ES estimation\n", str(bckt_es_residual))
# $ expected.exceed: num 10
# $ actual.exceed  : int 50
# $ H1             : chr "Mean of Excess Violations of VaR is greater than zero"
# $ boot.p.value   : num 0.00208
# $ p.value        : num 0.000193
# $ Decision       : chr "Reject H0"



#***:
#these are the code-level changes that must be made to obtain the 
#results that refer to the specific example presented in section 4.4 
#(Omega_3={International Business Machines Corporation (IBM), 
#Intel Corporation (INTC), Johnson & Johnson (JNJ), JPMorgan Chase & Co. (JPM)} 
#each weighing 1/4, I_3={American International Group Inc. (AIG), 
#The Coca-Cola Company (KO)}, training set with 1000 observations 
#between 2003-09-24 and 2007-09-13, test set from 2008-04-21 to 2009-02-03)
