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













#PROBABILITY EQUIVALENT LEVEL ANALYSIS
#initializations
train <- dji30ret %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>% 
  filter(date < as.Date("2008-07-01")) %>%
  tail(750) 
test <- dji30ret %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  tail(150)
train_test=rbind(train,test)
cat("Dimensione train set:", nrow(train))
cat("Dimensione test set:", nrow(test))
cat("Dimensione dataset:", nrow(train_test))

marg_settings <- marginal_settings(
  train_size = 750, 
  refit_size = 50, #length of forecasting window of marginal models
  default_spec = default_garch_spec()
)

uncond_vine_settings <- vine_settings(
  train_size = 750,
  refit_size =25, #how many times use the same copula
)

cond_vine_settings <- vine_settings(
  train_size = 750, 
  refit_size = 25, #how many times use the same copula
  family_set = c("parametric"),
  vine_type = "dvine")

#definition of risk levels at which pelcov should be detected
risk_levels_v=c(0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07)
final_pel_uv=array(NA,dim=c(4,length(risk_levels_v)))
col_index <- setdiff(2:ncol(train_test), 8)
col_sampled <- c(2,3,4,5,6,7) #portfolio 1) 
#col_sampled <- c(9,10,11,12,13,14) #portfolio 6)
#col_sampled <- c(15,16,17,18,19,20) #portfolio 7)
asset_names=names(train_test)[col_sampled]
weights_portaf=setNames(c(rep(1/(length(asset_names)-2), length(asset_names) - 2), 0, 0), asset_names)
realized=rowSums(test[,col_sampled[1:4]]/4)
pel_plots <- vector("list", length(risk_levels_v)) #vector containing all plots to see graphically PELs

for (v in 1:length(risk_levels_v)){
  #unconditional rolling window estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
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
  
  #ONE CONDITIONAL ASSET
  #find the correct interval to search for pel
  pel_trial= seq(0.01, 0.91, by = 0.1)
  intersections_trial=array(0,dim=c(4,length(pel_trial)))
  pel_trial_str=as.character(pel_trial)
  df_cond_list_trial=vector("list", length = length(pel_trial))
  
  #one conditional asset rolling window estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
  for (i in 1:length(pel_trial)){
    cond_risk_roll <- estimate_risk_roll(
      data = train_test[col_sampled[1:5]],
      weights =weights_portaf[1:5],
      marginal_settings = marg_settings,
      vine_settings = cond_vine_settings,
      alpha = c(risk_levels_v[v]),
      risk_measures = c("VaR", "ES_mean"),
      n_samples = 500,
      cond_vars = c(names(train_test)[col_sampled[5]]),
      cond_u=pel_trial[i], 
      prior_resid_strategy = TRUE,
      trace = TRUE
    )
    df_cond_list_trial[[i]]=risk_estimates(cond_risk_roll,exceeded = TRUE)
  }
  
  #search intersections for pelcov e peloces
  dfs_var_trial<- list()
  dfs_es_trial <- list()
  for (i in 1:10) {
    dfs_var_trial[[i]] <- subset(df_cond_list_trial[[i]], risk_measure == "VaR" & cond_u == pel_trial_str[i])
    dfs_es_trial[[i]] <- subset(df_cond_list_trial[[i]], risk_measure == "ES_mean" & cond_u == pel_trial_str[i])
    diff_series1 <- dfs_var_trial[[i]][,2] - uncond_var
    diff_series2 <- dfs_es_trial[[i]][,2] - uncond_es
    for (j in 2:length(diff_series1)) {
      if (diff_series1[j - 1] * diff_series1[j] < 0) {
        intersections_trial[1,i] <- intersections_trial[1,i] + 1
      }
      if (diff_series2[j - 1] * diff_series2[j] < 0) {
        intersections_trial[2,i] <- intersections_trial[2,i] + 1
      }
    }
  }
  
  #take 10 equally spaced points in the interval where intersections with pelcov 
  #and pelcoes are !=0 (minimum index associated with the first non-zero value 
  #for both pelcov and pelcoes is taken so that the for loop only has to be 
  #executed once in subsequent rows and does not increase too much computational cost)
  min_intersect=min(head(which(intersections_trial[1,]!=0),1), head(which(intersections_trial[2,]!=0),1))
  max_intersect=max(tail(which(intersections_trial[1,]!=0),1), tail(which(intersections_trial[2,]!=0),1))
  pelcov_1d= round(seq(pel_trial[min_intersect], pel_trial[max_intersect], length = 10),2) #0.1, 0.5, by = 0.05***
  pelcov_1d_str=as.character(pelcov_1d)
  df_cond_list_1d=vector("list", length = length(pelcov_1d))
  intersections=array(0,dim=c(4,10)) #dim=c(4,9)***
  
  #one conditional asset rolling window estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
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
  }
  
  #search intersections for pelcov e peloces
  dfs_var_1d<- list()
  dfs_es_1d <- list()
  for (i in 1:10) { #1:9***
    dfs_var_1d[[i]] <- subset(df_cond_list_1d[[i]], risk_measure == "VaR" & cond_u == pelcov_1d_str[i])
    dfs_es_1d[[i]] <- subset(df_cond_list_1d[[i]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[i])
    diff_series1 <- dfs_var_1d[[i]][,2] - uncond_var
    diff_series2 <- dfs_es_1d[[i]][,2] - uncond_es
    for (j in 2:length(diff_series1)) {
      if (diff_series1[j - 1] * diff_series1[j] < 0) {
        intersections[1,i] <- intersections[1,i] + 1
      }
      if (diff_series2[j - 1] * diff_series2[j] < 0) {
        intersections[2,i] <- intersections[2,i] + 1
      }
      
    }
  }
  
  #stores graphical results
  #pelcov
  first_elements_var_1d <- sapply(dfs_var_1d, function(df) df[1, "risk_est"])
  order_index_var_1d <- order(first_elements_var_1d)
  colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#6a3d9a", "#a6cee3", "#b2df8a", "#fdbf6f", "#cab2d6", "#fb9a99", "gold") #remove gold***
  legend_data_var_1d <- data.frame(
    labels = c(paste0("Var_", pelcov_1d_str[order_index_var_1d])),
    colors = colors[order_index_var_1d]
  )
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "VaR"), 
              aes(x = row_num, y = risk_est), color="black", linewidth=0.75)+
    geom_line(data = subset(df_cond_list_1d[[1]], risk_measure == "VaR" & cond_u == pelcov_1d_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[1])))+
    geom_line(data = subset(df_cond_list_1d[[2]], risk_measure == "VaR" & cond_u == pelcov_1d_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[2])))+
    geom_line(data = subset(df_cond_list_1d[[3]], risk_measure == "VaR" & cond_u == pelcov_1d_str[3]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[3])))+
    geom_line(data = subset(df_cond_list_1d[[4]], risk_measure == "VaR" & cond_u == pelcov_1d_str[4]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[4])))+
    geom_line(data = subset(df_cond_list_1d[[5]], risk_measure == "VaR" & cond_u == pelcov_1d_str[5]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[5])))+
    geom_line(data = subset(df_cond_list_1d[[6]], risk_measure == "VaR" & cond_u == pelcov_1d_str[6]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[6])))+
    geom_line(data = subset(df_cond_list_1d[[7]], risk_measure == "VaR" & cond_u == pelcov_1d_str[7]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[7])))+
    geom_line(data = subset(df_cond_list_1d[[8]], risk_measure == "VaR" & cond_u == pelcov_1d_str[8]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[8])))+
    geom_line(data = subset(df_cond_list_1d[[9]], risk_measure == "VaR" & cond_u == pelcov_1d_str[9]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[9])))+
    geom_line(data = subset(df_cond_list_1d[[10]], risk_measure == "VaR" & cond_u == pelcov_1d_str[10]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_1d_str[10])))+ #comment these 2 lines***
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcov research graphically, 1 conditional asset",
         subtitle = paste0("Unconditional VaR in black, VaR conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_var_1d$colors, labels = legend_data_var_1d$labels)
  pel_plots[[v]][[1]] <- plot
  #pelcoes
  first_elements_es_1d <- sapply(dfs_es_1d, function(df) df[1, "risk_est"])
  order_index_es_1d <- order(first_elements_es_1d)
  legend_data_es_1d <- data.frame(
    labels = c(paste0("ES_mean_", pelcov_1d_str[order_index_es_1d])),
    colors = colors[order_index_es_1d]
  )
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "ES_mean"), 
              aes(x = row_num, y = risk_est), color="black", linewidth=0.75)+
    geom_line(data = subset(df_cond_list_1d[[1]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[1])))+
    geom_line(data = subset(df_cond_list_1d[[2]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[2])))+
    geom_line(data = subset(df_cond_list_1d[[3]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[3]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[3])))+
    geom_line(data = subset(df_cond_list_1d[[4]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[4]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[4])))+
    geom_line(data = subset(df_cond_list_1d[[5]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[5]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[5])))+
    geom_line(data = subset(df_cond_list_1d[[6]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[6]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[6])))+
    geom_line(data = subset(df_cond_list_1d[[7]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[7]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[7])))+
    geom_line(data = subset(df_cond_list_1d[[8]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[8]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[8])))+
    geom_line(data = subset(df_cond_list_1d[[9]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[9]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[9])))+
    geom_line(data = subset(df_cond_list_1d[[10]], risk_measure == "ES_mean" & cond_u == pelcov_1d_str[10]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_1d_str[10])))+ #comment these 2 lines***
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcoes research graphically, 1 conditional asset",
         subtitle = paste0("Unconditional ES in black, ES conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_es_1d$colors, labels = legend_data_es_1d$labels)
  pel_plots[[v]][[2]] <- plot
  
  
  
  #TWO CONDITIONAL ASSETS
  #find the correct interval to search for pel
  #two conditional assets rolling window estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
  for (i in 1:length(pel_trial)){
    cond_risk_roll <- estimate_risk_roll(
      data = train_test[col_sampled],
      weights =weights_portaf,
      marginal_settings = marg_settings,
      vine_settings = cond_vine_settings,
      alpha = c(risk_levels_v[v]),
      risk_measures = c("VaR", "ES_mean"),
      n_samples = 500,
      cond_vars = c(names(train_test)[col_sampled[5]],names(train_test)[col_sampled[6]]),
      cond_u=pel_trial[i], 
      prior_resid_strategy = TRUE,
      trace = TRUE
    )
    df_cond_list_trial[[i]]=risk_estimates(cond_risk_roll,exceeded = TRUE)
  }
  
  #search intersections for pelcov e pelcoes
  for (i in 1:10) {
    dfs_var_trial[[i]] <- subset(df_cond_list_trial[[i]], risk_measure == "VaR" & cond_u == pel_trial_str[i])
    dfs_es_trial[[i]] <- subset(df_cond_list_trial[[i]], risk_measure == "ES_mean" & cond_u == pel_trial_str[i])
    diff_series1 <- dfs_var_trial[[i]][,2] - uncond_var
    diff_series2 <- dfs_es_trial[[i]][,2] - uncond_es
    for (j in 2:length(diff_series1)) {
      if (diff_series1[j - 1] * diff_series1[j] < 0) {
        intersections_trial[3,i] <- intersections_trial[3,i] + 1
      }
      if (diff_series2[j - 1] * diff_series2[j] < 0) {
        intersections_trial[4,i] <- intersections_trial[4,i] + 1
      }
    }
  }
  
  #take 10 equally spaced points in the interval where intersections with pelcov 
  #and pelcoes are !=0 (minimum index associated with the first non-zero value 
  #for both pelcov and pelcoes is taken so that the for loop only has to be 
  #executed once in subsequent rows and does not increase too much computational cost)
  min_intersect=min(head(which(intersections_trial[3,]!=0),1), head(which(intersections_trial[4,]!=0),1))
  max_intersect=max(tail(which(intersections_trial[3,]!=0),1), tail(which(intersections_trial[4,]!=0),1))
  pelcov_2d_a= round(seq(pel_trial[min_intersect], pel_trial[max_intersect], length = 10),2) #0.1, 0.5, by = 0.05***
  pelcov_2d_a_str=as.character(pelcov_2d_a)
  df_cond_list_2d_a=vector("list", length = length(pelcov_2d_a))
  
  #two conditional assets rolling window estimation approach: fit ARMA-GARCH, take residuals, apply PIT, fit vine
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
  }
  
  #search intersections for pelcov e pelcoes
  dfs_var_2d_a<- list()
  dfs_es_2d_a <- list()
  for (i in 1:10) { #1:9***
    dfs_var_2d_a[[i]] <- subset(df_cond_list_2d_a[[i]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[i])
    dfs_es_2d_a[[i]] <- subset(df_cond_list_2d_a[[i]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[i])
    diff_series1 <- dfs_var_2d_a[[i]][,2] - uncond_var
    diff_series2 <- dfs_es_2d_a[[i]][,2] - uncond_es
    for (j in 2:length(diff_series1)) {
      if (diff_series1[j - 1] * diff_series1[j] < 0) {
        intersections[3,i] <- intersections[3,i] + 1
      }
      if (diff_series2[j - 1] * diff_series2[j] < 0) {
        intersections[4,i] <- intersections[4,i] + 1
      }
    }
  }
  
  #stores graphical results
  #pelcov
  first_elements_var_2d_a <- sapply(dfs_var_2d_a, function(df) df[1, "risk_est"])
  order_index_var_2d_a <- order(first_elements_var_2d_a)
  legend_data_var_2d_a <- data.frame(
    labels = c(paste0("Var_", pelcov_2d_a_str[order_index_var_2d_a])),
    colors = colors[order_index_var_2d_a]
  )
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "VaR"), 
              aes(x = row_num, y = risk_est), color="black", linewidth=0.75)+
    geom_line(data = subset(df_cond_list_2d_a[[1]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[1])))+
    geom_line(data = subset(df_cond_list_2d_a[[2]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[2])))+
    geom_line(data = subset(df_cond_list_2d_a[[3]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[3]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[3])))+
    geom_line(data = subset(df_cond_list_2d_a[[4]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[4]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[4])))+
    geom_line(data = subset(df_cond_list_2d_a[[5]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[5]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[5])))+
    geom_line(data = subset(df_cond_list_2d_a[[6]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[6]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[6])))+
    geom_line(data = subset(df_cond_list_2d_a[[7]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[7]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[7])))+
    geom_line(data = subset(df_cond_list_2d_a[[8]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[8]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[8])))+
    geom_line(data = subset(df_cond_list_2d_a[[9]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[9]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[9])))+
    geom_line(data = subset(df_cond_list_2d_a[[10]], risk_measure == "VaR" & cond_u == pelcov_2d_a_str[10]),
              aes(x = row_num, y = risk_est,col=paste0("Var_",pelcov_2d_a_str[10])))+ #comment these 2 lines***
    labs(x = "trading day",
         y = "portfolio VaR",
         col = "Risk measure",
         title = "Pelcov research graphically, 2 conditional assets with same value",
         subtitle = paste0("Unconditional VaR in black, VaR conf. level ",risk_levels_v[v]*100, "%"))+
    scale_color_manual(name = "Risk Measure", values = legend_data_var_2d_a$colors, labels = legend_data_var_2d_a$labels)
  pel_plots[[v]][[3]] <- plot
  #pelcoes
  first_elements_es_2d_a <- sapply(dfs_es_2d_a, function(df) df[1, "risk_est"])
  order_index_es_2d_a <- order(first_elements_es_2d_a)
  legend_data_es_2d_a <- data.frame(
    labels = c(paste0("ES_mean_", pelcov_2d_a_str[order_index_es_2d_a])),
    colors = colors[order_index_es_2d_a]
  )
  plot<- ggplot() +
    geom_line(data = subset(df_risk, risk_measure == "ES_mean"), 
              aes(x = row_num, y = risk_est), color="black", linewidth=0.75)+
    geom_line(data = subset(df_cond_list_2d_a[[1]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[1]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[1])))+
    geom_line(data = subset(df_cond_list_2d_a[[2]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[2]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[2])))+
    geom_line(data = subset(df_cond_list_2d_a[[3]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[3]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[3])))+
    geom_line(data = subset(df_cond_list_2d_a[[4]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[4]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[4])))+
    geom_line(data = subset(df_cond_list_2d_a[[5]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[5]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[5])))+
    geom_line(data = subset(df_cond_list_2d_a[[6]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[6]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[6])))+
    geom_line(data = subset(df_cond_list_2d_a[[7]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[7]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[7])))+
    geom_line(data = subset(df_cond_list_2d_a[[8]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[8]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[8])))+
    geom_line(data = subset(df_cond_list_2d_a[[9]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[9]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[9])))+
    geom_line(data = subset(df_cond_list_2d_a[[10]], risk_measure == "ES_mean" & cond_u == pelcov_2d_a_str[10]),
              aes(x = row_num, y = risk_est,col=paste0("ES_mean_",pelcov_2d_a_str[10])))+ #comment these 2 lines***
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

#display pelcov and pelcoes
#titles1 <- c("1D PELCoV", "2D PELCoV", "1D PELCoES", "2D PELCoES") #uncomment if you want more verbose graphs
titles2 <- lapply(1:4, function(i) {asset_names[1:4]})
titles3 <- lapply(1:4, function(i) {
  if (i%%2==1) {asset_names[5]}
  else {asset_names[5:6]}})
interp <- lapply(1:nrow(final_pel_uv), function(i) {
  #compute polynomial regression of degree 3
  df <- data.frame(risk_levels_v, final_pel_uv[i, ])
  fit <- lm(final_pel_uv[i, ] ~ poly(risk_levels_v, degree = 3), data = df)
  df$fit <- predict(fit, newdata = df)
  
  ggplot(df, aes(x = risk_levels_v, y = final_pel_uv[i, ])) +
    geom_point() +
    geom_line(aes(y = fit), color = "blue", linewidth= 0.8) + 
    #uncomment if you want more verbose graphs
    labs(#title = bquote(Omega ~ "={" * .(paste(titles2[[i]], collapse = ", ")) * "}, I={" * .(paste(titles3[[i]], collapse = ", ")) *"}"), #.(titles1[i]) * ": " *
         #subtitle = "Polynomial regression line of degree 3 in blue",
         x = "conf. level v",
         y = "PEL u_v")
})
interp
final_pel_uv
# 1) AA AXP BA BAC C CAT (last 2 are conditional assets)
#0.1315929 0.1660231 0.1778372 0.1778481 0.1649654 0.1872483 0.1906231 0.2068702 0.2038284 0.2165409 0.2177612 0.2236431 0.2227600
#0.1424242 0.1730055 0.1614634 0.1558475 0.1446237 0.1641935 0.1593125 0.1647778 0.1740670 0.1628826 0.1559259 0.1676025 0.1716779
#0.1991707 0.2256502 0.2386801 0.1976724 0.1855499 0.2674615 0.2794323 0.2229091 0.2261692 0.2394510 0.2390196 0.1909870 0.1939024
#0.1838813 0.1942647 0.2005169 0.1598492 0.1554198 0.2227960 0.2114919 0.1788485 0.1850932 0.1830472 0.1764336 0.1351351 0.1438542
# 6) DD DIS GE GM HD HPQ (last 2 are conditional assets)
#0.2625907 0.2513230 0.2600990 0.2710265 0.2724731 0.2582703 0.2813208 0.2835163 0.2676515 0.2819536 0.2853943 0.2952299 0.2984783
#0.2869123 0.2855908 0.2680608 0.2736190 0.2642260 0.2459919 0.2558029 0.2801931 0.2413611 0.2387772 0.2598615 0.2600721 0.2849837
#0.4172696 0.3310995 0.2238636 0.2810654 0.3400000 0.3877723 0.3779633 0.2991948 0.2538919 0.2716556 0.3130928 0.2772805 0.4365810
#0.3790141 0.3332824 0.2434351 0.2575058 0.3163871 0.3429487 0.3007859 0.2758795 0.2167840 0.2176060 0.2589209 0.2385388 0.3943970
# 7) IBM INTC JNJ JPM AIG KO (last 2 are conditional assets)
#0.1547107 0.1679472 0.1942857 0.1891977 0.1987059 0.2121538 0.2196500 0.2337767 0.2199415 0.2363700 0.2447399 0.2461094 0.2481570
#0.1363043 0.1498187 0.1664935 0.1620625 0.1524918 0.1765347 0.1767130 0.1716349 0.1822581 0.1672222 0.1815167 0.1832440 0.1900000
#0.1583371 0.2443167 0.2674737 0.2266436 0.2326786 0.2451887 0.2524026 0.2596215 0.2544361 0.2680000 0.2384028 0.2374672 0.2427226
#0.1284085 0.2119363 0.2175200 0.1891223 0.1790966 0.2027203 0.1973913 0.1940000 0.2088222 0.1983946 0.1701660 0.1741622 0.1762500
# 8) SPY 1329 ETFMIB GDAXIEX GCJ4 CCK4 (last 2 are conditional assets)
#0.4330632 0.4506940 0.4584478 0.4824329 0.4600778 0.4607028 0.4611606 0.4433006 0.4677973 0.4743519 0.4629703 0.4695897 0.4830257
#0.4662500 0.4890334 0.4480587 0.4923675 0.4370400 0.4294943 0.4650104 0.4129321 0.4423387 0.4410939 0.4499849 0.4560852 0.4409484
#0.4379221 0.4585175 0.4658549 0.4730480 0.4496754 0.4700928 0.4643601 0.4537870 0.4671429 0.4755140 0.4552888 0.4743959 0.4911050
#0.4688369 0.4788561 0.4554041 0.4971166 0.4409750 0.4489094 0.4571217 0.4215634 0.4504274 0.4589161 0.4456067 0.4564160 0.4462035




#save plots in pdfs
for (i in 1:length(risk_levels_v)) {
  pdf_name <- paste("risk_level_", risk_levels_v[i], ".pdf", sep = "")
  pdf(pdf_name)
  for (j in 1:length(pel_plots[[i]])) {
    print(pel_plots[[i]][[j]])
  }
  dev.off()
}



#all results below refers to the setting *** of the commented parts:
#change the code where you find *** in order to obtain same results.
# 1) AA AXP BA BAC C CAT (last 2 are conditional assets)
#final_pel_uv
#0.1771676 0.1828829 0.1597222 0.1973684 0.1896341 0.1829861 0.1892617 0.2022581 0.2091837 0.2235849 0.2416084 0.2266055 0.2330882
#0.1645455 0.1747253 0.1557692 0.1576923 0.1650000 0.1503106 0.1433735 0.1576531 0.1733831 0.1578125 0.1752874 0.1589744 0.1610465
#0.1787975 0.1993007 0.1891447 0.2049180 0.2178947 0.2099237 0.2117925 0.2223485 0.2355670 0.2446429 0.2437500 0.2417647 0.2461165
#0.1401840 0.1626582 0.1635870 0.1641566 0.1671642 0.1588816 0.1624031 0.1745161 0.1803922 0.1719697 0.1940299 0.1942308 0.1838028
# 6) DD DIS GE GM HD HPQ (last 2 are conditional assets)
#final_pel_uv
#0.2719585 0.2560976 0.2748188 0.2698606 0.2558712 0.2779116 0.2453917 0.2662996 0.2694656 0.2699134 0.3027108 0.2852941 0.2943182
#0.2688735 0.2742729 0.2592025 0.2698880 0.2458092 0.2678571 0.2431343 0.2339844 0.2529304 0.2468750 0.2751938 0.2613043 0.2619792
#0.2943478 0.2841564 0.2770000 0.2770202 0.2895722 0.2898734 0.2838028 0.2983240 0.3080537 0.3057325 0.3338415 0.3043333 0.3171429
#0.2624332 0.2893072 0.2823256 0.2746032 0.2718062 0.2834862 0.2798701 0.2730263 0.2859694 0.2704969 0.3158854 0.2881910 0.2959239
# 7) IBM INTC JNJ JPM AIG KO (last 2 are conditional assets)
#final_pel_uv
#0.2035714 0.2084746 0.1991489 0.2051351 0.1965839 0.2185484 0.2271429 0.2494220 0.2487805 0.2608974 0.2461765 0.2604938 0.2574830
#0.1821739 0.1916045 0.1632597 0.1842995 0.1596875 0.1761905 0.1783163 0.2109005 0.2110169 0.2142292 0.1902778 0.2046763 0.2154255
#0.2034483 0.2241259 0.2136076 0.2392638 0.2378261 0.2536765 0.2588816 0.2677632 0.2756000 0.2841667 0.2791262 0.3008621 0.3056818
#0.1708520 0.2002008 0.1769036 0.1829480 0.2000000 0.1964912 0.2046729 0.2206897 0.2170854 0.2339623 0.2133333 0.2370968 0.2378378
