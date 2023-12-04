odds_output <- function(dir, sha, chapuy, schmitz, architecture, out){

  ##################################################################################################################################################
  # LOADING PATH
  ##################################################################################################################################################
  
  results_path <- paste0(dir, "/RESULTS"); 
  setwd(results_path); 
  
  ##################################################################################################################################################
  # DATA PREPARATION
  ##################################################################################################################################################
    
  if(grepl(paste("nonorm", collapse = NULL), dir) == TRUE){
    prefix <- "NN_"
  } else   if(grepl(paste("minmax", collapse = NULL), dir) == TRUE){
    prefix <- "N_"
  } else   if(grepl(paste("standard", collapse = NULL), dir) == TRUE){
    prefix <- "STD_"
  } 
  
  if(grepl(paste("OS12", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "OS12")
  } else if (grepl(paste("OS36", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "OS36")
  } else if (grepl(paste("OS60", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "OS60")
  } else if (grepl(paste("PFS12", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "PFS12")
  } else if (grepl(paste("PFS36", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "PFS36")
  } else if (grepl(paste("PFS60", collapse = NULL), dir) == TRUE){
    name_row <- paste0(prefix, architecture, "PFS60")
  }
  
  #Delete NA values of each clinical prognostic factor
  sha <- subset(sha, rev_ipi != "NA")
  chapuy <- subset(chapuy, rev_ipi != "NA")
  schmitz <- subset(schmitz, rev_ipi != "NA")

  # browser()
  
  #transform each prognostic factors with as-factor
  sha$rev_ipi <- as.factor(sha$rev_ipi); chapuy$rev_ipi <- as.factor(chapuy$rev_ipi); schmitz$rev_ipi <- as.factor(schmitz$rev_ipi)
  sha$coo <- as.factor(sha$coo); chapuy$coo <- as.factor(chapuy$coo); schmitz$coo <- as.factor(schmitz$coo)
  sha$Risk <- as.factor(sha$Risk); chapuy$Risk <- as.factor(chapuy$Risk); schmitz$Risk <- as.factor(schmitz$Risk)
  
  #Relvel Risk for output purposes
  table(sha$Risk); levels(sha$Risk) <- list("0" = "1", "1" = "0"); table(sha$Risk)
  table(chapuy$Risk); levels(chapuy$Risk) <- list("0" = "1", "1" = "0"); table(chapuy$Risk)
  table(schmitz$Risk); levels(schmitz$Risk) <- list("0" = "1", "1" = "0"); table(schmitz$Risk)
  
  if(grepl(paste("OS", collapse = NULL), dir) == TRUE){
    sha <- subset(sha, select=c("Risk", "rev_ipi","coo", "OS", "OS_time"))
    chapuy <- subset(chapuy, select=c("Risk", "rev_ipi","coo", "OS", "OS_time"))
    schmitz <- subset(schmitz, select=c("Risk", "rev_ipi","coo", "OS", "OS_time"))
    colnames(sha)[colnames(sha) == "OS_time"] <- "time_train"; sha$time_train <- as.numeric(sha$time_train)
    colnames(sha)[colnames(sha) == "OS"] <- "exitus_train"; sha$exitus_train <- as.numeric(sha$exitus_train)
    colnames(chapuy)[colnames(chapuy) == "OS_time"] <- "time_val"; chapuy$time_val <- as.numeric(chapuy$time_val)
    colnames(chapuy)[colnames(chapuy) == "OS"] <- "exitus_val"; chapuy$exitus_val <- as.numeric(chapuy$exitus_val)
    colnames(schmitz)[colnames(schmitz) == "OS_time"] <- "time_test"; schmitz$time_test <- as.numeric(schmitz$time_test)
    colnames(schmitz)[colnames(schmitz) == "OS"] <- "exitus_test"; schmitz$exitus_test <- as.numeric(schmitz$exitus_test)
  } else if (grepl(paste("PFS", collapse = NULL), dir) == TRUE){
    sha <- subset(sha, select=c("Risk", "rev_ipi", "coo","PFS", "PFS_time"))
    chapuy <- subset(chapuy, select=c("Risk", "rev_ipi","coo", "PFS", "PFS_time"))
    schmitz <- subset(schmitz, select=c("Risk", "rev_ipi","coo", "PFS", "PFS_time"))
    colnames(sha)[colnames(sha) == "PFS_time"] <- "time_train"; sha$time_train <- as.numeric(sha$time_train)
    colnames(sha)[colnames(sha) == "PFS"] <- "exitus_train"; sha$exitus_train <- as.numeric(sha$exitus_train)
    colnames(chapuy)[colnames(chapuy) == "PFS_time"] <- "time_val"; chapuy$time_val <- as.numeric(chapuy$time_val)
    colnames(chapuy)[colnames(chapuy) == "PFS"] <- "exitus_val"; chapuy$exitus_val <- as.numeric(chapuy$exitus_val)
    colnames(schmitz)[colnames(schmitz) == "PFS_time"] <- "time_test"; schmitz$time_test <- as.numeric(schmitz$time_test)
    colnames(schmitz)[colnames(schmitz) == "PFS"] <- "exitus_test"; schmitz$exitus_test <- as.numeric(schmitz$exitus_test)
  }
 
  ##################################################################################################################################################
  # LOGISTIC REGRESSION
  ##################################################################################################################################################

  # Select the right row
  index <- which((rownames(out) == name_row) == TRUE)
  
  # #sha
  # #rev_ipi
  # logistic_model <- glm(exitus_train ~ rev_ipi, data = sha, family = "binomial")
  # out[index,1] <- exp(summary(logistic_model)$coefficients[2])
  # out[index,2] <- summary(logistic_model)$coefficients[8]
  # #coo
  # logistic_model <- glm(exitus_train ~ coo, data = sha, family = "binomial")
  # out[index,3] <- exp(summary(logistic_model)$coefficients[2])
  # out[index,4] <- summary(logistic_model)$coefficients[8]
  # #AE Risk
  # logistic_model <- glm(exitus_train ~ Risk, data = sha, family = "binomial")
  # out[index,5] <- exp(summary(logistic_model)$coefficients[2])
  # out[index,6] <- summary(logistic_model)$coefficients[8]

  #chapuy
  #rev_ipi
  logistic_model <- glm(exitus_val ~ rev_ipi, data = chapuy, family = "binomial")
  out[index,1] <- exp(summary(logistic_model)$coefficients[2])
  out[index,2] <- summary(logistic_model)$coefficients[8]
  # #coo
  # logistic_model <- glm(exitus_val ~ coo, data = chapuy, family = "binomial")
  # out[index,3] <- exp(summary(logistic_model)$coefficients[2])
  # out[index,4] <- summary(logistic_model)$coefficients[8]
  #AE Risk
  logistic_model <- glm(exitus_val ~ Risk, data = chapuy, family = "binomial")
  out[index,3] <- exp(summary(logistic_model)$coefficients[2])
  out[index,4] <- summary(logistic_model)$coefficients[8]
  ## Multivariate
  logistic_model <- glm(exitus_val ~ rev_ipi+Risk, data = chapuy, family = "binomial")
  res <- plot_model(logistic_model)
  ggsave(paste(name_row,"chapuy_log_reg.pdf",sep="_"), res, width = 2000, height = 1000, units = "px")
  
  #schmitz
  #rev_ipi
  logistic_model <- glm(exitus_test ~ rev_ipi, data = schmitz, family = "binomial")
  out[index,5] <- exp(summary(logistic_model)$coefficients[2])
  out[index,6] <- summary(logistic_model)$coefficients[8]
  # #coo
  # logistic_model <- glm(exitus_test ~ coo, data = schmitz, family = "binomial")
  # out[index,9] <- exp(summary(logistic_model)$coefficients[2])
  # out[index,10] <- summary(logistic_model)$coefficients[8]
  #AE Risk
  logistic_model <- glm(exitus_test ~ Risk, data = schmitz, family = "binomial")
  out[index,7] <- exp(summary(logistic_model)$coefficients[2])
  out[index,8] <- summary(logistic_model)$coefficients[8]
  ## Multivariate
  logistic_model <- glm(exitus_test ~ rev_ipi+coo+Risk, data = schmitz, family = "binomial")
  res <- plot_model(logistic_model)
  ggsave(paste(name_row,"schmitz_log_reg.pdf",sep="_"), res, width = 2000, height = 1000, units = "px")
  
  setwd(dir)
return(out)
}