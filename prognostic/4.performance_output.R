performance_output <- function(dir, sha, chapuy, schmitz, architecture, out){

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

  #transform each prognostic factors with as-factor
  sha$rev_ipi <- as.factor(sha$rev_ipi); chapuy$rev_ipi <- as.factor(chapuy$rev_ipi); schmitz$rev_ipi <- as.factor(schmitz$rev_ipi)
  sha$coo <- as.factor(sha$coo); chapuy$coo <- as.factor(chapuy$coo); schmitz$coo <- as.factor(schmitz$coo)
  sha$Risk <- as.factor(sha$Risk); chapuy$Risk <- as.factor(chapuy$Risk); schmitz$Risk <- as.factor(schmitz$Risk)

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
  # PERFORMANCE ASSESSMENT WITH MCC
  ##################################################################################################################################################
  
  # Select the right row
  index <- which((rownames(out) == name_row) == TRUE)
  
  #revipi
  logistic_model <- glm(exitus_train ~ rev_ipi, data = sha, family = "binomial")
  #sha
  predicted <- predict(logistic_model, sha, type = "response")
  
  # browser()
  
  ROCPred <- prediction(predicted, sha$exitus_train)
  mat <- performance(ROCPred, measure = "mat")
  # out[index,1] <- mat@y.values[[1]][2]
  #chapuy
  predicted <- predict(logistic_model, chapuy, type = "response")
  ROCPred <- prediction(predicted, chapuy$exitus_val)
  mat <- performance(ROCPred, measure = "mat")
  out[index,1] <- mat@y.values[[1]][2]
  #schmitz
  predicted <- predict(logistic_model, schmitz, type = "response")
  ROCPred <- prediction(predicted, schmitz$exitus_test)
  mat <- performance(ROCPred, measure = "mat")
  out[index,4] <- mat@y.values[[1]][2]
  
  #coo
  logistic_model <- glm(exitus_train ~ coo, data = sha, family = "binomial")
  #sha
  predicted <- predict(logistic_model, sha, type = "response")
  ROCPred <- prediction(predicted, sha$exitus_train)
  mat <- performance(ROCPred, measure = "mat")
  # out[index,2] <- mat@y.values[[1]][2]
  #chapuy
  predicted <- predict(logistic_model, chapuy, type = "response")
  ROCPred <- prediction(predicted, chapuy$exitus_val)
  mat <- performance(ROCPred, measure = "mat")
  out[index,2] <- mat@y.values[[1]][2]
  #schmitz
  predicted <- predict(logistic_model, schmitz, type = "response")
  ROCPred <- prediction(predicted, schmitz$exitus_test)
  mat <- performance(ROCPred, measure = "mat")
  out[index,5] <- mat@y.values[[1]][2]
  
  #AE Risk
  logistic_model <- glm(exitus_train ~ Risk, data = sha, family = "binomial")
  #sha
  predicted <- predict(logistic_model, sha, type = "response")
  ROCPred <- prediction(predicted, sha$exitus_train)
  mat <- performance(ROCPred, measure = "mat")
  # out[index,3] <- mat@y.values[[1]][2]
  #chapuy
  predicted <- predict(logistic_model, chapuy, type = "response")
  ROCPred <- prediction(predicted, chapuy$exitus_val)
  mat <- performance(ROCPred, measure = "mat")
  out[index,3] <- mat@y.values[[1]][2]
  #schmitz
  predicted <- predict(logistic_model, schmitz, type = "response")
  ROCPred <- prediction(predicted, schmitz$exitus_test)
  mat <- performance(ROCPred, measure = "mat")
  out[index,6] <- mat@y.values[[1]][2]
  
return(out)
  
}