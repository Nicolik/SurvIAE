gene_sig <- function(dir, sha, chapuy, schmitz, architecture, best_genes){

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

  #transform each prognostic factors with as-factor
  sha$rev_ipi <- as.factor(sha$rev_ipi); chapuy$rev_ipi <- as.factor(chapuy$rev_ipi); schmitz$rev_ipi <- as.factor(schmitz$rev_ipi)
  sha$coo <- as.factor(sha$coo); chapuy$coo <- as.factor(chapuy$coo); schmitz$coo <- as.factor(schmitz$coo)
  sha$Risk <- as.factor(sha$Risk); chapuy$Risk <- as.factor(chapuy$Risk); schmitz$Risk <- as.factor(schmitz$Risk)
  
  #Gene name standardization
  best_genes[grepl(paste("-", collapse = NULL),best_genes) == TRUE] <- chartr("-", ".", best_genes[grepl(paste("-", collapse = NULL),best_genes) == TRUE])
  
  if(grepl(paste("OS", collapse = NULL), dir) == TRUE){
    sha <- subset(sha, select=c(best_genes, "rev_ipi","coo", "OS", "OS_time"))
    chapuy <- subset(chapuy, select=c(best_genes, "rev_ipi","coo",  "OS", "OS_time"))
    schmitz <- subset(schmitz, select=c(best_genes, "rev_ipi","coo", "OS", "OS_time"))
    colnames(sha)[colnames(sha) == "OS_time"] <- "time_train"; sha$time_train <- as.numeric(sha$time_train)
    colnames(sha)[colnames(sha) == "OS"] <- "exitus_train"; sha$exitus_train <- as.numeric(sha$exitus_train)
    colnames(chapuy)[colnames(chapuy) == "OS_time"] <- "time_val"; chapuy$time_val <- as.numeric(chapuy$time_val)
    colnames(chapuy)[colnames(chapuy) == "OS"] <- "exitus_val"; chapuy$exitus_val <- as.numeric(chapuy$exitus_val)
    colnames(schmitz)[colnames(schmitz) == "OS_time"] <- "time_test"; schmitz$time_test <- as.numeric(schmitz$time_test)
    colnames(schmitz)[colnames(schmitz) == "OS"] <- "exitus_test"; schmitz$exitus_test <- as.numeric(schmitz$exitus_test)
  } else if (grepl(paste("PFS", collapse = NULL), dir) == TRUE){
    sha <- subset(sha, select=c(best_genes, "rev_ipi","coo",  "PFS", "PFS_time"))
    chapuy <- subset(chapuy, select=c(best_genes,"rev_ipi","coo", "PFS", "PFS_time"))
    schmitz <- subset(schmitz, select=c(best_genes, "rev_ipi","coo", "PFS", "PFS_time"))
    colnames(sha)[colnames(sha) == "PFS_time"] <- "time_train"; sha$time_train <- as.numeric(sha$time_train)
    colnames(sha)[colnames(sha) == "PFS"] <- "exitus_train"; sha$exitus_train <- as.numeric(sha$exitus_train)
    colnames(chapuy)[colnames(chapuy) == "PFS_time"] <- "time_val"; chapuy$time_val <- as.numeric(chapuy$time_val)
    colnames(chapuy)[colnames(chapuy) == "PFS"] <- "exitus_val"; chapuy$exitus_val <- as.numeric(chapuy$exitus_val)
    colnames(schmitz)[colnames(schmitz) == "PFS_time"] <- "time_test"; schmitz$time_test <- as.numeric(schmitz$time_test)
    colnames(schmitz)[colnames(schmitz) == "PFS"] <- "exitus_test"; schmitz$exitus_test <- as.numeric(schmitz$exitus_test)
  }

  ##################################################################################################################################################
  # GENE CATHEROGIZATION ACCORDING TO THE OUTCOME
  ##################################################################################################################################################
  
  ## GENES
  
  output <- matrix(0, 10, 6)
  row.names(output) <- best_genes
  colnames(output) <- c("Sha_ODDs", "Sha_pval", "Chapuy_ODDs", "Chapuy_pval", "Schmits_ODDs", "Schmits_pval")
  
  for (i in best_genes) {
    
    print(i)
    
    res.cut2 <- surv_cutpoint(sha, time = "time_train", event = "exitus_train", variables = i);
    cutoff <- res.cut2$cutpoint[1,1]
    
    sha$X <- 0; sha$X[sha[i] >= cutoff] <- 1; table(sha$X);
    chapuy$X <- 0;  chapuy$X[chapuy[i] >= cutoff] <- 1; table(chapuy$X)
    schmitz$X <- 0;  schmitz$X[schmitz[i] >= cutoff] <- 1; table(schmitz$X)
    
    #sha
    logistic_model <- glm(exitus_train ~ X, data = sha, family = "binomial")
    output[i,1] <- exp(summary(logistic_model)$coefficients[2])
    output[i,2] <- summary(logistic_model)$coefficients[8]
    #chapuy
    logistic_model <- glm(exitus_val ~ X, data = chapuy, family = "binomial")
    output[i,3] <- exp(summary(logistic_model)$coefficients[2])
    output[i,4] <- summary(logistic_model)$coefficients[8]
    #schmitz
    logistic_model <- glm(exitus_test ~ X, data = schmitz, family = "binomial")
    output[i,5] <- exp(summary(logistic_model)$coefficients[2])
    output[i,6] <- summary(logistic_model)$coefficients[8]
    
    colnames(sha)[colnames(sha) == "X"] <- paste0(i, "_cate");
    colnames(chapuy)[colnames(chapuy) == "X"] <- paste0(i, "_cate");
    colnames(schmitz)[colnames(schmitz) == "X"] <- paste0(i, "_cate");
    
  }
  
  # browser()
  
  output <- round(output, digits = 2)
  # write.csv(output, paste0(results_path,"/",name_row,".csv"), row.names=TRUE)
  return(dir)

}