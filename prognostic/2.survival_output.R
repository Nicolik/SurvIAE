survival_output <- function(dir, sha, chapuy, schmitz, architecture){

  ##################################################################################################################################################
  # LOADING PATH
  ##################################################################################################################################################
  
  results_path <- paste0(dir, "/RESULTS"); 
  setwd(results_path); 
  
  ##################################################################################################################################################
  # DATA PREPARATION
  ##################################################################################################################################################
    
  # browser()
  
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
  # K-M AND COX REGRESSION
  ##################################################################################################################################################

  #Definition of lists for each panel construction
  splots_revipi <- list(); splots_coo <- list(); splots_risk <- list()

  #sha
  # #Analisi rev_ipi UNI
  sopravv <- survfit(Surv(time_train,exitus_train)~as.factor(rev_ipi),data=sha); #sopravv
  splots_revipi[[1]] <- ggsurvplot(sopravv, data=sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 70), pval=TRUE, ylab="Probability",
                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
                                   palette= c("darkgreen","darkred"), legend.labs = c("Poor","Good+VeryGood"),
                                   legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12, font.y=12,
                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); # splots_revipi[[1]]

  # #Analisi coo UNI
  sopravv <- survfit(Surv(time_train,exitus_train)~as.factor(coo),data=sha); #sopravv
  splots_coo[[1]] <- ggsurvplot(sopravv, data=sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 70), pval=TRUE, ylab="Probability",
                                risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
                                palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
                                legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12, font.y=12,
                                font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2);

  # #Analisi Latent Risk UNI
  sopravv <- survfit(Surv(time_train,exitus_train)~as.factor(Risk),data=sha); sopravv
  splots_risk[[1]] <- ggsurvplot(sopravv, data=sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 70), pval=TRUE, ylab="Probability",
                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
                                 palette= c("darkgray","darkblue"), legend.labs = c("C1","C2"),
                                 legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12, font.y=12,
                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); #splots_risk[[1]]

  #chapuy
  # #Analisi rev_ipi UNI
  sopravv <- survfit(Surv(time_val,exitus_val)~as.factor(rev_ipi),data=chapuy); sopravv
  splots_revipi[[2]] <- ggsurvplot(sopravv, data=chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 130), pval=TRUE, ylab="Probability",
                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
                                   palette= c("darkgreen","darkred"), legend.labs = c("Poor","Good+VeryGood"),
                                   legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12, font.y=12,
                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); #splots_revipi[[2]]

  # #Analisi coo UNI
  sopravv <- survfit(Surv(time_val,exitus_val)~as.factor(coo),data=chapuy); #sopravv
  splots_coo[[2]] <- ggsurvplot(sopravv, data=chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 130), pval=TRUE, ylab="Probability",
                                risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
                                palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
                                legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12, font.y=12,
                                font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2);

  # #Analisi Latent Risk UNI
  sopravv <- survfit(Surv(time_val,exitus_val)~as.factor(Risk),data=chapuy); sopravv
  splots_risk[[2]] <- ggsurvplot(sopravv, data=chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 130), pval=TRUE, ylab="Probability",
                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
                                 palette= c("darkgray","darkblue"), legend.labs = c("C1","C2"),
                                 legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12, font.y=12,
                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); #splots_risk[[2]]

  # cox <- coxph(formula=Surv(time_val,exitus_val) ~Risk+rev_ipi+coo, data=chapuy)
  # res <- ggforest(cox, data=chapuy, main = "Probability - Chapuy et al.")
  # ggsave(paste(name_row,"multi_Chapuy.pdf",sep="_"), res, width = 4000, height = 3000, units = "px")

  #schmitz
  # #Analisi rev_ipi UNI
  sopravv <- survfit(Surv(time_test,exitus_test)~as.factor(rev_ipi),data=schmitz); sopravv
  splots_revipi[[3]] <- ggsurvplot(sopravv, data=schmitz, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 200), pval=TRUE, ylab="Probability",
                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
                                   palette= c("darkgreen","darkred"), legend.labs = c("Poor","Good+VeryGood"),
                                   legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12, font.y=12,
                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); #splots_revipi[[3]]
  # #Analisi coo UNI
  sopravv <- survfit(Surv(time_test,exitus_test)~as.factor(coo),data=schmitz); #sopravv
  splots_coo[[3]] <- ggsurvplot(sopravv, data=schmitz, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 200), pval=TRUE, ylab="Probability",
                                risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
                                palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
                                legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12, font.y=12,
                                font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2);


  # #Analisi Latent Risk UNI
  sopravv <- survfit(Surv(time_test,exitus_test)~as.factor(Risk),data=schmitz); sopravv
  splots_risk[[3]] <- ggsurvplot(sopravv, data=schmitz, risk.table=TRUE, conf.int=TRUE, xlim=c(0, 200), pval=TRUE, ylab="Probability",
                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
                                 palette= c("darkgray","darkblue"), legend.labs = c("C1","C2"),
                                 legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12, font.y=12,
                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); #splots_risk[[3]]

  # cox <- coxph(formula=Surv(time_test,exitus_test) ~Risk+rev_ipi+coo, data=schmitz)
  # res <- ggforest(cox, data=schmitz, main = "Probability - Schmitz et al.")
  # ggsave(paste(name_row,"multi_Schmitz.pdf",sep="_"), res, width = 4000, height = 3000, units = "px")

  #Stampo pannello revipi
  res <- arrange_ggsurvplots(splots_revipi, print = FALSE, title = "Revised-IPI",
                             ncol = 3, nrow = 1, risk.table.height = 0.2)
  ggsave(paste(name_row,"Rev_IPI.pdf",sep="_"), res, width = 4000, height = 3000, units = "px")

  #Stampo pannello coo
  res <- arrange_ggsurvplots(splots_coo, print = FALSE, title = "COO",
                             ncol = 3, nrow = 1, risk.table.height = 0.2)
  ggsave(paste(name_row,"COO.pdf",sep="_"), res, width = 4000, height = 3000, units = "px")

  #Stampo pannello Risk
  res <- arrange_ggsurvplots(splots_risk, print = FALSE, title = "AE-Latent Clusters",
                             ncol = 3, nrow = 1, risk.table.height = 0.2)

  ggsave(paste(name_row,"Risk.pdf",sep="_"), res, width = 4000, height = 3000, units = "px")
  
  setwd(dir)

}