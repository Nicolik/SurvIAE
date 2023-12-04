#import libraries
library(survival)
library(survminer)
library(forestmodel)
library(limma)
library(ROCR)
library(data.table)
library(partykit)
library(caret)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(htmlTable)
library(magrittr)
library(formattable)
library(caTools)
library(epiDisplay)
library(sjPlot)
library(mltools)

#LOADING STARTING PATH
path <- "path/to/data"
setwd(path)

##################################################################################################################################################
# CLINICAL DATA
##################################################################################################################################################

# LOADING DS
myfile <- "./original_datasets/sha_clinic_info.csv" 
mydata_sha <- read.csv(myfile)
mydata_sha <-mydata_sha[order(mydata_sha$Accession),]
rownames(mydata_sha) <- mydata_sha$Accession

myfile <- "./original_datasets/chapuy_clinic_info.csv" 
mydata_chapuy <- read.csv(myfile)
mydata_chapuy <-mydata_chapuy[order(mydata_chapuy$NAME_short),]
rownames(mydata_chapuy) <- mydata_chapuy$NAME_short

myfile <- "./original_datasets/smitz_clinic_info.csv" 
mydata_schmitz <- read.csv(myfile)
mydata_schmitz <-mydata_schmitz[order(mydata_schmitz$dbGaP.subject.ID),]
rownames(mydata_schmitz) <- mydata_schmitz$dbGaP.subject.ID

# Outcomes standardization
colnames(mydata_sha)[colnames(mydata_sha) == "OS"] <- "global_OS_time"
colnames(mydata_sha)[colnames(mydata_sha) == "OS_status"] <- "global_OS"
colnames(mydata_sha)[colnames(mydata_sha) == "PFS"] <- "global_PFS_time"
colnames(mydata_sha)[colnames(mydata_sha) == "PFS_status"] <- "global_PFS"
mydata_sha$global_OS_time <- as.numeric(mydata_sha$global_OS_time); mydata_sha$global_OS <- as.numeric(mydata_sha$global_OS)
mydata_sha$global_PFS_time <- as.numeric(mydata_sha$global_PFS_time); mydata_sha$global_PFS <- as.numeric(mydata_sha$global_PFS)

colnames(mydata_chapuy)[colnames(mydata_chapuy) == "OS"] <- "global_OS_time"
colnames(mydata_chapuy)[colnames(mydata_chapuy) == "OS_STAT"] <- "global_OS"
colnames(mydata_chapuy)[colnames(mydata_chapuy) == "PFS"] <- "global_PFS_time"
colnames(mydata_chapuy)[colnames(mydata_chapuy) == "PFS_STAT"] <- "global_PFS"
mydata_chapuy <- subset(mydata_chapuy, global_OS != "na"); mydata_chapuy <- subset(mydata_chapuy, global_PFS != "na")
mydata_chapuy$global_OS_time <- as.numeric(mydata_chapuy$global_OS_time); mydata_chapuy$global_OS <- as.numeric(mydata_chapuy$global_OS)
mydata_chapuy$global_PFS_time <- as.numeric(mydata_chapuy$global_PFS_time); mydata_chapuy$global_PFS <- as.numeric(mydata_chapuy$global_PFS)

colnames(mydata_schmitz)[colnames(mydata_schmitz) == "OS_months"] <- "global_OS_time"
colnames(mydata_schmitz)[colnames(mydata_schmitz) == "OS_status"] <- "global_OS"
colnames(mydata_schmitz)[colnames(mydata_schmitz) == "PFS_months"] <- "global_PFS_time"
colnames(mydata_schmitz)[colnames(mydata_schmitz) == "PFS_status"] <- "global_PFS"
mydata_schmitz$global_OS_time <- as.numeric(mydata_schmitz$global_OS_time); mydata_schmitz$global_OS <- as.numeric(mydata_schmitz$global_OS)
mydata_schmitz$global_PFS_time <- as.numeric(mydata_schmitz$global_PFS_time); mydata_schmitz$global_PFS <- as.numeric(mydata_schmitz$global_PFS)

#dependent variables filtering and standardization
mydata_sha <- subset(mydata_sha, select=c('retrospective_COO_class','IPI_SCORE','global_OS','global_OS_time', 'global_PFS', 'global_PFS_time'))
#rev_ipi
mydata_sha$rev_ipi <- "NA"
mydata_sha$rev_ipi[mydata_sha$IPI_SCORE == "0" | mydata_sha$IPI_SCORE == "1" | mydata_sha$IPI_SCORE == "2"] <- "2.Good_and_VeryGood"
mydata_sha$rev_ipi[mydata_sha$IPI_SCORE == "3" | mydata_sha$IPI_SCORE == "4" | mydata_sha$IPI_SCORE == "5" ] <- "1.Poor"
mydata_sha$rev_ipi <- as.factor(mydata_sha$rev_ipi); table(mydata_sha$rev_ipi)
mydata_sha$IPI_SCORE <- NULL
#coo
colnames(mydata_sha)[colnames(mydata_sha) == "retrospective_COO_class"] <- "coo"
mydata_sha$coo <- as.factor(mydata_sha$coo); table(mydata_sha$coo)
levels(mydata_sha$coo) <- list("1.ABC+UNC" = "ABC", "1.ABC+UNC" = "UNC", "2.GCB" = "GCB"); table(mydata_sha$coo)

mydata_chapuy <- subset(mydata_chapuy, select=c('IPI','COO_byGEP','global_OS','global_OS_time', 'global_PFS', 'global_PFS_time'))
#rev_ipi
mydata_chapuy$rev_ipi <- "NA"
mydata_chapuy$rev_ipi[mydata_chapuy$IPI == "0" | mydata_chapuy$IPI == "1" | mydata_chapuy$IPI == "2"] <- "2.Good_and_VeryGood"
mydata_chapuy$rev_ipi[mydata_chapuy$IPI == "3" | mydata_chapuy$IPI == "4" | mydata_chapuy$IPI == "5" ] <- "1.Poor"
mydata_chapuy$rev_ipi <- as.factor(mydata_chapuy$rev_ipi); table(mydata_chapuy$rev_ipi)
mydata_chapuy$IPI <- NULL
#coo
colnames(mydata_chapuy)[colnames(mydata_chapuy) == "COO_byGEP"] <- "coo"
mydata_chapuy$coo <- as.factor(mydata_chapuy$coo); table(mydata_chapuy$coo)
levels(mydata_chapuy$coo) <- list("1.ABC+UNC" = "ABC", "1.ABC+UNC" = "UNC","2.GCB" = "GCB"); table(mydata_chapuy$coo)

mydata_schmitz <- subset(mydata_schmitz, select=c('IPI.Group','Gene.Expression.Subgroup','global_OS','global_OS_time', 'global_PFS', 'global_PFS_time'))
#rev_ipi
mydata_schmitz$rev_ipi <- "NA"
mydata_schmitz$rev_ipi[mydata_schmitz$IPI.Group == "Low"] <- "2.Good_and_VeryGood"
mydata_schmitz$rev_ipi[mydata_schmitz$IPI.Group == "Intermediate"] <- "2.Good_and_VeryGood"
mydata_schmitz$rev_ipi[mydata_schmitz$IPI.Group == "High"] <- "1.Poor"
table(mydata_schmitz$rev_ipi)
mydata_schmitz$IPI.Group <- NULL
#coo
colnames(mydata_schmitz)[colnames(mydata_schmitz) == "Gene.Expression.Subgroup"] <- "coo"
mydata_schmitz$coo <- as.factor(mydata_schmitz$coo); table(mydata_schmitz$coo)
levels(mydata_schmitz$coo) <- list("1.ABC+UNC" = "ABC","1.ABC+UNC" = "UC","2.GCB" = "GCB"); table(mydata_schmitz$coo)
# 
# #Pre-survival analysis - OVERALL SURVIVAL
# #Sha
# # rev-ipi
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(rev_ipi),data=mydata_sha); sopravv
# ggsurvplot(sopravv, data=mydata_sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                               palette= c("darkgreen","darkred"), legend.labs = c("P","G+VG"),
#                               legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12,font.y=12,
#                               font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(coo),data=mydata_sha); sopravv
# ggsurvplot(sopravv, data=mydata_sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# 
# #Chapuy
# # rev-ipi
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(rev_ipi),data=mydata_chapuy); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,130), pval=TRUE,ylab="Probability OS",
#                               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
#                               palette= c("darkgreen","darkred","gray"), legend.labs = c("P","G+VG","NA"),
#                               legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12,font.y=12,
#                               font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(coo),data=mydata_chapuy); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,130), pval=TRUE,ylab="Probability OS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# 
# #Schmitz
# # rev-ipi
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(rev_ipi),data=mydata_schmitz); sopravv
# ggsurvplot(sopravv, data=mydata_schmitz, risk.table=TRUE, conf.int=TRUE, xlim=c(0,200), pval=TRUE,ylab="Probability OS",
#                               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
#                               palette= c("darkgreen","darkred","gray"), legend.labs = c("P","G+VG","NA"),
#                               legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12,font.y=12,
#                               font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(OS_time,OS)~as.factor(coo),data=mydata_schmitz); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,200), pval=TRUE,ylab="Probability OS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# 
# #Pre-survival analysis - PROGRESSION FREE SURIVVAL
# #Sha
# # rev-ipi
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(rev_ipi),data=mydata_sha); sopravv
# ggsurvplot(sopravv, data=mydata_sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("darkgreen","darkred"), legend.labs = c("P","G+VG"),
#            legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(coo),data=mydata_sha); sopravv
# ggsurvplot(sopravv, data=mydata_sha, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Sha et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# 
# #Chapuy
# # rev-ipi
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(rev_ipi),data=mydata_chapuy); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,130), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
#            palette= c("darkgreen","darkred","gray"), legend.labs = c("P","G+VG","NA"),
#            legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(coo),data=mydata_chapuy); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,130), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Chapuy et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# 
# #Schmitz
# # rev-ipi
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(rev_ipi),data=mydata_schmitz); sopravv
# ggsurvplot(sopravv, data=mydata_schmitz, risk.table=TRUE, conf.int=TRUE, xlim=c(0,200), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=20,
#            palette= c("darkgreen","darkred","gray"), legend.labs = c("P","G+VG","NA"),
#            legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);
# # coo
# sopravv <- survfit(Surv(PFS_time,PFS)~as.factor(coo),data=mydata_schmitz); sopravv
# ggsurvplot(sopravv, data=mydata_chapuy, risk.table=TRUE, conf.int=TRUE, xlim=c(0,200), pval=TRUE,ylab="Probability PFS",
#            risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#            palette= c("green","pink"), legend.labs = c("ABC+UNC","GCB"),
#            legend="top", risk.table.title="N. at risk", title="Schmitz et al.", font.x=12,font.y=12,
#            font.tickslab=12,font.legend=14,pval.size=6.5,risk.table.height = 0.2);

#Remove general outcomes
# mydata_sha$PFS_time <- NULL; mydata_sha$PFS <- NULL; mydata_sha$OS <- NULL; mydata_sha$OS_time <- NULL
# mydata_chapuy$PFS_time <- NULL; mydata_chapuy$PFS <- NULL; mydata_chapuy$OS <- NULL; mydata_chapuy$OS_time <- NULL
# mydata_schmitz$PFS_time <- NULL; mydata_schmitz$PFS <- NULL; mydata_schmitz$OS <- NULL; mydata_schmitz$OS_time <- NULL

##################################################################################################################################################
# DEGS DATA
##################################################################################################################################################

# Creation of the output tables
output_odds <- matrix(0, 72, 8)
row.names(output_odds) <- c("NN_L_OS12", "NN_L_OS36", "NN_L_OS60", "NN_L_PFS12", "NN_L_PFS36", "NN_L_PFS60",
                       "NN_M_OS12", "NN_M_OS36", "NN_M_OS60", "NN_M_PFS12", "NN_M_PFS36", "NN_M_PFS60",
                       "NN_S_OS12", "NN_S_OS36", "NN_S_OS60", "NN_S_PFS12", "NN_S_PFS36", "NN_S_PFS60",
                       "NN_XS_OS12", "NN_XS_OS36", "NN_XS_OS60", "NN_XS_PFS12", "NN_XS_PFS36", "NN_XS_PFS60",
                       "N_L_OS12", "N_L_OS36", "N_L_OS60", "N_L_PFS12", "N_L_PFS36", "N_L_PFS60",
                       "N_M_OS12", "N_M_OS36", "N_M_OS60", "N_M_PFS12", "N_M_PFS36", "N_M_PFS60",
                       "N_S_OS12", "N_S_OS36", "N_S_OS60", "N_S_PFS12", "N_S_PFS36", "N_S_PFS60",
                       "N_XS_OS12", "N_XS_OS36", "N_XS_OS60", "N_XS_PFS12", "N_XS_PFS36", "N_XS_PFS60",
                       "STD_L_OS12", "STD_L_OS36", "STD_L_OS60", "STD_L_PFS12", "STD_L_PFS36", "STD_L_PFS60",
                       "STD_M_OS12", "STD_M_OS36", "STD_M_OS60", "STD_M_PFS12", "STD_M_PFS36", "STD_M_PFS60",
                       "STD_S_OS12", "STD_S_OS36", "STD_S_OS60", "STD_S_PFS12", "STD_S_PFS36", "STD_S_PFS60",
                       "STD_XS_OS12", "STD_XS_OS36", "STD_XS_OS60", "STD_XS_PFS12", "STD_XS_PFS36","STD_XS_PFS60")

output_performance <- matrix(0, 72, 6)
row.names(output_performance) <- c("NN_L_OS12", "NN_L_OS36", "NN_L_OS60", "NN_L_PFS12", "NN_L_PFS36", "NN_L_PFS60",
                            "NN_M_OS12", "NN_M_OS36", "NN_M_OS60", "NN_M_PFS12", "NN_M_PFS36", "NN_M_PFS60",
                            "NN_S_OS12", "NN_S_OS36", "NN_S_OS60", "NN_S_PFS12", "NN_S_PFS36", "NN_S_PFS60",
                            "NN_XS_OS12", "NN_XS_OS36", "NN_XS_OS60", "NN_XS_PFS12", "NN_XS_PFS36", "NN_XS_PFS60",
                            "N_L_OS12", "N_L_OS36", "N_L_OS60", "N_L_PFS12", "N_L_PFS36", "N_L_PFS60",
                            "N_M_OS12", "N_M_OS36", "N_M_OS60", "N_M_PFS12", "N_M_PFS36", "N_M_PFS60",
                            "N_S_OS12", "N_S_OS36", "N_S_OS60", "N_S_PFS12", "N_S_PFS36", "N_S_PFS60",
                            "N_XS_OS12", "N_XS_OS36", "N_XS_OS60", "N_XS_PFS12", "N_XS_PFS36", "N_XS_PFS60",
                            "STD_L_OS12", "STD_L_OS36", "STD_L_OS60", "STD_L_PFS12", "STD_L_PFS36", "STD_L_PFS60",
                            "STD_M_OS12", "STD_M_OS36", "STD_M_OS60", "STD_M_PFS12", "STD_M_PFS36", "STD_M_PFS60",
                            "STD_S_OS12", "STD_S_OS36", "STD_S_OS60", "STD_S_PFS12", "STD_S_PFS36", "STD_S_PFS60",
                            "STD_XS_OS12", "STD_XS_OS36", "STD_XS_OS60", "STD_XS_PFS12", "STD_XS_PFS36","STD_XS_PFS60")

#Move to subpath
subpath <- paste0(getwd(), "/genes-clf-v13")
setwd(subpath)

#Move to subpath
subdir_pre <- dir(path = ".", full.names = FALSE, recursive = FALSE)
subdir_pre <- subdir_pre[grepl(paste("genes-clf-", collapse = NULL),subdir_pre) == TRUE]
# subdir_pre <- subdir_pre[3] #For debug purposes
# j <- subdir_pre #For debug purposes

for (j in subdir_pre) {
  
  print(j)
  
  pre_proc_path <- paste0(subpath, "/", j) 
  if(pre_proc_path != getwd()){
    setwd(pre_proc_path)
  }
  
  subdirs <- dir(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  # subdirs <- subdirs[6] #For debug purposes
  
  for (i in subdirs) {
    
    print(i)
    
    outcome_path <- paste0(pre_proc_path, "/",i);
    if(outcome_path != getwd()){
      setwd(outcome_path)
    }
    
    #Function that load DEGs data for each pre-processing and for each outcome
    source(paste0(path,'/loading_deg_datasets.R'))
    deg_df_list <- loading_deg_datasets(outcome_path)
    deg_sha <- deg_df_list[[1]]; deg_chapuy <- deg_df_list[[2]]; deg_schmitz <- deg_df_list[[3]]
    remove(deg_df_list)
    
    if (!dir.exists("RESULTS")){
      dir.create("RESULTS")
    }

    files_risk <- list.files(path = ".", pattern = "results_RISK", all.files = FALSE, #Filtering of files with AE Risk
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    
    files_shap <- list.files(path = ".", pattern = "results_SHAP", all.files = FALSE, #Filtering of files with SHAP
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

    if(sum(grepl(paste("DimensionL", collapse = NULL),files_risk) == TRUE, na.rm=TRUE) == 3){
      subfiles_risk <- files_risk[grepl(paste("DimensionL", collapse = NULL),files_risk) == TRUE]
      subfiles_shap <- files_shap[grepl(paste("DimensionL", collapse = NULL),files_shap) == TRUE]
      arch <- "L_"; print(arch);

      # Sha
      myfile <- subfiles_risk[grepl(paste("sha", collapse = NULL),subfiles_risk) == TRUE]
      df_train <- read.csv(myfile)
      df_train$X <- NULL #Erase 1st column
      df_train <-df_train[order(df_train$Patient),]
      rownames(df_train) <- df_train$Patient
      df_train$Patient <- NULL #Erase 1st column
      if(nrow(df_train) == nrow(deg_sha) & nrow(df_train) ==  nrow(mydata_sha)) {
        df_train <- cbind(df_train,deg_sha, mydata_sha, make.row.names=T)
      } else {
        df_train$Patient <- rownames(df_train)
        deg_sha$Patient <- rownames(deg_sha)
        mydata_sha$Patient <- rownames(mydata_sha)
        df_train <- inner_join(df_train, deg_sha, by = "Patient")
        df_train <- inner_join(df_train, mydata_sha, by = "Patient")
        df_train$Patient <- NULL; deg_sha$Patient <- NULL; mydata_sha$Patient <- NULL
      }

      #Chapuy
      myfile <- subfiles_risk[grepl(paste("chapuy", collapse = NULL),subfiles_risk) == TRUE]
      df_val <- read.csv(myfile)
      df_val$X <- NULL #Erase 1st column
      df_val <-df_val[order(df_val$Patient),]
      rownames(df_val) <- df_val$Patient
      df_val$Patient <- NULL #Erase 1st column
      if(nrow(df_val) == nrow(deg_chapuy) & nrow(df_val) == nrow(mydata_chapuy)) {
        df_val <- cbind(df_val,deg_chapuy, make.row.names=T)
      } else {
        df_val$Patient <- rownames(df_val)
        deg_chapuy$Patient <- rownames(deg_chapuy)
        mydata_chapuy$Patient <- rownames(mydata_chapuy)
        df_val <- inner_join(df_val, deg_chapuy, by = "Patient")
        df_val <- inner_join(df_val, mydata_chapuy, by = "Patient")
        df_val$Patient <- NULL; deg_chapuy$Patient <- NULL; mydata_chapuy$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_1", collapse = NULL),subfiles_shap) == TRUE]
      shap_chapuy <- read.csv(myfile); shap_chapuy$X <- NULL #Erase 1st column

      #Schmitz
      myfile <- subfiles_risk[grepl(paste("schmitz", collapse = NULL),subfiles_risk) == TRUE]
      df_test <- read.csv(myfile)
      df_test$X <- NULL #Erase 1st column
      df_test <-df_test[order(df_test$Patient),]
      rownames(df_test) <- df_test$Patient
      df_test$Patient <- NULL #Erase 1st column
      if(nrow(df_test) == nrow(deg_schmitz) & nrow(df_test) == nrow(mydata_schmitz)) {
        df_test <- cbind(df_test,deg_schmitz, mydata_schmitz, make.row.names=T)
      } else {
        df_test$Patient <- rownames(df_test)
        deg_schmitz$Patient <- rownames(deg_schmitz)
        mydata_schmitz$Patient <- rownames(mydata_schmitz)
        df_test <- inner_join(df_test, deg_schmitz, by = "Patient")
        df_test <- inner_join(df_test, mydata_schmitz, by = "Patient")
        df_test$Patient <- NULL; deg_schmitz$Patient <- NULL; mydata_schmitz$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_2", collapse = NULL),subfiles_shap) == TRUE]
      shap_schmitz <- read.csv(myfile); shap_schmitz$X <- NULL #Erase 1st column

      if(setequal(shap_chapuy$Gene, shap_schmitz$Gene)==TRUE){
        shap_added <- shap_chapuy
        shap_added$SHAP.Value <- shap_chapuy$SHAP.Value + shap_schmitz$SHAP.Value
      }
      shap_added <-shap_added[order(shap_added$SHAP.Value, decreasing = TRUE),]
      list_best_genes <- shap_added$Gene[1:10]

      # source(paste0(path,'/survival_output.R'))
      # survival_output(outcome_path, df_train, df_val, df_test, arch)
      source(paste0(path,'/odds_output.R'))
      output_odds <- odds_output(outcome_path, df_train, df_val, df_test, arch, output_odds)
      source(paste0(path,'/performance_output.R'))
      output_performance <- performance_output(outcome_path, df_train, df_val, df_test, arch, output_performance)
      source(paste0(path,'/gene_sig.R'))
      outcome_path <- gene_sig(outcome_path, df_train, df_val, df_test, arch, list_best_genes)

      if(outcome_path != getwd()){
        setwd(outcome_path)
      }
    }

    if(sum(grepl(paste("DimensionM", collapse = NULL),files_risk) == TRUE, na.rm=TRUE) == 3){
      subfiles_risk <- files_risk[grepl(paste("DimensionM", collapse = NULL),files_risk) == TRUE]
      subfiles_shap <- files_shap[grepl(paste("DimensionM", collapse = NULL),files_shap) == TRUE]
      arch <- "M_"; print(arch);

      # Sha
      myfile <- subfiles_risk[grepl(paste("sha", collapse = NULL),subfiles_risk) == TRUE]
      df_train <- read.csv(myfile)
      df_train$X <- NULL #Erase 1st column
      df_train <-df_train[order(df_train$Patient),]
      rownames(df_train) <- df_train$Patient
      df_train$Patient <- NULL #Erase 1st column
      if(nrow(df_train) == nrow(deg_sha) & nrow(df_train) ==  nrow(mydata_sha)) {
        df_train <- cbind(df_train,deg_sha, mydata_sha, make.row.names=T)
      } else {
        df_train$Patient <- rownames(df_train)
        deg_sha$Patient <- rownames(deg_sha)
        mydata_sha$Patient <- rownames(mydata_sha)
        df_train <- inner_join(df_train, deg_sha, by = "Patient")
        df_train <- inner_join(df_train, mydata_sha, by = "Patient")
        df_train$Patient <- NULL; deg_sha$Patient <- NULL; mydata_sha$Patient <- NULL
      }

      #Chapuy
      myfile <- subfiles_risk[grepl(paste("chapuy", collapse = NULL),subfiles_risk) == TRUE]
      df_val <- read.csv(myfile)
      df_val$X <- NULL #Erase 1st column
      df_val <-df_val[order(df_val$Patient),]
      rownames(df_val) <- df_val$Patient
      df_val$Patient <- NULL #Erase 1st column
      if(nrow(df_val) == nrow(deg_chapuy) & nrow(df_val) == nrow(mydata_chapuy)) {
        df_val <- cbind(df_val,deg_chapuy, make.row.names=T)
      } else {
        df_val$Patient <- rownames(df_val)
        deg_chapuy$Patient <- rownames(deg_chapuy)
        mydata_chapuy$Patient <- rownames(mydata_chapuy)
        df_val <- inner_join(df_val, deg_chapuy, by = "Patient")
        df_val <- inner_join(df_val, mydata_chapuy, by = "Patient")
        df_val$Patient <- NULL; deg_chapuy$Patient <- NULL; mydata_chapuy$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_1", collapse = NULL),subfiles_shap) == TRUE]
      shap_chapuy <- read.csv(myfile); shap_chapuy$X <- NULL #Erase 1st column

      #Schmitz
      myfile <- subfiles_risk[grepl(paste("schmitz", collapse = NULL),subfiles_risk) == TRUE]
      df_test <- read.csv(myfile)
      df_test$X <- NULL #Erase 1st column
      df_test <-df_test[order(df_test$Patient),]
      rownames(df_test) <- df_test$Patient
      df_test$Patient <- NULL #Erase 1st column
      if(nrow(df_test) == nrow(deg_schmitz) & nrow(df_test) == nrow(mydata_schmitz)) {
        df_test <- cbind(df_test,deg_schmitz, mydata_schmitz, make.row.names=T)
      } else {
        df_test$Patient <- rownames(df_test)
        deg_schmitz$Patient <- rownames(deg_schmitz)
        mydata_schmitz$Patient <- rownames(mydata_schmitz)
        df_test <- inner_join(df_test, deg_schmitz, by = "Patient")
        df_test <- inner_join(df_test, mydata_schmitz, by = "Patient")
        df_test$Patient <- NULL; deg_schmitz$Patient <- NULL; mydata_schmitz$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_2", collapse = NULL),subfiles_shap) == TRUE]
      shap_schmitz <- read.csv(myfile); shap_schmitz$X <- NULL #Erase 1st column

      if(setequal(shap_chapuy$Gene, shap_schmitz$Gene)==TRUE){
        shap_added <- shap_chapuy
        shap_added$SHAP.Value <- shap_chapuy$SHAP.Value + shap_schmitz$SHAP.Value
      }
      shap_added <-shap_added[order(shap_added$SHAP.Value, decreasing = TRUE),]
      list_best_genes <- shap_added$Gene[1:10]

      source(paste0(path,'/survival_output.R'))
      survival_output(outcome_path, df_train, df_val, df_test, arch)
      source(paste0(path,'/odds_output.R'))
      output_odds <- odds_output(outcome_path, df_train, df_val, df_test, arch, output_odds)
      source(paste0(path,'/performance_output.R'))
      output_performance <- performance_output(outcome_path, df_train, df_val, df_test, arch, output_performance)
      source(paste0(path,'/gene_sig.R'))
      outcome_path <- gene_sig(outcome_path, df_train, df_val, df_test, arch, list_best_genes)

      if(outcome_path != getwd()){
        setwd(outcome_path)
      }
    }

    if(sum(grepl(paste("DimensionS", collapse = NULL),files_risk) == TRUE, na.rm=TRUE) == 3){
      subfiles_risk <- files_risk[grepl(paste("DimensionS", collapse = NULL),files_risk) == TRUE]
      subfiles_shap <- files_shap[grepl(paste("DimensionS", collapse = NULL),files_shap) == TRUE]
      arch <- "S_"; print(arch);
      
      # Sha
      myfile <- subfiles_risk[grepl(paste("sha", collapse = NULL),subfiles_risk) == TRUE]
      df_train <- read.csv(myfile)
      df_train$X <- NULL #Erase 1st column
      df_train <-df_train[order(df_train$Patient),]
      rownames(df_train) <- df_train$Patient
      df_train$Patient <- NULL #Erase 1st column
      if(nrow(df_train) == nrow(deg_sha) & nrow(df_train) ==  nrow(mydata_sha)) {
        df_train <- cbind(df_train,deg_sha, mydata_sha, make.row.names=T)
      } else {
        df_train$Patient <- rownames(df_train)
        deg_sha$Patient <- rownames(deg_sha)
        mydata_sha$Patient <- rownames(mydata_sha)
        df_train <- inner_join(df_train, deg_sha, by = "Patient")
        df_train <- inner_join(df_train, mydata_sha, by = "Patient")
        df_train$Patient <- NULL; deg_sha$Patient <- NULL; mydata_sha$Patient <- NULL
      }
      
      #Chapuy
      myfile <- subfiles_risk[grepl(paste("chapuy", collapse = NULL),subfiles_risk) == TRUE]
      df_val <- read.csv(myfile)
      df_val$X <- NULL #Erase 1st column
      df_val <-df_val[order(df_val$Patient),]
      rownames(df_val) <- df_val$Patient
      df_val$Patient <- NULL #Erase 1st column
      if(nrow(df_val) == nrow(deg_chapuy) & nrow(df_val) == nrow(mydata_chapuy)) {
        df_val <- cbind(df_val,deg_chapuy, make.row.names=T)
      } else {
        df_val$Patient <- rownames(df_val)
        deg_chapuy$Patient <- rownames(deg_chapuy)
        mydata_chapuy$Patient <- rownames(mydata_chapuy)
        df_val <- inner_join(df_val, deg_chapuy, by = "Patient")
        df_val <- inner_join(df_val, mydata_chapuy, by = "Patient")
        df_val$Patient <- NULL; deg_chapuy$Patient <- NULL; mydata_chapuy$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_1", collapse = NULL),subfiles_shap) == TRUE]
      shap_chapuy <- read.csv(myfile); shap_chapuy$X <- NULL #Erase 1st column
      
      #Schmitz
      myfile <- subfiles_risk[grepl(paste("schmitz", collapse = NULL),subfiles_risk) == TRUE]
      df_test <- read.csv(myfile)
      df_test$X <- NULL #Erase 1st column
      df_test <-df_test[order(df_test$Patient),]
      rownames(df_test) <- df_test$Patient
      df_test$Patient <- NULL #Erase 1st column
      if(nrow(df_test) == nrow(deg_schmitz) & nrow(df_test) == nrow(mydata_schmitz)) {
        df_test <- cbind(df_test,deg_schmitz, mydata_schmitz, make.row.names=T)
      } else {
        df_test$Patient <- rownames(df_test)
        deg_schmitz$Patient <- rownames(deg_schmitz)
        mydata_schmitz$Patient <- rownames(mydata_schmitz)
        df_test <- inner_join(df_test, deg_schmitz, by = "Patient")
        df_test <- inner_join(df_test, mydata_schmitz, by = "Patient")
        df_test$Patient <- NULL; deg_schmitz$Patient <- NULL; mydata_schmitz$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_2", collapse = NULL),subfiles_shap) == TRUE]
      shap_schmitz <- read.csv(myfile); shap_schmitz$X <- NULL #Erase 1st column
      
      if(setequal(shap_chapuy$Gene, shap_schmitz$Gene)==TRUE){
        shap_added <- shap_chapuy
        shap_added$SHAP.Value <- shap_chapuy$SHAP.Value + shap_schmitz$SHAP.Value
      }
      shap_added <-shap_added[order(shap_added$SHAP.Value, decreasing = TRUE),]
      list_best_genes <- shap_added$Gene[1:10]
      
      source(paste0(path,'/survival_output.R'))
      survival_output(outcome_path, df_train, df_val, df_test, arch)
      source(paste0(path,'/odds_output.R'))
      output_odds <- odds_output(outcome_path, df_train, df_val, df_test, arch, output_odds)
      source(paste0(path,'/performance_output.R'))
      output_performance <- performance_output(outcome_path, df_train, df_val, df_test, arch, output_performance)
      source(paste0(path,'/gene_sig.R'))
      outcome_path <- gene_sig(outcome_path, df_train, df_val, df_test, arch, list_best_genes)
      
      if(outcome_path != getwd()){
        setwd(outcome_path)
      }
    }

    if(sum(grepl(paste("DimensionXS", collapse = NULL),files_risk) == TRUE, na.rm=TRUE) == 3){
      subfiles_risk <- files_risk[grepl(paste("DimensionXS", collapse = NULL),files_risk) == TRUE]
      subfiles_shap <- files_shap[grepl(paste("DimensionXS", collapse = NULL),files_shap) == TRUE]
      arch <- "XS_"; print(arch);

      # Sha
      myfile <- subfiles_risk[grepl(paste("sha", collapse = NULL),subfiles_risk) == TRUE]
      df_train <- read.csv(myfile)
      df_train$X <- NULL #Erase 1st column
      df_train <-df_train[order(df_train$Patient),]
      rownames(df_train) <- df_train$Patient
      df_train$Patient <- NULL #Erase 1st column
      if(nrow(df_train) == nrow(deg_sha) & nrow(df_train) ==  nrow(mydata_sha)) {
        df_train <- cbind(df_train,deg_sha, mydata_sha, make.row.names=T)
      } else {
        df_train$Patient <- rownames(df_train)
        deg_sha$Patient <- rownames(deg_sha)
        mydata_sha$Patient <- rownames(mydata_sha)
        df_train <- inner_join(df_train, deg_sha, by = "Patient")
        df_train <- inner_join(df_train, mydata_sha, by = "Patient")
        df_train$Patient <- NULL; deg_sha$Patient <- NULL; mydata_sha$Patient <- NULL
      }

      #Chapuy
      myfile <- subfiles_risk[grepl(paste("chapuy", collapse = NULL),subfiles_risk) == TRUE]
      df_val <- read.csv(myfile)
      df_val$X <- NULL #Erase 1st column
      df_val <-df_val[order(df_val$Patient),]
      rownames(df_val) <- df_val$Patient
      df_val$Patient <- NULL #Erase 1st column
      if(nrow(df_val) == nrow(deg_chapuy) & nrow(df_val) == nrow(mydata_chapuy)) {
        df_val <- cbind(df_val,deg_chapuy, make.row.names=T)
      } else {
        df_val$Patient <- rownames(df_val)
        deg_chapuy$Patient <- rownames(deg_chapuy)
        mydata_chapuy$Patient <- rownames(mydata_chapuy)
        df_val <- inner_join(df_val, deg_chapuy, by = "Patient")
        df_val <- inner_join(df_val, mydata_chapuy, by = "Patient")
        df_val$Patient <- NULL; deg_chapuy$Patient <- NULL; mydata_chapuy$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_1", collapse = NULL),subfiles_shap) == TRUE]
      shap_chapuy <- read.csv(myfile); shap_chapuy$X <- NULL #Erase 1st column

      #Schmitz
      myfile <- subfiles_risk[grepl(paste("schmitz", collapse = NULL),subfiles_risk) == TRUE]
      df_test <- read.csv(myfile)
      df_test$X <- NULL #Erase 1st column
      df_test <-df_test[order(df_test$Patient),]
      rownames(df_test) <- df_test$Patient
      df_test$Patient <- NULL #Erase 1st column
      if(nrow(df_test) == nrow(deg_schmitz) & nrow(df_test) == nrow(mydata_schmitz)) {
        df_test <- cbind(df_test,deg_schmitz, mydata_schmitz, make.row.names=T)
      } else {
        df_test$Patient <- rownames(df_test)
        deg_schmitz$Patient <- rownames(deg_schmitz)
        mydata_schmitz$Patient <- rownames(mydata_schmitz)
        df_test <- inner_join(df_test, deg_schmitz, by = "Patient")
        df_test <- inner_join(df_test, mydata_schmitz, by = "Patient")
        df_test$Patient <- NULL; deg_schmitz$Patient <- NULL; mydata_schmitz$Patient <- NULL
      }
      myfile <- subfiles_shap[grepl(paste("val_2", collapse = NULL),subfiles_shap) == TRUE]
      shap_schmitz <- read.csv(myfile); shap_schmitz$X <- NULL #Erase 1st column

      if(setequal(shap_chapuy$Gene, shap_schmitz$Gene)==TRUE){
        shap_added <- shap_chapuy
        shap_added$SHAP.Value <- shap_chapuy$SHAP.Value + shap_schmitz$SHAP.Value
      }
      shap_added <-shap_added[order(shap_added$SHAP.Value, decreasing = TRUE),]
      list_best_genes <- shap_added$Gene[1:10]

      source(paste0(path,'/survival_output.R'))
      survival_output(outcome_path, df_train, df_val, df_test, arch)
      source(paste0(path,'/odds_output.R'))
      output_odds <- odds_output(outcome_path, df_train, df_val, df_test, arch, output_odds)
      source(paste0(path,'/performance_output.R'))
      output_performance <- performance_output(outcome_path, df_train, df_val, df_test, arch, output_performance)
      source(paste0(path,'/gene_sig.R'))
      outcome_path <- gene_sig(outcome_path, df_train, df_val, df_test, arch, list_best_genes)

      if(outcome_path != getwd()){
        setwd(outcome_path)
      }
    }

  } # Closing second for (outcomes)
} # Closing first for (pre-processing)

source(paste0(path,'/output_tables.R'))
output_tables(output_odds, output_performance)



