loading_deg_datasets <- function(dir){
  
  setwd(dir)
  #Loading of DEG datasets
  #Sha 1
  #DEG
  myfile <- "./input_files/X_train.csv"
  sha <- read.csv(myfile)
  # sha <-sha[order(sha$PatientID),]
  rownames(sha) <- sha$PatientID
  #time
  myfile <- "./input_files/t_train.csv"
  time_sha <- read.csv(myfile)
  # time_sha <-time_sha[order(time_sha$PatientID),]
  rownames(time_sha) <- time_sha$PatientID
  time_sha$PatientID <- NULL
  #outcome
  myfile <- "./input_files/y_train.csv"
  outcome_sha <- read.csv(myfile)
  # outcome_sha <-outcome_sha[order(outcome_sha$PatientID),]
  rownames(outcome_sha) <- outcome_sha$PatientID
  outcome_sha$PatientID <- NULL
  #Concatenatenation of 3 datasets
  sha <- cbind(sha, time_sha, outcome_sha)
  remove(time_sha,outcome_sha)
  #Sha 2
  myfile <- "./input_files/X_val.csv"
  sha_2 <- read.csv(myfile)
  # sha_2 <-sha_2[order(sha_2$PatientID),]
  rownames(sha_2) <- sha_2$PatientID
  #time
  myfile <- "./input_files/t_val.csv"
  time_sha_2 <- read.csv(myfile)
  # time_sha_2 <-time_sha_2[order(time_sha_2$PatientID),]
  rownames(time_sha_2) <- time_sha_2$PatientID
  time_sha_2$PatientID <- NULL
  #outcome
  myfile <- "./input_files/y_val.csv"
  outcome_sha_2 <- read.csv(myfile)
  # outcome_sha_2 <-outcome_sha_2[order(outcome_sha_2$PatientID),]
  rownames(outcome_sha_2) <- outcome_sha_2$PatientID
  outcome_sha_2$PatientID <- NULL
  #Concatenatenation of 3 datasets
  sha_2 <- cbind(sha_2, time_sha_2, outcome_sha_2)
  remove(time_sha_2,outcome_sha_2)
  #Concatenation of Sha and Sha_2
  sha <- rbind(sha, sha_2)
  remove(sha_2)
  sha <-sha[order(sha$PatientID),]
  sha$PatientID <- NULL
  
  
  #Chapuy
  #DEG
  myfile <- "./input_files/X_test.csv"
  chapuy <- read.csv(myfile)
  chapuy <-chapuy[order(chapuy$PatientID),]
  rownames(chapuy) <- chapuy$PatientID
  chapuy$PatientID <- NULL
  #time
  myfile <- "./input_files/t_test.csv"
  time_chapuy <- read.csv(myfile)
  time_chapuy <-time_chapuy[order(time_chapuy$PatientID),]
  rownames(time_chapuy) <- time_chapuy$PatientID
  time_chapuy$PatientID <- NULL
  #outcome
  myfile <- "./input_files/y_test.csv"
  outcome_chapuy <- read.csv(myfile)
  outcome_chapuy <-outcome_chapuy[order(outcome_chapuy$PatientID),]
  rownames(outcome_chapuy) <- outcome_chapuy$PatientID
  outcome_chapuy$PatientID <- NULL
  #Concatenatenation of 3 datasets
  chapuy <- cbind(chapuy, time_chapuy, outcome_chapuy)
  remove(time_chapuy,outcome_chapuy)
  
  #Schmitz
  #DEG
  myfile <- "./input_files/X_schmitz.csv"
  schmitz <- read.csv(myfile)
  colnames(schmitz)[1] <- "PatientID"
  schmitz <-schmitz[order(schmitz$PatientID),]
  rownames(schmitz) <- schmitz$PatientID
  schmitz$PatientID <- NULL
  #time
  myfile <- "./input_files/t_schmitz.csv"
  time_schmitz <- read.csv(myfile)
  colnames(time_schmitz)[1] <- "PatientID"
  time_schmitz <-time_schmitz[order(time_schmitz$PatientID),]
  rownames(time_schmitz) <- time_schmitz$PatientID
  time_schmitz$PatientID <- NULL
  #outcome
  myfile <- "./input_files/y_schmitz.csv"
  outcome_schmitz <- read.csv(myfile)
  colnames(outcome_schmitz)[1] <- "PatientID"
  outcome_schmitz <-outcome_schmitz[order(outcome_schmitz$PatientID),]
  rownames(outcome_schmitz) <- outcome_schmitz$PatientID
  outcome_schmitz$PatientID <- NULL
  #Concatenatenation of 3 datasets
  schmitz <- cbind(schmitz, time_schmitz, outcome_schmitz)
  remove(time_schmitz,outcome_schmitz)
  
  my_list <- list(sha, chapuy, schmitz)
  return(my_list)
}