output_tables <- function(out_odds, out_perf){
  
  # setHtmlTableTheme(theme = "Google docs")
  
  #Table formatting
  out_odds <- round(out_odds, digits = 3)
  out_odds<- as.character(out_odds)
  out_odds[out_odds == '0'] <- "<0.001"
  out_odds[nchar(out_odds)>6] <- "Inf"
  
  
  # Formatted out_odds table 
  out_odds %>%
    addHtmlTableStyle(css.cgroup = "font-family: Arial, Helvetica, sans-serif; font-size: 115%",
                      css.rgroup = "font-family: Arial, Helvetica, sans-serif; font-style: italic; font-size: 105%",
                      css.header = "font-family: Arial, Helvetica, sans-serif; font-size: 105%",
                      css.cell = "font-family: Arial, Helvetica, sans-serif; font-size: 75%") %>% 
    htmlTable(header = c("OR", "p val", "OR", "p val", "OR", "p val", "OR", "p val"),
              rnames = c("L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S-OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60",
                         "L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S-OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60",
                         "L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S- OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60"),
              rgroup = c("Data not Normalized", "Normalized Data", "Standardized Data"),
              n.rgroup = c(24, 24, 24),
              cgroup = rbind(c("","Val Set","","Test Set"),
                             c("R-IPI", "AE-Risk", "R-IPI", "AE-Risk")),
              n.cgroup = rbind(c(2, 2, 2, 2),
                               c(2, 2, 2, 2)),
              caption = "",
              tfoot = "")
  
  #Table formatting
  out_perf <- round(out_perf, digits = 2)
  

  # Formatted out_odds table 
  out_perf %>%
    addHtmlTableStyle(css.cgroup = "font-family: Arial, Helvetica, sans-serif; font-size: 115%",
                      css.rgroup = "font-family: Arial, Helvetica, sans-serif; font-style: italic; font-size: 105%",
                      css.header = "font-family: Arial, Helvetica, sans-serif; font-size: 105%",
                      css.cell = "font-family: Arial, Helvetica, sans-serif; font-size: 75%") %>% 
    htmlTable(header = c("R-IPI", "COO", "AE-Risk", "R-IPI", "COO", "AE-Risk"),
              rnames = c("L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S-OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60",
                         "L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S-OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60",
                         "L-OS12", "L-OS36", "L-OS60", "L-PFS12","L-PFS36","L-PFS60",
                         "M-OS12", "M-OS36", "M-OS60", "M-PFS12","M-PFS36","M-PFS60",
                         "S-OS12", "S-OS36", "S- OS60", "S-PFS12","S-PFS36","S-PFS60",
                         "XS-OS12", "XS-OS36", "XS-OS60", "XS-PFS12","XS-PFS36","XS-PFS60"),
              rgroup = c("Data not Normalized", "Normalized Data", "Standardized Data"),
              n.rgroup = c(24, 24, 24),
              cgroup = c("Val Set","Test Set"),
              n.cgroup = c(3,3),
              caption = "",
              tfoot = "")
  
  browser()
  
  # SIGNIFICATIVE GENES FROM BEST MODELS
  d <- as.data.frame(output_performance)
  best_models <- NULL
  for (r in 1:nrow(output_performance)){
    if (d$V3[r]>d$V1[r] && d$V6[r]>d$V4[r]){
      best_models <- c(best_models, row.names(d)[r])
    }
  }
  
  remove(d)
  genes_tab <- data.frame(matrix(nrow = 10, ncol = length(best_models)))
  colnames(genes_tab) <- best_models
  color_matrix <- matrix(data=NA, nrow=10, ncol=length(best_models))
  splits <- NULL
  
  for (i in 1:length(best_models)){
    a <- strsplit(best_models[i], split = "_")
    splits <- c(splits, a)
    remove(a)
  }
  
  best_models <- data.frame(matrix(unlist(splits), ncol=3, byrow=T))
  remove(splits)
  
  idx <- 0
  
  for (i in colnames(genes_tab)){
    idx <- idx + 1
    
    if(best_models$X1[idx] == "NN"){
      setwd(paste0(subpath, "/genes-clf-nonorm"))
    } else   if(best_models$X1[idx] == "N"){
      setwd(paste0(subpath, "/genes-clf-minmax"))
    } else   if(best_models$X1[idx] == "STD"){
      setwd(paste0(subpath, "/genes-clf-standard"))
    } 
    
    if(best_models$X3[idx] == "OS12"){
      setwd(paste0(getwd(), "/OS12-DEG"))
    } else   if(best_models$X3[idx] == "OS36"){
      setwd(paste0(getwd(), "/OS36-DEG"))
    } else   if(best_models$X3[idx] == "OS60"){
      setwd(paste0(getwd(), "/OS60-DEG"))
    } else   if(best_models$X3[idx] == "PFS12"){
      setwd(paste0(getwd(), "/PFS12-DEG"))
    } else   if(best_models$X3[idx] == "PFS36"){
      setwd(paste0(getwd(), "/PFS36-DEG"))
    } else   if(best_models$X3[idx] == "PFS60"){
      setwd(paste0(getwd(), "/PFS60-DEG"))
    }
    
    if(best_models$X2[idx] == "L"){
      dimension <- "DimensionL"
    } else   if(best_models$X2[idx] == "M"){
      dimension <- "DimensionM"
    } else   if(best_models$X2[idx] == "S"){
      dimension <- "DimensionS"
    } else   if(best_models$X2[idx] == "XS"){
      dimension <- "DimensionXS"
    }
    subfiles_shap <- files_shap[grepl(paste(dimension, collapse = NULL),files_shap) == TRUE]
    
    myfile <- subfiles_shap[grepl(paste("val_1", collapse = NULL),subfiles_shap) == TRUE]
    shap_chapuy <- read.csv(myfile); shap_chapuy$X <- NULL
    
    myfile <- subfiles_shap[grepl(paste("val_2", collapse = NULL),subfiles_shap) == TRUE]
    shap_schmitz <- read.csv(myfile); shap_schmitz$X <- NULL
    
    if(setequal(shap_chapuy$Gene, shap_schmitz$Gene)==TRUE){
      shap_added <- shap_chapuy
      shap_added$SHAP.Value <- shap_chapuy$SHAP.Value + shap_schmitz$SHAP.Value
    }
    shap_added <-shap_added[order(shap_added$SHAP.Value, decreasing = TRUE),]
    list_best_genes <- shap_added$Gene[1:10]
    
    genes_tab[[i]] <- list_best_genes
    
    output <- read.csv(paste0("RESULTS/", i, ".csv"))
    
    for (k in 1:10){
      if(!(is.na(output$Chapuy_pval[k]) | is.na(output$Sha_pval[k]) | is.na(output$Sha_pval[k]))){
        if (output$Sha_pval[k]<0.1 && output$Chapuy_pval[k]<0.1 && output$Schmits_pval[k]<0.1){
          color_matrix[k,idx] <- 1
        }
        else
        {
          color_matrix[k,idx] <- 0
        }
      }
    }
    
  }
  
  color_formatter <- function(x, format = "color") {
    formatter("span",
              style = function(x) style(
                display = "block",
                "padding-top" = "2px",
                "padding-bottom" = "2px",
                "padding-left" = "5px",
                "padding-right" = "5px",
                "border-radius" = "2px",
                "color" = "black",
                "background-color" = ifelse(x == 1, "green",
                                            ifelse(x == 0, "red",
                                                   ifelse(is.na(x), "grey"))),
                "font-weight" = "bold"
              )
    )
  }
  
  # Applica la formattazione alle stringhe nel data.frame
  formatted_df <- genes_tab
  for (i in 1:nrow(genes_tab)) {
    for (j in 1:ncol(genes_tab)) {
      formatted_df[i, j] <- color_formatter(color_matrix[i, j])(color_matrix[i, j])
    }
  }
  
  # Visualizza il data.frame formattato
  formatted_df
  htmlTable(formatted_df)
  formatted_table <- formatted_df
  
  for (i in 1:nrow(genes_tab)) {
    for (j in 1:ncol(genes_tab)) {
      init <- strsplit(formatted_df[i, j], ">")[[1]][1]
      formatted_table[[i, j]] <- paste0(init, ">", genes_tab[i, j], "</span")
    }
  }
  
  formatted_table <- t(formatted_table)
  htmlTable(formatted_table)
  
}
  