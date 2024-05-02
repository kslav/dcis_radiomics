library(readxl)
library(readr)
library(modelr)
library(purrr)
library(caret)
library(pROC)
library(MASS)
library(splitstackshape)
library(epiR)


# data cleaning
# bsl_feature, size, surgery, and upgrade are data frames with relavent information
# cluster and pca are the phynotype information
data_cleaning <- function(bsl_feature, size, surgery, upgrade, cluster, pca, outcome){
  df <- merge(bsl_feature, surgery, by='case')
  df <- merge(df, size, by='case')
  df <- merge(df, upgrade, by='case')
  df <- merge(df, cluster, by='case')
  df <- merge(df, pca, by='case')
  
  df$ETHNICITY[df$ETHNICITY=='Unknown' | df$ETHNICITY=='Not Reported'] <- 'Unknown'
  df$RACE[df$RACE=='Unknown' | df$RACE=='Not Reported'] <- 'Unknown'
  df$ELIG_BIOPSY_DCIS_GRADE[df$ELIG_BIOPSY_DCIS_GRADE=='Missing' | df$ELIG_BIOPSY_DCIS_GRADE=='Grade cannot be assessed'] <- 'Unknown'
  df$PARENCHYMAL_ENHANCEMENT[df$PARENCHYMAL_ENHANCEMENT=='Missing'] <- 'Unknown'
  df$MRI_MORPHOLOGY[df$MRI_MORPHOLOGY=='Missing'] <- 'Unknown'
  
  if (outcome == 'upgrade'){
    df$UPGRADE_TO_INVASION[df$UPGRADE_TO_INVASION=='No pathology data'] <- NA
    df$UPGRADE_TO_INVASION <- ifelse(df$UPGRADE_TO_INVASION=='Yes', 1, 0)
    
    df_out <- df[complete.cases(df$UPGRADE_TO_INVASION), ]
    df_out$dcis_interrisk <- NULL
    df_out$dcis_highrisk <- NULL
    
    df_out$`k=2` <- as.factor(df_out$`k=2`)
    df_out$`k=3` <- as.factor(df_out$`k=3`)
  }else{
    df$dcis_interrisk <- ifelse(df$DCIS_SCORE>=39, 1, 0)
    df$dcis_highrisk <- ifelse(df$DCIS_SCORE>=55, 1, 0)
    
    df_out <- df[complete.cases(df$dcis_interrisk), ]
    df_out$UPGRADE_TO_INVASION <- NULL
    
    df_out$K <- as.factor(df_out$K)
  }
  
  df_out[sapply(df_out, is.character)] <- lapply(df_out[sapply(df_out, is.character)], as.factor)
  levels(df_out$RACE) <- c("Other", "Other", "Black or African American", "Other", "Unknown", "White")
  
  print(summary(df_out))
  
  return(df_out)
}


# marginal association analysis
# chi square test for marginal association between column x and y
association_analysis <- function(data, x_col, y_col){
  print(table(data[, x_col], data[, y_col]))
  chisq.test(data[, x_col], data[, y_col])
}


# model comparison
# compare two nested models using likelihood ratio test
model_comparison <- function(train_data, test_data, formula0, y_col, formula1, sens_threshold){
  true_outcome <- test_data[[y_col]]
  full_data <- rbind(train_data, test_data)
  
  # train and test models with training and testing data
  # threshold is selected such that the resulting sensitivities is the closest to sens_threshold
  model_training <- function(my_formula, train_data, test_data, true_outcome, sens_threshold){
    model <- glm(my_formula, data = train_data, family = 'binomial')
    pred <- predict(model, newdata=test_data, type = 'response')
    
    roc_obj <- roc(true_outcome, pred, quiet=TRUE)
    threshold <- mean(roc_obj$thresholds[which(abs(roc_obj$sensitivities - sens_threshold) ==
                                                 min(abs(roc_obj$sensitivities - sens_threshold)))])
    sens <- roc_obj$sensitivities[which(abs(roc_obj$sensitivities - sens_threshold) == 
                                          min(abs(roc_obj$sensitivities - sens_threshold)))][1]
    spec <- mean(roc_obj$specificities[which(abs(roc_obj$sensitivities - sens_threshold) == 
                                               min(abs(roc_obj$sensitivities - sens_threshold)))])
    ci <- round(ci.auc(roc_obj), 4)
    print(paste0('Threshold ', round(threshold, 4), ', auc (ci) ', ci[2], ' (', ci[1], ', ', ci[3], ')'))
    
    print(cm <- confusionMatrix(data=as.factor(ifelse(pred > threshold, 1, 0)), reference = as.factor(true_outcome), positive='1')$table)
    print(epi.tests(as.table(matrix(c(cm[2,2], cm[2,1], cm[1,2], cm[1,1]), nrow = 2, byrow = TRUE)), conf.level = 0.95, digits = 4))
    
    return(roc_obj)
  }
  
  print('Testing model 0')
  roc_obj0 <- model_training(paste(y_col, formula0), train_data, test_data, true_outcome, sens_threshold)

  print('Testing model 1')
  roc_obj1 <- model_training(paste(y_col, formula0, formula1), train_data, test_data, true_outcome, sens_threshold)
  
  # model comparison using all available data
  print(anova(glm(paste(y_col, formula0), data = full_data, family = 'binomial'), 
              glm(paste(y_col, formula0, formula1), data = full_data, family = 'binomial'), test='LRT'))
  
  return(list(roc_obj0, roc_obj1))
}


# example of comparing model C+M with model C+M+K for the primary outcome
association_analysis(data, 'k=2', 'UPGRADE_TO_INVASION')

set.seed(1)
obj <- stratified(data, 'DCIS_SCORE', 0.6, bothSets=TRUE) # training and testing set split using stratified sampling
train_data <- obj[[1]]
test_data <- obj[[2]]

out <- model_comparisont(train_data, test_data, 
                         ' ~ AGE+RACE+ETHNICITY+MRI_MORPHOLOGY+PARENCHYMAL_ENHANCEMENT+ELIG_BIOPSY_DCIS_GRADE+ELIG_BIOPSY_NECROSIS+ELIG_BIOPSY_ER_STATUS+ELIG_BIOPSY_PR_STATUS', 
                         'UPGRADE_TO_INVASION', ' + `k=2`', 0.9)
