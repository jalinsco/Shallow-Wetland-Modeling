#############################################X
#------ Identify Shallow, Open Water --------X
#############################################X
#------ tunes & fits probability forests ----X
#############################################X


# Load packages ------------------

# required R packages
required_pkgs <- c("dplyr", 
                   "ggplot2", 
                   "readr", 
                   "sf", 
                   "blockCV", 
                   "ranger", 
                   "caret", 
                   "stringr")

# load 
missing <- required_pkgs[!vapply(required_pkgs, 
                                 requireNamespace, 
                                 logical(1), 
                                 quietly = TRUE)]

# check if missing
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}


# Load CV folds ------------------

# load spatial block CV info
load("data/Shallow Wetlands/checkerboard_fold_info_APR8.Rdata")
folds <- scv$folds_list

# inspect 
cv_plot(cv = scv) + geom_sf(color = "blue", fill = NA)


# Load training data ------------------

# un-comment below to choose the data frame for the appropriate model

# MODEL #1: the main model >>>
# uses custom superpixels as units of analysis
# includes PlanetScope-derived variables
sdm_rf <- read_csv("data/Shallow Wetlands/sw_model1_df.csv")

# MODEL #2: competing model >>>
# uses 30-m pixels (from Harmonized Landsat-Sentinel) as units of analysis
# same variable set, including PlanetScope-derived variables
#sdm_rf <- read_csv("data/Shallow Wetlands/sw_model2_df.csv")

# MODEL #3: competing model >>>
# uses 30-m pixels (from Harmonized Landsat-Sentinel) as units of analysis
# new variable set, excluding PlanetScope-derived variables
#sdm_rf <- read_csv("data/Shallow Wetlands/sw_model3_df.csv")


# Tuning & Cross-Validation ------------------

# this loop fits models with different hypertuning parameters
# predictions from the best model are used for cross-validation

# remove the lat/lon info 
sdm_rf <- sdm_rf %>% dplyr::select(-x, -y)

# create a data frame to store prediction of each CV fold
test_table <- data.frame(test_col = sdm_rf$test_col)
test_table$preds <- numeric(nrow(test_table))

# initialize a data frames to collect results acrpss folds 
all_hyperparams <- data.frame()  # to collect results across folds
imp_df <- data.frame(vars = colnames(sdm_rf)[2:(ncol(sdm_rf))])

# initialize a list to store models with different hypertuning 
mods <- list() 

# create vectors to store train/test split info for model tuning
fold_data <- data.frame(fold = seq(1:2), 
                        test.n = numeric(length(folds)),
                        train.n = numeric(length(folds)),
                        train.nSuitable = numeric(length(folds)),
                        test.nSuitable = numeric(length(folds)),
                        mtry = numeric(length(folds)),
                        min.node.size = numeric(length(folds)),
                        num.trees = numeric(length(folds)),
                        F1 = numeric(length(folds))
)

# function to summarize model fit
sum_fit <- function(data, lev = NULL, model = NULL) {
  confusion <- caret::confusionMatrix(data$pred, data$obs, positive = lev[1])
  accuracy <- confusion$overall["Accuracy"]
  precision <- confusion$byClass["Pos Pred Value"]
  recall <- confusion$byClass["Sensitivity"]
  f1 <- (2 * precision * recall) / (precision + recall)
  out <- c(Accuracy = accuracy, Precision = precision, Recall = recall, F1 = f1)
  names(out) <- c("Accuracy", "Precision", "Recall", "F1")
  return(out)
}


# loop to fit
# iter tuning & cross-validation 
# (approx 30 mins on 32GB RAM system)
set.seed(100124)
for(k in seq_len(length(folds))) {
  
  # ---------------X
  # prepare data >>>
  # ---------------X
  
  # extract the training and testing indices 
  train_idx <- unlist(folds[[k]][1]) # training set indices
  test_idx <- unlist(folds[[k]][2]) # testing set indices
  
  # create train and test sets 
  trainData <- sdm_rf[train_idx, ]
  testData <- sdm_rf[test_idx, ]
  
  # ---------------X
  # tune hyparams >>>
  # ---------------X
  
  # create a tuning grid
  tune_grid <- expand.grid(
    mtry = c(4, 5, 6, 7, 8, 9), # num var randomly sampled at each split
    splitrule = "hellinger", # splitting rule
    min.node.size = c(1, 2, 3, 4, 5, 6) # minimum size of terminal nodes
  )
  
  # trainControl
  train_control <- trainControl(
    method = "cv",
    number = 5,
    search = "grid",
    classProbs = TRUE,
    summaryFunction = sum_fit)
  
  # remake the dataset for tuning
  trainData_caret <- trainData %>% 
    mutate(test_col = as.factor(ifelse(test_col == 1, "suitable", "unsuitable"))) 
  
  # initialize
  num_trees_values = c(1400, 1600, 1800, 2000)
  best_f1 <- -Inf
  best_params <- NULL
  
  # run
  for (num_trees in num_trees_values) {
    rf_tuned <- caret::train(
      test_col ~ .,
      data = trainData_caret,
      method = "ranger",
      trControl = train_control,
      tuneGrid = tune_grid,
      num.trees = num_trees,  
      metric = "F1")  # performance metric
    
    # store all hyperparameter combinations and their performance
    fold_results <- rf_tuned$results
    fold_results$num.trees <- num_trees
    fold_results$fold <- k  
    all_hyperparams <- rbind(all_hyperparams, fold_results)
    
    # check if this model has the best F1 score so far
    new_f1 <- max(rf_tuned$results$F1, na.rm = TRUE)
    if (new_f1 > best_f1) {
      best_f1 <- new_f1
      best_params <- rf_tuned$bestTune
      best_params$num.trees <- num_trees  
    }
  }
  
  # ---------------X
  # fit model >>> 
  # ---------------X
  
  rf_mod <- ranger::ranger(test_col ~., 
                           data = trainData,
                           mtry = best_params$mtry,
                           splitrule = best_params$splitrule,
                           min.node.size = best_params$min.node.size,
                           num.trees = best_params$num.trees,
                           probability = TRUE,
                           importance = "permutation"
  ) 
  
  
  # ---------------X
  # save info >>>
  # ---------------X
  
  # predict probability of '1' in test set >>>
  # via ranger, with probability forest
  test_table$preds[test_idx] <- predict(rf_mod, data = testData)$predictions[,2] 
  
  # save variable importance
  imps <- ranger::importance(rf_mod) 
  imp_df[,(k+1)] <- imps
  
  # save fold data
  mods[[k]] <- rf_mod
  fold_data$train.n[k] <- nrow(trainData)
  fold_data$train.nSuitable[k] <- nrow(trainData[which(trainData$test_col == 1),])
  fold_data$test.n[k] <- nrow(testData)
  fold_data$test.nSuitable[k] <- nrow(trainData[which(testData$test_col == 1),])
  fold_data$mtry[k] <- best_params$mtry
  fold_data$min.node.size[k] <- best_params$min.node.size
  fold_data$num.trees[k] <- best_params$num.trees
  fold_data$F1[k] <- best_f1
  
  # clear
  rm(trainData_caret, trainData, train_idx, testData, test_idx)
  gc()
  
}

# print results for final fold
rf_mod

# check variable importance
colnames(imp_df)[2:3] <- c('fold.1', 'fold.2')
imp_df <- imp_df %>% 
  rowwise() %>%
  mutate(avg = (fold.1+fold.2)/2,
         sd = sd(c_across(fold.1:fold.2))) %>% 
  ungroup()
ggplot(imp_df, aes(x = reorder(vars, avg), y = avg)) +
  geom_bar(stat = "identity", fill = "#048BA8") +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Mean Importance (±sd)"
  ) +
  theme_minimal()

# set decision threshold 
op_th <- 0.47  

# confusion matrix to assess predictive performance
test_table <- test_table %>% 
  mutate(preds.dt = as.factor(ifelse(preds >= op_th, 1, 0)), 
         test_col = as.factor(test_col))
cM.dt <- caret::confusionMatrix(test_table$preds.dt, test_table$test_col, positive = '1')
head(test_table)

cM_df.dt <- as.data.frame(cM.dt$table)
cM_df.dt <- cM_df.dt %>%
  mutate(Prediction = factor(Prediction, levels = c(1, 0))) %>%
  group_by(Reference) %>% 
  mutate(
    total = sum(Freq),
    frac_fill = if_else(Prediction == Reference, Freq / total, 0),
    frac = Freq / total,
    lb = str_c(Freq, " (", round(frac * 100), "%)"))

ggplot(cM_df.dt, aes(Prediction, Reference, fill = frac_fill)) +
  geom_tile() +
  geom_text(aes(label = lb), size = 8) +
  scale_fill_gradient(low = "white", high = "#048BA8") +
  scale_x_discrete(position = "bottom") +
  geom_tile(color = "black", fill = "black", alpha = 0) +
  theme(axis.text = element_text(size = 10),
        title = element_text(size = 14),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(x = "Predicted Class", y = "True Class")


# Fit the final model ------------------

# this section fits a final version of model #1 for prediction
# uses all training data 

# identify the best parameters 
fold_data
idx <- which(fold_data$F1 == max(fold_data$F1))
best <- fold_data[idx,]
mtry <- best$mtry
min.node.size <- best$min.node.size
num.trees <- best$num.trees

# check parameters
mtry # 5
min.node.size # 1
num.trees # 2000

# train using all available data, tuned
rf_final <- ranger::ranger(test_col ~., 
                           data = sdm_rf,
                           num.trees = num.trees,
                           splitrule = "hellinger",
                           probability = TRUE,
                           mtry = mtry,
                           min.node.size = min.node.size,
                           importance = "permutation")

# print results
rf_final



