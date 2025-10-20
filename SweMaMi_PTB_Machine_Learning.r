library(ANCOMBC)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(DT)
library(mia)
library(caret)
library(microbiome)
library(microViz)
library(patchwork)
library(ggplot2)
library(mixOmics)
library(devtools)
library(vegan)
library(scutr)
library(BiocManager)
library(decontam)
library(mixOmics)
library(compositions)
library(MASS)
library(MLeval)

ancombc_input_files_csv <- list.files(path="/PATH", pattern=".csv", all.files=T, full.names=T)

for (files in ancombc_input_files_csv) {
    
    file_name <-  str_split(files,"/")
    file_name <- data.frame(file_name)
    file_name_2 <- str_split(file_name[10,],".csv")
    file_name_2 <- data.frame(file_name_2)
    name <- paste(file_name_2[1,],"",sep="")
    print(name)
    csv_file <- read.csv(files,sep=",")
    csv_file <- csv_file[,-c(1)]
    assign(name,csv_file)
}


metadata_Q1_for_ML <- read.csv("/PATH/metadata_Q1_for_ML.csv",sep=",")
metadata_Q2_for_ML <- read.csv("/PATH/metadata_Q2_for_ML.csv",sep=",")

metadata_Q2_for_ML[,c(1:38)] <- sapply(metadata_Q2_for_ML[,c(1:38)],function(x) as.character(x))

metadata_Q1_for_ML[,c(1:47)] <- sapply(metadata_Q1_for_ML[,c(1:47)],function(x) as.character(x))


# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(F2_case_control_clr)[-c(1, 2)]

# Add "_f" to the selected column names
names(F2_case_control_clr)[match(columns_to_modify, names(F2_case_control_clr))] <- paste0(columns_to_modify, "_f")


# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(F1_case_control_clr)[-c(1, 2)]

# Add "_f" to the selected column names
names(F1_case_control_clr)[match(columns_to_modify, names(F1_case_control_clr))] <- paste0(columns_to_modify, "_f")


# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(V1_case_control_clr)[-c(1, 2)]

# Add "_v" to the selected column names
names(V1_case_control_clr)[match(columns_to_modify, names(V1_case_control_clr))] <- paste0(columns_to_modify, "_v")

# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(V2_case_control_clr)[-c(1, 2)]

# Add "_v" to the selected column names
names(V2_case_control_clr)[match(columns_to_modify, names(V2_case_control_clr))] <- paste0(columns_to_modify, "_v")




for (i in 2:ncol(metadata_Q2_for_ML)) {
    print(i)
    tab_num <- table(metadata_Q2_for_ML[i])
    print(tab_num)
    }

## vaginal species + fecal species all T1
F1_case_control_clr_2 <- F1_case_control_clr[,-c(2)]
V1_F1_all_species <- merge(V1_case_control_clr,F1_case_control_clr_2,by="Studienummer") 
V1_F1_all_species <- na.omit(V1_F1_all_species)
V1_F1_all_species_2 <- V1_F1_all_species[,-c(1)]

## vaginal species + fecal species all + metadata Q1 T1
V1_F1_all_species_metadata <- merge(V1_F1_all_species,metadata_Q1_for_ML,by="Studienummer") 
V1_F1_all_species_metadata <- na.omit(V1_F1_all_species_metadata)
V1_F1_all_species_metadata_2 <- V1_F1_all_species_metadata[,-c(1)]

### metadata alone
T1_metadata_alone <- merge(V1_F1_all_species[,c(1:2)],metadata_Q1_for_ML,by="Studienummer")
T1_metadata_alone <- na.omit(T1_metadata_alone)
T1_metadata_alone_2 <- T1_metadata_alone[,-c(1)]

#write.csv(V1_F1_all_species_2,"/PATH/V1_F1_all_species.csv",row.names=FALSE)
#write.csv(V1_F1_all_species_metadata_2,"/PATH/V1_F1_all_species_metadata.csv",row.names=FALSE)
#write.csv(T1_metadata_alone_2,"/PATH/T1_metadata_alone.csv",row.names=FALSE)


## vaginal species + fecal species all T2
F2_case_control_clr_2 <- F2_case_control_clr[,-c(2)]
V2_F2_all_species <- merge(V2_case_control_clr,F2_case_control_clr_2,by="Studienummer") 
V2_F2_all_species <- na.omit(V2_F2_all_species)
V2_F2_all_species_2 <- V2_F2_all_species[,-c(1)]

## vaginal species + fecal species all + metadata Q2 T2
V2_F2_all_species_metadata <- merge(V2_F2_all_species,metadata_Q2_for_ML,by="Studienummer") 
V2_F2_all_species_metadata <- na.omit(V2_F2_all_species_metadata)
V2_F2_all_species_metadata_2 <- V2_F2_all_species_metadata[,-c(1)]

## metadata only for Q2
T2_metadata_alone <- merge(V2_F2_all_species[,c(1:2)],metadata_Q2_for_ML,by="Studienummer")
T2_metadata_alone <- na.omit(T2_metadata_alone)
T2_metadata_alone_2 <- T2_metadata_alone[,-c(1)]

#write.csv(V2_F2_all_species_2,"/PATH/V2_F2_all_species.csv",row.names=FALSE)
#write.csv(V2_F2_all_species_metadata_2,"/PATH/V2_F2_all_species_metadata.csv",row.names=FALSE)
#write.csv(T2_metadata_alone_2,"/PATH/T2_metadata_alone.csv",row.names=FALSE)


library(mixOmics) 
library(vctrs) 
library(caret)  
library(xgboost, Ckmeans.1d.dp) 
library(dplyr) 
library(ggplot2) 
library(stringi) 
library(arsenal) 
library(e1071) 
library(kernlab) 
library(nnet) 
library(randomForest) 
library(pROC) 
library(MLeval) 
library(skimr) 
library(RANN) 
library(tidyverse)
library(scutr) 
library(cluster) 
library(mclust) 

### tuning methods

train.control_1 <- trainControl(method = "repeatedcv", repeats=5, search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)
train.control_2 <- trainControl(method = "LOOCV", search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)
train.control_3 <- trainControl(method = "boot", search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)


## For metadata alone
case_control_clr_ML <- read.csv("/PATH/T1_metadata_alone.csv",sep=",")

### Repeated CV train.control_1

for (i in 1:10) {
    
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",115)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_1)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T1_metadata_only_training_rf_repeatcv_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T1_metadata_only_test_rf_repeatcv_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T1_metadata_only_rf_repeatcv_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T1_metadata_only_rf_repeatcv_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T1_metadata_only_loop_repeatedcv_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}

### LOOCV train.control_2

for (i in 1:10) {
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    random_number
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",115)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_2)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T1_metadata_only_training_rf_loocv_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T1_metadata_only_test_rf_loocv_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T1_metadata_only_rf_loocv_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T1_metadata_only_rf_loocv_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T1_metadata_only_loop_loocv_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}

### boot train.control_3

for (i in 1:10) {
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    random_number
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",115)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_3)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T1_metadata_only_training_rf_boot_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T1_metadata_only_test_rf_boot_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T1_metadata_only_rf_boot_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T1_metadata_only_rf_boot_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T1_metadata_only_loop_boot_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}


### time point 2 metadata only

case_control_clr_ML <- read.csv("/PATH/T2_metadata_alone.csv",sep=",")

### Repeated CV train.control_1

for (i in 1:10) {
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    random_number
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",90)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_1)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T2_metadata_only_training_rf_repeatcv_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T2_metadata_only_test_rf_repeatcv_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T2_metadata_only_rf_repeatcv_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T2_metadata_only_rf_repeatcv_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T2_metadata_only_loop_repeatedcv_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}

### LOOCV train.control_2

for (i in 1:10) {
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    random_number
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",90)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_2)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T2_metadata_only_training_rf_loocv_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T2_metadata_only_test_rf_loocv_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T2_metadata_only_rf_loocv_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T2_metadata_only_rf_loocv_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T2_metadata_only_loop_loocv_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}

### boot train.control_3

for (i in 1:10) {
    print(i)
    tryCatch(
    {
    random_number <- sample(1000:5000, 1)
    random_number
    set.seed(random_number)  
    control_subset <- resample_random(case_control_clr_ML, "Control", "key",90)

    case_set <- filter(case_control_clr_ML, key == "Case")
    subsampled <- rbind(control_subset,case_set)
    indexes <- createDataPartition(subsampled$key, #Column which contains class variable 
                               times = 1,
                               p = 0.8,
                               list = FALSE)
    train_dataset <- subsampled[indexes,]
    test_dataset <- subsampled[-indexes,]
    rf_training <- train(key ~ .,data = train_dataset, 
                             method = "rf", 
                             tuneLength = 7,
                             replace=F,
                             importance=TRUE,
                             trControl = train.control_3)


    prediction_rfMod<- predict(rf_training, newdata = test_dataset, type = "prob")
    
    train_name_dataset <- paste("/PATH/T2_metadata_only_training_rf_boot_dataframe_",i,".csv",sep="")
    test_name_dataset <- paste("/PATH/T2_metadata_only_test_rf_boot_dataframe_",i,".csv",sep="")
    training_name <- paste("/PATH/T2_metadata_only_rf_boot_training_model_",i,".rds",sep="")
    prediction_name <- paste("/PATH/T2_metadata_only_rf_boot_predict_",i,".rds",sep="")
    
    write.csv(train_dataset,train_name_dataset)
    write.csv(test_dataset,test_name_dataset)    
    saveRDS(rf_training,training_name)
    saveRDS(prediction_rfMod,prediction_name)


    print(paste("T2_metadata_only_loop_boot_",i,sep=""))
    },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}



input_files_rds <- list.files(path="/PATH", pattern=".rds", all.files=T, full.names=T)

for (files in input_files_rds) {
    
    file_name <-  str_split(files,"/")
    file_name <- data.frame(file_name)
    file_name_2 <- str_split(file_name[10,],".rds")
    file_name_2 <- data.frame(file_name_2)
    name <- paste(file_name_2[1,],"",sep="")
    print(name)
    rds_file <- readRDS(files)
    assign(name,rds_file)
}



input_files_csv <- list.files(path="/PATH", pattern=".csv", all.files=T, full.names=T)

for (files in input_files_csv) {
    
    file_name <-  str_split(files,"/")
    file_name <- data.frame(file_name)
    file_name_2 <- str_split(file_name[10,],".csv")
    file_name_2 <- data.frame(file_name_2)
    name <- paste(file_name_2[1,],"",sep="")
    print(name)
    csv_file <- read.csv(files,sep=",")
    csv_file <- csv_file[,-c(1)]
    assign(name,csv_file)
}

#### for metadata repeatcv

T1_metadata_only_repeatedcv_AUC_ROC <- NULL
for (i in 1:10) {
     tryCatch(
    {
    prediction_model <- paste("T1_metadata_only_rf_repeatcv_predict_",i,sep="")
    prediction_model_df <- get(prediction_model)
    test_set <- paste("T1_metadata_only_test_rf_repeatcv_dataframe_",i,sep="")
    test_set_df <- get(test_set)
    eval_test <- evalm(data.frame(prediction_model_df,test_set_df$key))
    AUC_values <- c(i,eval_test$optres$Group1[c('AUC-ROC'),c('Score')])
    AUC_vales <- data.frame(AUC_values)
    T1_metadata_only_repeatedcv_AUC_ROC <- rbind(T1_metadata_only_repeatedcv_AUC_ROC,AUC_values) },
    error = function(e) {
      cat("Error in iteration ", i, ": ", conditionMessage(e), "\n")
    }
  )
}
T1_metadata_only_repeatedcv_AUC_ROC <- data.frame(T1_metadata_only_repeatedcv_AUC_ROC)
T1_metadata_only_repeatedcv_AUC_ROC <- T1_metadata_only_repeatedcv_AUC_ROC[order(T1_metadata_only_repeatedcv_AUC_ROC$X2, decreasing = TRUE),]

best_iter <- T1_metadata_only_repeatedcv_AUC_ROC[1,1]
best_prediction_model <- paste("T1_metadata_only_rf_repeatcv_predict_",best_iter,sep="")
best_prediction_model_df <- get(best_prediction_model)
best_test_set <- paste("T1_metadata_only_test_rf_repeatcv_dataframe_",best_iter,sep="")
best_test_set_df <- get(best_test_set)
best_eval_species_T1_metadata_only_repeatedcv <- evalm(data.frame(best_prediction_model_df,best_test_set_df$key))


eval_T1_species_metadata <- evalm(data.frame(T1_all_species_metadata_rf_repeatcv_predict_4,T1_all_species_metadata_test_rf_repeatcv_dataframe_4))
T1_species_metadata_df <- as.data.frame(eval_T1_species_metadata$probs$Group1)
x_1 <- 1-T1_species_metadata_df$SPEC
y_1 <- T1_species_metadata_df$SENS

T1_metadata_only <- evalm(data.frame(T1_metadata_only_rf_repeatcv_predict_10,T1_metadata_only_test_rf_repeatcv_dataframe_10$key))
T1_metadata_only_df <- as.data.frame(T1_metadata_only$probs$Group1)
x_2 <- 1-T1_metadata_only_df$SPEC
y_2 <- T1_metadata_only_df$SENS



eval_V1_species_alone <- evalm(data.frame(T1_all_species_only_rf_repeatcv_predict_8,T1_all_species_only_test_rf_repeatcv_dataframe_8$key))
V1_species_alone_df <- as.data.frame(eval_V1_species_alone$probs$Group1)
x_4 <- 1-V1_species_alone_df$SPEC
y_4 <- V1_species_alone_df$SENS


eval_F1_species_alone <- evalm(data.frame(T1_all_species_only_rf_repeatcv_predict_1,T1_all_species_only_test_rf_repeatcv_dataframe_1$key))
F1_species_alone_df <- as.data.frame(eval_F1_species_alone$probs$Group1)
x_3 <- 1-F1_species_alone_df$SPEC
y_3 <- F1_species_alone_df$SENS

#### making figure

Figure <- ggplot() +
  geom_line(aes(x_1, y_1, color = "1. TP1 all taxa and metadata 0.77"), linewidth = 1) +
  geom_line(aes(x_2, y_2, color = "2. TP1 metadata only 0.74"), linewidth = 1) +
  geom_line(aes(x_4, y_4, color = "4. TP1 vaginal only 0.51"), linewidth = 1) +
  geom_line(aes(x_3, y_3, color = "3. TP1 fecal only 0.59"), linewidth = 1) +
  labs(title = "Initial AUC-ROC curves TP1", x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = c("1. TP1 all taxa and metadata 0.77" = "purple", "2. TP1 metadata only 0.74" = "red", "4. TP1 vaginal only 0.51" = "blue","3. TP1 fecal only 0.59" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1,size=18),# Adjust axis label size
    axis.title = element_text(size = 18),         # Adjust axis title size
    plot.title = element_text(size = 20, face = "bold"),  # Adjust plot title size
    panel.background = element_rect(fill = "white"),
       legend.text = element_text(size = 12))

#ggsave("/PATH/ML_TP1_1.pdf", plot = Figure, width = 12, height = 12, dpi = 600)



species_T1_model <- delta_vaginal_fecal_rf_repeatcv_training_model_4
species_T1_model_var_imp <- varImp(species_T1_model)
species_T1_model_var_imp_to_use <- species_T1_model_var_imp$importance[rowSums(species_T1_model_var_imp$importance != 0, na.rm = TRUE) > 0, ]


fecal_T1_model_var_imp_to_use <- species_T1_model_var_imp_to_use[grepl("_f",rownames(species_T1_model_var_imp_to_use)),]
fecal_T1_model_var_imp_to_use <- fecal_T1_model_var_imp_to_use[order(fecal_T1_model_var_imp_to_use$Case,decreasing=TRUE),]

vaginal_T1_model_var_imp_to_use <- species_T1_model_var_imp_to_use[grepl("_v",rownames(species_T1_model_var_imp_to_use)),]
vaginal_T1_model_var_imp_to_use <- vaginal_T1_model_var_imp_to_use[order(vaginal_T1_model_var_imp_to_use$Case,decreasing=TRUE),]
  
species_T1_model_var_imp_to_use_ordered <- species_T1_model_var_imp_to_use[order(species_T1_model_var_imp_to_use$Case,decreasing=TRUE),]
species_T1_model_var_imp_to_use_ordered_2 <- species_T1_model_var_imp_to_use_ordered[c(1:100),]
#write.csv(species_T1_model_var_imp_to_use_ordered_2,"/PATH/species_T1_model_var_imp_to_use_ordered.csv")


#fecal_T1_model_var_imp_to_use

fecal_T1_varimp <- fecal_T1_model_var_imp_to_use[c(1:10),]
colnames(fecal_T1_varimp) <- c("vars","value")

F1_varimp_plot <- ggplot(fecal_T1_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "fecal species", title = "fecal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines


#ggsave("/PATH/F1_varimp_plot.pdf", plot = F1_varimp_plot, width = 12, height = 12, dpi = 600)

#vaginal_T1_model_var_imp_to_use


vaginal_T1_varimp <- vaginal_T1_model_var_imp_to_use[c(1:10),]
colnames(vaginal_T1_varimp) <- c("vars","value")


V1_varimp_plot <- ggplot(vaginal_T1_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "vaginal species", title = "vaginal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

#ggsave("/PATH/V1_varimp_plot.pdf", plot = V1_varimp_plot, width = 12, height = 12, dpi = 600)


T1_all_species_metadata_test_rf_repeatcv_dataframe_4[,c(2281:2326)] 

vars_to_change <- c("Prev_IUFD", "Prev_PTB", "preeclampsia", "Any_pre_exist_diag", "Subfertility", "Age_groups", 
                    "Primipara", "treated_to_get_preg", "Ever_HPV_pos", "Ever_dysplasia", "HPV_vaccination", 
                    "antibiotic_last_3_months", "antibiotics_during_pregnancy", "cancer_before_pregnancy", 
                    "Diagnosed_eating_disorder", "Pregnancy_problems", "vitamins", "probiotics", "eat_fish", 
                    "Q1_prior_gi_medication", "Q1_during_gi_medication", "Q1_prior_neuro_medication", 
                    "Q1_during_neuro_medication", "Q1_any_drugs", "Q1_multiple_drugs", "miscarriage_previously")

T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2 <- T1_all_species_metadata_test_rf_repeatcv_dataframe_4
T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change] <- sapply(T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change],
                                                                                function(x) as.factor(x))
                                                                                  
T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change]                                                                                  

## T1 species metadata
eval_T1_species_metadata <- evalm(data.frame(T1_all_species_metadata_rf_repeatcv_predict_4,T1_all_species_metadata_test_rf_repeatcv_dataframe_4))

T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2 <- T1_all_species_metadata_test_rf_repeatcv_dataframe_4
T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change] <- sapply(T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change],
                                                                                function(x) as.factor(x))
                                                                                  
T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2[,vars_to_change]                                                                                  
CM_prediction_species_meta_T1 <- predict(T1_all_species_metadata_rf_repeatcv_training_model_4,newdata = T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2)

CM_species_meta_T1 <- confusionMatrix(CM_prediction_species_meta_T1,as.factor(T1_all_species_metadata_test_rf_repeatcv_dataframe_4_2$key))

## T1 metadata alone
T1_metadata_only <- evalm(data.frame(T1_metadata_only_rf_repeatcv_predict_10,T1_metadata_only_test_rf_repeatcv_dataframe_10$key))

CM_prediction_MD1 <- predict(T1_metadata_only_rf_repeatcv_training_model_10,newdata = T1_metadata_only_test_rf_repeatcv_dataframe_10)

CM_species_MD1 <- confusionMatrix(CM_prediction_MD1,as.factor(T1_metadata_only_test_rf_repeatcv_dataframe_10$key))
 

## V1 alone
eval_V1_species_alone <- evalm(data.frame(T1_all_species_only_rf_repeatcv_predict_8,T1_all_species_only_test_rf_repeatcv_dataframe_8$key))
 
CM_prediction_MD1 <- predict(T1_all_species_only_rf_repeatcv_training_model_8,newdata = T1_all_species_only_test_rf_repeatcv_dataframe_8)

CM_species_meta_V1 <- confusionMatrix(CM_prediction_MD1,as.factor(T1_all_species_only_test_rf_repeatcv_dataframe_8$key))


## F1 alone

eval_F1_species_alone <- evalm(data.frame(T1_all_species_only_rf_repeatcv_predict_1,T1_all_species_only_test_rf_repeatcv_dataframe_1$key))

CM_prediction_F1 <- predict(T1_all_species_only_rf_repeatcv_training_model_1 ,newdata = T1_all_species_only_test_rf_repeatcv_dataframe_1)

CM_species_meta_F1 <- confusionMatrix(CM_prediction_F1,as.factor(T1_all_species_only_test_rf_repeatcv_dataframe_1$key))


## V1 F1 alpha metadata
eval_V1_F1_clr_alpha_metadata <- evalm(data.frame(delta_vaginal_rf_repeatcv_predict_6,delta_vaginal_test_rf_repeatcv_dataframe_6$key))
delta_vaginal_test_rf_repeatcv_dataframe_6_2 <- delta_vaginal_test_rf_repeatcv_dataframe_6
delta_vaginal_test_rf_repeatcv_dataframe_6_2$Age_groups <- as.factor(delta_vaginal_test_rf_repeatcv_dataframe_6_2$Age_groups)
CM_prediction_TP1_alpha <- predict(delta_vaginal_rf_repeatcv_training_model_6 ,newdata = delta_vaginal_test_rf_repeatcv_dataframe_6_2)

CM_TP1_alpha <- confusionMatrix(CM_prediction_TP1_alpha,as.factor(delta_vaginal_test_rf_repeatcv_dataframe_6_2$key))

## T1 narrowed and metadata without alpha

eval_T1_var_imp_species_metadata_loocv <- evalm(data.frame(delta_vaginal_fecal_rf_repeatcv_predict_6,delta_vaginal_fecal_test_rf_loocv_dataframe_6$key))

delta_vaginal_fecal_test_rf_loocv_dataframe_6_2 <- delta_vaginal_fecal_test_rf_loocv_dataframe_6
delta_vaginal_fecal_test_rf_loocv_dataframe_6_2$Age_groups <- as.factor(delta_vaginal_fecal_test_rf_loocv_dataframe_6_2$Age_groups)


CM_prediction_TP1_varimp <- predict(delta_vaginal_fecal_rf_repeatcv_training_model_6 ,newdata = delta_vaginal_fecal_test_rf_loocv_dataframe_6_2)

CM_TP1_varimp <- confusionMatrix(CM_prediction_TP1_varimp,as.factor(delta_vaginal_fecal_test_rf_loocv_dataframe_6_2$key))

## All T2 taxa with metadata

eval_T2_species_metadata <- evalm(data.frame(T2_all_species_only_rf_loocv_predict_1,T2_all_species_only_test_rf_loocv_dataframe_1$key))

CM_prediction_T2_species_metadata <- predict(T2_all_species_only_rf_loocv_training_model_1 ,newdata = T2_all_species_only_test_rf_loocv_dataframe_1)

CM_T2_species_metadata_varimp <- confusionMatrix(CM_prediction_T2_species_metadata,as.factor(T2_all_species_only_test_rf_loocv_dataframe_1$key))

## V2 alone without metadata

species_only_eval_V2 <- evalm(data.frame(T2_all_species_metadata_rf_loocv_predict_5,T2_all_species_metadata_test_rf_loocv_dataframe_5$key))

CM_prediction_V2 <- predict(T2_all_species_metadata_rf_loocv_training_model_5 ,newdata = T2_all_species_metadata_test_rf_loocv_dataframe_5)

CM_V2 <- confusionMatrix(CM_prediction_V2,as.factor(T2_all_species_metadata_test_rf_loocv_dataframe_5$key))

## F2 alone without metadata

species_only_eval_F2 <- evalm(data.frame(T2_all_species_metadata_rf_loocv_predict_1,T2_all_species_metadata_test_rf_loocv_dataframe_1$key))

CM_prediction_F2 <- predict(T2_all_species_metadata_rf_loocv_training_model_1 ,newdata = T2_all_species_metadata_test_rf_loocv_dataframe_1)

CM_F2 <- confusionMatrix(CM_prediction_F2,as.factor(T2_all_species_metadata_test_rf_loocv_dataframe_1$key))

## T2 metadata alone

T2_metadata_only <- evalm(data.frame(T2_metadata_only_rf_boot_predict_7,T2_metadata_only_test_rf_boot_dataframe_7$key))

CM_prediction_MD2 <- predict(T2_metadata_only_rf_boot_training_model_1 ,newdata = T2_metadata_only_test_rf_boot_dataframe_7)

CM_MD2 <- confusionMatrix(CM_prediction_MD2,as.factor(T2_metadata_only_test_rf_boot_dataframe_7$key))

### combined var imp metadata and species T2 without alpha
T2_all_species_only_var_imp_no_alpha <- evalm(data.frame(delta_vaginal_fecal_rf_repeatcv_predict_1,delta_vaginal_fecal_test_rf_repeatcv_dataframe_1$key))

delta_vaginal_fecal_test_rf_repeatcv_dataframe_1_2 <- delta_vaginal_fecal_test_rf_repeatcv_dataframe_1
delta_vaginal_fecal_test_rf_repeatcv_dataframe_1_2$Age_groups <- as.factor(delta_vaginal_fecal_test_rf_repeatcv_dataframe_1_2$Age_groups)



CM_prediction_T2_varimp <- predict(delta_vaginal_fecal_rf_repeatcv_training_model_1 ,newdata = delta_vaginal_fecal_test_rf_repeatcv_dataframe_1_2)

CM_T2_varimp <- confusionMatrix(CM_prediction_T2_varimp,as.factor(delta_vaginal_fecal_test_rf_repeatcv_dataframe_1_2$key))

### combined var imp metadata and species and alpha diversity T2
V2_F2_clr_alpha_metadata <- evalm(data.frame(delta_vaginal_fecal_rf_repeatcv_predict_4,delta_vaginal_fecal_test_rf_repeatcv_dataframe_4$key))

delta_vaginal_fecal_test_rf_repeatcv_dataframe_4_2 <- delta_vaginal_fecal_test_rf_repeatcv_dataframe_4
delta_vaginal_fecal_test_rf_repeatcv_dataframe_4_2$Age_groups <- as.factor(delta_vaginal_fecal_test_rf_repeatcv_dataframe_4_2$Age_groups)


CM_prediction_T2_alpha <- predict(delta_vaginal_fecal_rf_repeatcv_training_model_4 ,newdata = delta_vaginal_fecal_test_rf_repeatcv_dataframe_4_2)

CM_T2_alpha <- confusionMatrix(CM_prediction_T2_alpha,as.factor(delta_vaginal_fecal_test_rf_repeatcv_dataframe_4_2$key))

