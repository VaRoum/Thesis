library(bnlearn)
library(yardstick)
library(caret)

data<-openxlsx::read.xlsx("C:/Users/Vaggelis/Desktop/Vagelis_Bayesian networks for NM Risk Assessment/MarvinsData/Copy of NM_data.xlsx",sheet=6,colNames=TRUE)

#convert all columns to factors
data[] <- lapply(data, as.factor) #lapply returns list so [] keeps the initial data format (data.frame)
drops<-c("Dissolution","Immunological.effects")
data<-data[,!names(data)%in%drops]
data<-data[duplicated(data)==F,]
str(data)

set.seed(1298)
data <- data[sample(nrow(data)),]
N_folds <- 10
N_data <- dim(data)[1]
step <- round(N_data/N_folds)


#five times cross-validation
#total_acc <- rep(0,5) 
#total_mcc <-rep(0,5)
#total_confmatList <-list()

#for(k in 1:5){

acc_best <- rep(0, N_folds)
mcc_best <- rep(0, N_folds)
confmatList <- list()
start <- 1


for(i in 1: N_folds){
  
  select <- start:(start+step-1)
  
  test <- data[select,]
  observed <- test$NM.Hazard
  #test[,12:18] <- NA # In most cases the results don't change with or without them
  train <- data[-select,]
  
  learned <- bnlearn::gs(train)
  learned1 <- cextend(learned, strict = TRUE) 
  
  training.fit <- bnlearn::bn.fit(learned1, train, method="bayes")
  
  train_imputed <- impute(training.fit, train)
  
  learned_imputed <- gs(train_imputed)
  learned_imputed1 <- cextend(learned, strict = TRUE)
  
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(learned_imputed1))
  graph::nodeRenderInfo(g) <- list(fontsize=50)
  Rgraphviz::renderGraph(g)
  
  training.fit1 <- bnlearn::bn.fit(learned_imputed1, train_imputed, method="bayes")
  
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(pred) <- levels(observed)
  
  for (j in 1:N_test){
    pred[j] <- as.character(predict(training.fit1, data=test[j,!is.na(test[j,])], 
                                    
                                    node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
  }
  
  acc_best[i] <- sum(pred==observed)/N_test 
  mcc_best[i] <- yardstick::mcc_vec(observed, pred)
  confmatList[[i]] <- caret::confusionMatrix(pred, observed)
  
  start <- start+step
  
}

number_of_na_acc <- sum(is.na(acc_best)==TRUE)
acc_best_without_na <- na.omit(acc_best)
number_of_na_mcc <- sum(is.na(mcc_best)==TRUE)
mcc_best_without_na <- na.omit(mcc_best)

mean_acc <- mean(acc_best_without_na)
mean_mcc <- mean(mcc_best_without_na)

#total_acc[k] <-mean(mean_acc)
#total_mcc[k] <-mean(mean_mcc)
#total_confmatList[[k]] <-confmatList # A list of lists. For example, for k=1 we have the 10 confusion matrices of the 1st cross
# validation

#}

#total_acc_mean <- mean(total_acc)
#total_mcc_mean <- mean(total_mcc)

# The results are the same either with 1 cross validation or with 5 cross validations