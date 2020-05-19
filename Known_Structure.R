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

ug <- bnlearn::empty.graph(names(data))

bnlearn::arcs(ug, check.cycles = TRUE) = matrix(c( "NM.Hazard","Shape", "NM.Hazard","Nanoparticle", 
                                                   
                                                   "NM.Hazard", "Surface.area", "NM.Hazard",
                                                   
                                                   "Surface.charge", "NM.Hazard","Surface.coatings", "NM.Hazard", 
                                                   
                                                   "Surface.reactivity", "NM.Hazard", "Aggregation", "NM.Hazard", 
                                                   
                                                   "Particle.size", "NM.Hazard","Administration.route", "NM.Hazard",
                                                   
                                                   "Study.type", "NM.Hazard","Cytotoxicity", "NM.Hazard",
                                                   
                                                   "Neurological.effects", "NM.Hazard", "Pulmonary.effects", "NM.Hazard",
                                                   
                                                   "Fibrosis", "NM.Hazard", "RCNS.effects", "NM.Hazard", "Genotoxicity", "NM.Hazard",
                                                   
                                                   "Inflammation", "Nanoparticle","Shape",
                                                   
                                                   "Nanoparticle", "Surface.reactivity", "Nanoparticle", "Neurological.effects",
                                                   
                                                   "Nanoparticle", "Surface.coatings","Nanoparticle", "Surface.charge",
                                                   
                                                   "Nanoparticle", "Administration.route", "Nanoparticle", "Fibrosis",
                                                   
                                                   "Shape","Genotoxicity","Surface.area", "Neurological.effects",
                                                   
                                                   "Surface.coatings", "Surface.area", "Surface.coatings", "Particle.size", 
                                                   
                                                   "Surface.coatings", "Cytotoxicity", "Surface.coatings", "Pulmonary.effects", 
                                                   
                                                   "Surface.coatings", "Aggregation", "Surface.coatings", "Study.type",
                                                   
                                                   "Pulmonary.effects", "Inflammation","Inflammation","RCNS.effects"),
                                                
                                                ncol = 2, byrow = TRUE, dimnames = list(c(), c("from", "to")))


# five times cross-validation
total_acc <- rep(0,5) 
total_mcc <-rep(0,5)
total_confmatList <-list()

for(k in 1:5){

acc_best <- rep(0, N_folds)
mcc_best <- rep(0, N_folds)
confmatList <- list()
start <- 1

for(i in 1: N_folds){
  
  select <- start:(start+step-1)
  
  test <- data[select,]
  observed <- test$NM.Hazard
  #test[,12:18] <- NA  #without all these NAs the acc and mcc are better
  train <- data[-select,]
  
  training.fit <- bnlearn::bn.fit(ug, train, method="bayes") #method = "bayes"/"mle"
  
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(pred) <- levels(observed)
  
  for (j in 1:N_test){
    pred[j] <- as.character(predict(training.fit, data=test[j,!is.na(test[j,])], 
                                    
                                    node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
  }
  
  acc_best[i] <- sum(pred==observed)/N_test 
  mcc_best[i] <- yardstick::mcc_vec(observed, pred)
  confmatList[[i]] <- caret::confusionMatrix(pred, observed)
  
  #confmatcaret::confusionMatrix(pred, observed) #It holds the accuracy of the last prediction. This is the reason
  # why the last accuracy of acc_best is the same with the accuracy of confmat and why the mean of acc_best is not the same
  # with the accuracy of confmat
  
  start <- start+step
  
}

total_acc[k] <- mean(acc_best)
total_mcc[k] <- mean(mcc_best)
total_confmatList[[k]] <-confmatList # A list of lists. For example, for k=1 we have the 10 confusion matrices of the 1st cross
                                     # validation

}

total_acc_mean <- mean(total_acc)
total_mcc_mean <- mean(total_mcc)
