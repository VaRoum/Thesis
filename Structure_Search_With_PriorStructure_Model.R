library(bnlearn)
library(yardstick)

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

for(i in 1:N_folds){
  
  select <- start:(start+step-1)
  test <- data[select,]
  
  observed <- test$NM.Hazard
  #test[,12:18] <- NA
  train <- data[-select,]
  
  learned <- bnlearn::structural.em(train, maximize = "tabu", maximize.args=list(score="bic"), 
                                    
                                    start=ug, fit="bayes", impute ="bayes-lw", impute.args = list(n=1000))
  
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(learned))
  graph::nodeRenderInfo(g) <- list(fontsize=50)
  Rgraphviz::renderGraph(g)
  
  training.fit <- bnlearn::bn.fit(learned, train, method="bayes") #method = "bayes"/"mle"
  
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
  start <- start+step
  
}

number_of_na_acc <- sum(is.na(acc_best)==TRUE)
acc_best_without_na <- na.omit(acc_best)
number_of_na_mcc <- sum(is.na(mcc_best)==TRUE)
mcc_best_without_na <- na.omit(mcc_best)

total_acc[k] <-mean(acc_best_without_na)
total_mcc[k] <-mean(mcc_best_without_na)
total_confmatList[[k]] <-confmatList # A list of lists. For example, for k=1 we have the 10 confusion matrices of the 1st cross
# validation

}

total_acc_mean <- mean(total_acc)
total_mcc_mean <- mean(total_mcc)

