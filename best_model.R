library(bnlearn)
library(yardstick)

# Here we use Mathews Correlation Coefficient to compare the different methodologies.
# It is regarded as a balanced measure which can be used even if the classes are of 
# very different sizes.  A coefficient of +1 represents a perfect prediction, 0 no better 
# than random prediction and -1 indicates total disagreement between prediction and observation.
# When there are more than two labels the MCC will no longer range between -1 and +1. Instead the
# minimum  value will be between -1 and 0 depending on the true distribution.
# The maximum value is always +1. Another possible comparison metric is the Cohen's Kappa.

data <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx", 
                        sheetIndex=6, header=T)
# Remove variables with only one category
drops <- c("Dissolution","Immunological.effects")
data <- data[ , !(names(data) %in% drops)]
# Remove duplicates
data <- data[duplicated(data)==F,]
# Shuffle data
set.seed(1298)
data <- data[sample(nrow(data)),]
N_folds <- 10 
N_data  <- dim(data)[1]
step  <- round(N_data/N_folds)

############################3
### Random structure
#############################
acc_best <- rep(0, N_folds) 
start <- 1
# Test cross-val accuracy of best model
for(i in 1: N_folds){
  select <- start:(start+step-1)
  test <- data[select,] 
  #Extract observations
  observed <- test$NM.Hazard
  # Delete all observed values
  test[,12:18] <- NA
  train <- data[-select,]
  # fit the model
  training.fit  <- bnlearn::bn.fit(best.model, train, method="bayes") #method = "bayes"/"mle"
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  # Convert all results to factors in order to calculate the confusion matrix
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(pred) <- levels(observed)
  for (j in 1:N_test){
    # n : number of samples used to average
    pred[j] <-  as.character(predict(training.fit, data=test[j,!is.na(test[j,])], 
                                     node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
    #method = "bayes-lw" / "parents"
  }
  acc_best[i]  <- yardstick::mcc_vec(observed, pred)#sum(pred==observed)/N_test
  # Update starting position of test
  start <- start+step
}


##################################################
### Marvin model approach (fixed model structure)
##################################################
#create an empty network
ug <- bnlearn::empty.graph(names(data))

#add the arcs by assigning a two-column matrix containing the labels of their end-nodes.
#Undirected arcs are represented as their two possible orientations
bnlearn::arcs(ug, check.cycles = TRUE) = matrix(c( "NM.Hazard","Shape", "NM.Hazard","Nanoparticle", 
                                                   "NM.Hazard", "Surface.area",  "NM.Hazard",
                                                   "Surface.charge",  "NM.Hazard","Surface.coatings",  "NM.Hazard", 
                                                   "Surface.reactivity",  "NM.Hazard", "Aggregation",  "NM.Hazard", 
                                                   "Particle.size",  "NM.Hazard","Administration.route",  "NM.Hazard",
                                                   "Study.type",  "NM.Hazard","Cytotoxicity",  "NM.Hazard",
                                                   "Neurological.effects",  "NM.Hazard", "Pulmonary.effects",  "NM.Hazard",
                                                   "Fibrosis",  "NM.Hazard", "RCNS.effects",   "NM.Hazard", "Genotoxicity",  "NM.Hazard",
                                                   "Inflammation", "Nanoparticle","Shape",
                                                   "Nanoparticle",  "Surface.reactivity", "Nanoparticle", "Neurological.effects",
                                                   "Nanoparticle", "Surface.coatings","Nanoparticle", "Surface.charge",
                                                   "Nanoparticle", "Administration.route", "Nanoparticle",  "Fibrosis",
                                                   "Shape","Genotoxicity","Surface.area", "Neurological.effects",
                                                   "Surface.coatings", "Surface.area",  "Surface.coatings",  "Particle.size", 
                                                   "Surface.coatings", "Cytotoxicity", "Surface.coatings", "Pulmonary.effects", 
                                                   "Surface.coatings", "Aggregation",   "Surface.coatings",  "Study.type",
                                                   "Pulmonary.effects", "Inflammation","Inflammation","RCNS.effects"),
                                                ncol = 2, byrow = TRUE, dimnames = list(c(), c("from", "to")))
acc_marvin <- rep(0, N_folds) 
start <- 1
# Test cross-val accuracy of best model
for(i in 1: N_folds){
  select <- start:(start+step-1)
  test <- data[select,] 
  #Extract observations
  observed <- test$NM.Hazard
  # Delete all observed values
  test[,12:18] <- NA
  train <- data[-select,]
  # fit the model
  training.fit  <- bnlearn::bn.fit(ug, train, method="bayes") #method = "bayes"/"mle"
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  # Convert all results to factors in order to calculate the confusion matrix
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(pred) <- levels(observed)
  for (j in 1:N_test){
    # n : number of samples used to average
    pred[j] <-  as.character(predict(training.fit, data=test[j,!is.na(test[j,])], 
                                     node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
    #method = "bayes-lw" / "parents"
  }
  acc_marvin[i]  <- sum(pred==observed)/N_test#yardstick::mcc_vec(observed, pred)#
  # Update starting position of test
  start <- start+step
}


################################
### Structure search with prior
################################
prior.struct.search <- function(max,score, fit, impute, n, fit.method, pred.method){ 
acc_struct <- rep(0, N_folds) 
start <- 1
# Test cross-val accuracy of best model
for(i in 1:N_folds){
  select <- start:(start+step-1)
  test <- data[select,] 
  #Extract observations
  observed <- test$NM.Hazard
  # Delete all observed values
  test[,12:18] <- NA
  train <- data[-select,]
  # Find best structure
  learned  <- bnlearn::structural.em(train, maximize = max, maximize.args=list(score= score), 
                               start=ug, fit=fit, impute =impute, impute.args = list(n=n))
  # fit the model
  training.fit  <- bnlearn::bn.fit(learned, train, method=fit.method) #method = "bayes"/"mle"
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  # Convert all results to factors in order to calculate the confusion matrix
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(pred) <- levels(observed)
  for (j in 1:N_test){
    # n : number of samples used to average
    pred[j] <-  as.character(predict(training.fit, data=test[j,!is.na(test[j,])], 
                                     node= "NM.Hazard",prob = T, method = pred.method, n=10000))
    #method = "bayes-lw" / "parents"
  }
  acc_struct[i]  <- yardstick::mcc_vec(observed, pred)#sum(pred==observed)/N_test
 # print(paste0(acc_struct,sep=" ", i))
  # Update starting position of test
  start <- start+step
}
return(acc_struct)
}

# algorithms: "hc" / "tabu"
#scores: "loglik", "aic", "bic" (default), "bde", "bds", "mbde", "bdla", "k2"
#fitting method = "bayes"/"mle"
# predict methods:   "bayes-lw", "parents"

str_res1 <- rep(0,10)
for (i in 1:length(str_res1)){
  str_res1[i] <- mean(prior.struct.search("tabu","bic", "bayes", "bayes-lw", 100,"bayes", "bayes-lw"))
}

str_res2 <- rep(0,10)
for (i in 1:length(str_res2)){
  str_res2[i] <- mean(prior.struct.search("tabu","loglik", "bayes", "bayes-lw", 500,"bayes", "bayes-lw"))
}

str_res3 <- rep(0,10)
for (i in 1:length(str_res3)){
  str_res3[i] <- mean(prior.struct.search("tabu","bde", "bayes", "bayes-lw", 1000,"bayes", "bayes-lw"))
}

str_res4 <- rep(0,10)
for (i in 1:length(str_res4)){
  str_res4[i] <- mean(prior.struct.search("tabu","bds", "bayes", "bayes-lw", 5000,"bayes", "bayes-lw"))
}

res_mean <- cbind(mean(str_res1), mean(str_res2), mean(str_res3), mean(str_res4))
res_max <- cbind(max(str_res1), max(str_res2), max(str_res3), max(str_res4))
res_min <- cbind(min(str_res1), min(str_res2), min(str_res3), min(str_res4))

results <- data.frame(rbind(res_mean, res_max, res_min))
rownames(results) <- c("mean", "max", "min")
colnames(results) <- c("Approach1","Approach2","Approach3","Approach4")

#tabu is better than hc and BIC better than AIC
# n=1000 is better than 100, 500 and 5000
