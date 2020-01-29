library(bnlearn)

setwd("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data")

#read data from sheet 1
  train <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                      sheetIndex=1,header=T)
#Sheet 3 contains validation data
test.init  <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                                sheetIndex=3,header=T)
#remove all duplicates from test
test <- test.init[duplicated(test.init)==F,] 
# Sheet 3 contains observed hazard values of the test set
results  <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                           sheetIndex=4,header=T)
# convert categories into strings
results <- sapply(results, as.character)
#select only the rows that are not duplicates
results <- results[duplicated(test.init)==F]

# Results presented in Marvin et al
marvin.res  <-  xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                                sheetIndex=5,header=T)
marvin.res  <- sapply(marvin.res[,5], as.character)
marvin.res <- marvin.res[duplicated(test.init)==F]

#get the number of rows of train and test
len_train = dim(train)[1]
len_test = dim(test)[1]
# Combine train and test in a single dataset
complete <- rbind(train, test)    
# get the factors of the complete dataset (level differences caused errors in predict())
complete[] <- lapply(complete, as.factor)
# Resplit the two set (train and test)
train = complete[1:len_train, ]
l = len_train+1
lf = len_train + len_test
test = complete[l:lf, ]
# Drop columns that contain only one column
drops <- c("Dissolution","Immunological.effects")
train <- train[ , !(names(train) %in% drops)]
test <- test[ , !(names(test) %in% drops)]


#create an empty network
ug <- bnlearn::empty.graph(names(train))

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

#plot the resulting graph
g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(ug))
graph::nodeRenderInfo(g) <- list(fontsize=50)
Rgraphviz::renderGraph(g)
# Parametric estimation
ug.learned  <- bn.fit(ug, train, method="bayes") #method = "bayes"/"mle"
N_test <- dim(test)[1]
pred <- rep(0, N_test)
for (i in 1:N_test){
   # Predictions for the Hazard node
    pred[i] <-  as.character(predict(ug.learned, data=test[i,!is.na(test[i,])], 
                                  node= "NM.Hazard",prob = T, method = "bayes-lw"))
  
}
sum(pred==results)
# Calculate accuracies
my.acc  <- sum(pred==results)/N_test
marvin.accu  <- sum(marvin.res==results)/N_test

# Convert all results to factors in order to calculate the confusion matrix
pred <- as.factor(pred) 
results <- as.factor(results)
marvin.res <- as.factor(marvin.res)
# confusion matrix
my.conf <- caret::confusionMatrix(pred, results)
marvin.conf <- caret::confusionMatrix(marvin.res, results)



##############################
##Structure learning with EM
##############################
#Learn the structure
#scores: "loglik", "aic", "bic" (default), "bde", "bds", "mbde", "bdla", "k2"
# fitting methods:   "bayes-lw", "parents"
#set.seed(1126)
learned  <- structural.em(train, maximize = "tabu", maximize.args=list(score= "bds"), start=ug,
                          fit="bayes", impute = "bayes-lw", impute.args = list(n=1000))
g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(learned))
graph::nodeRenderInfo(g) <- list(fontsize=50)
Rgraphviz::renderGraph(g)
#best results: tabu, bds, set.seed(1126)


training.fit = bn.fit(learned, train, method="bayes")
N_test <- dim(test)[1]
pred <- rep(0, N_test)
# Convert all results to factors in order to calculate the confusion matrix
pred <- as.factor(pred) 
results <- as.factor(results)
marvin.res <- as.factor(marvin.res)
levels(pred) <- levels(results)
for (i in 1:N_test){
  # n : number of samples used to average
  pred[i] <-  as.character(predict(training.fit, data=test[i,!is.na(test[i,])], 
                                   node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
  #method = "bayes-lw" / "parents"
}
sum(pred==results)

# Calculate accuracies
my.acc  <- sum(pred==results)/N_test
marvin.accu  <- sum(marvin.res==results)/N_test

# confusion matrix
my.conf <- caret::confusionMatrix(pred, results)
marvin.conf <- caret::confusionMatrix(marvin.res, results)

###########################
## Data imputation
##########################
impute.data  <- impute(training.fit, data=train, method = "bayes-lw")
struct.imputed  <- gs(train)

g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(struct.imputed))
graph::nodeRenderInfo(g) <- list(fontsize=50)
Rgraphviz::renderGraph(g)

training.fit2 = bn.fit(struct.imputed, train, method="bayes")
N_test <- dim(test)[1]
pred <- rep(0, N_test)
# Convert all results to factors in order to calculate the confusion matrix
pred <- as.factor(pred) 
results <- as.factor(results)
marvin.res <- as.factor(marvin.res)
levels(pred) <- levels(results)
for (i in 1:N_test){
  # n : number of samples used to average
  pred[i] <-  as.character(predict(training.fit2, data=test[i,!is.na(test[i,])], 
                                   node= "NM.Hazard",prob = T, method = "bayes-lw", n=10000))
  #method = "bayes-lw" / "parents"
}
sum(pred==results)