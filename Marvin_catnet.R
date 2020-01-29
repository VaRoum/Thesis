library(bnlearn)
library(catnet)
library(xlsx)
library(igraph)
library(Rgraphviz)
setwd("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data")

#read data from sheet 1
df <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                  sheetIndex=1,header=T)
# Sheet 2 counts the number of non-NAs of each column 
NA.count <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                            sheetIndex=2,header=F)
#Sheet 3 contains validation data
df.validate  <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                                sheetIndex=3,header=T)
df.res  <- xlsx::read.xlsx("/home/periklis/Desktop/PhD/Bayesian Networks/R scripts/Marvin_data/NM_data.xlsx",
                                sheetIndex=4,header=T)
df.res <- sapply(df.res, as.character)
str(df)

# Drop columns that contain only one column
#drops <- c("Dissolution","Immunological.effects", "Nanoparticle")
#df <- df[ , !(names(df) %in% drops)]

############################3
# UKNOWN ORDERING
############################

cnets <- catnet::cnSearchSA(data = df, selectMode = 'AIC')
#find network with maximum AIC 
aicnet <- catnet::cnFindAIC(object=cnets)
bicnet <- catnet::cnFindBIC(object=cnets)

#Produce a dot file containing the graph visualization from 'Rgraphviz' package
# to plot run on terminal: xdot filename.dot
catnet::cnDot(aicnet, "aicnet")

predictions <- catnet::cnPredict(aicnet,df.validate)
sum(predictions[,20]==df.res)/91
#cbind(predictions[,20]==df.res, df.res) check for which cats the preds fail



###############################
## TRY WITH ORDERED NODES
###############################

ordering <- c("NM.Hazard","Nanoparticle","Surface.coatings", "Surface.area","Shape",
              "Pulmonary.effects", "Inflammation","Dissolution","Surface.reactivity", 
              "Surface.charge",  
              "Study.type","Administration.route",
               "Particle.size","Aggregation","Genotoxicity", "Neurological.effects",
              "Immunological.effects", "Cytotoxicity", "Fibrosis",
               "RCNS.effects")
nodeOrder <- match(ordering,names(df))

netlist <- catnet::cnSearchOrder(data=df, perturbations=NULL, maxParentSet=7, parentSizes=2,
                         maxComplexity=5000, nodeOrder=nodeOrder)

aicnet <- catnet::cnFindAIC(object=netlist)

aicnet <- catnet::cnFind(object = netlist,complexity = 5000)
#Produce a dot file containing the graph visualization from 'Rgraphviz' package
# to plot run on terminal: xdot filename.dot
catnet::cnDot(aicnet, "newnet")

predictions <- catnet::cnPredict(aicnet,df.validate)
sum(predictions[,20]==df.res)/91


##########################
### Prefix the network
##########################

# Prefix parents of the network
fixparPool <- vector("list", dim(df)[2])

#Shape
fixparPool[[1]] <-  c(match("Nanoparticle" ,names(df)),match("NM.Hazard" ,names(df)) )
# Nanoparticle
fixparPool[[2]] <- c(match("NM.Hazard" ,names(df)))
# Dissolution
fixparPool[[3]] <- c(match("Nanoparticle" ,names(df)),match("NM.Hazard" ,names(df)) ) 
# Surface.area
fixparPool[[4]] <- c(match("Surface.coatings" ,names(df)), match("NM.Hazard" ,names(df)))
# Surface.charge
fixparPool[[5]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df))) 
# Surface.coatingss
fixparPool[[6]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df)) ) 
# Surface.reactivity
fixparPool[[7]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df))) 
# Aggregation
fixparPool[[8]] <- c(match("NM.Hazard" ,names(df)), match("Surface.coatings" ,names(df))) 
# Particle.size
fixparPool[[9]] <- c(match("NM.Hazard" ,names(df)), match("Surface.coatings" ,names(df)))
# Administration.route
fixparPool[[10]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df))) 
# Study.type
fixparPool[[11]] <- c(match("NM.Hazard" ,names(df)), match("Surface.coatings" ,names(df)))
# Cytotoxicity
fixparPool[[12]] <- c(match("NM.Hazard" ,names(df)), match("Surface.coatings" ,names(df)))
# Neurological.effects
fixparPool[[13]] <- c(match("NM.Hazard" ,names(df)), match("Surface.area" ,names(df))) 
# Pulmonary.effects
fixparPool[[14]] <- c(match("NM.Hazard" ,names(df)), match("Surface.coatings" ,names(df)))
# Fibrosis
fixparPool[[15]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df))) 
# RCNS.effects
fixparPool[[16]] <- c(match("NM.Hazard" ,names(df)), match("Inflammation" ,names(df))) 
# Immunological.effects
fixparPool[[17]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df))) 
# Genotoxicity
fixparPool[[18]] <- c(match("NM.Hazard" ,names(df)), match("Shape" ,names(df))) 
# Inflammation
fixparPool[[19]] <- c(match("NM.Hazard" ,names(df)), match("Pulmonary.effects" ,names(df))) 

cnet <- catnet::cnNew(nodes =names(df), cats = sapply(df, unique), parents = fixparPool)
norder <- catnet::cnOrder(cnet)

parPool <- vector("list", 20)
for(i in 1:20) parPool[[i]] <- c(match("NM.Hazard" ,names(df)), match("Nanoparticle" ,names(df)),
                              match("Surface.coatings" ,names(df)), match("Surface.area" ,names(df)),
                              match("Shape" ,names(df)),  match("Inflammation" ,names(df)),
                               match("Pulmonary.effects" ,names(df)) )

eval <-  catnet::cnSearchOrder(data=df, parentSizes=2, nodeOrder=nodeOrder, fixedParents=fixparPool)
eval
aicnet <- catnet::cnFind(object = eval, complexity = 100)
catnet::cnDot(aicnet, "eval")


predictions <- catnet::cnPredict(aicnet,df.validate[1:10,])
sum(predictions[,20]==df.res)/91

# Get CPTs
catnet::cnProb(aicnet)