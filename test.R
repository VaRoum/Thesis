
library(bnlearn)
gr=empty.graph(LETTERS[1:5]) 
arc.set=matrix(c("A","B","B","C","B","D","C","E","D","E"),ncol=2,byrow=TRUE,
               dimnames=list(NULL,c("FROM","TO")))
e = gr
arcs(e)=arc.set 

par(mfrow=c(1,2)) 
plot(gr)
plot(e)
#this is a comment