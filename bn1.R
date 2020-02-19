library(bnlearn) # καλω το πακετο bnlearn
gr=empty.graph(LETTERS[1:12]) # καλω αδειο γραφο τυπου DAG με κομβους 12 κεφαλαια γραμματα
# το gr θα ανηκει στην class bn

gr
arc.set=matrix(c("A","B","B","C","B","D","C","E","D","E"),ncol=2,byrow=TRUE,dimnames=list(NULL,c("FROM","TO")))
# φτιαχνει εναν πινακα 2 στηλων δινοντας μια "θεωρητικη" μορφη στο μπεϋζιανο μου δικτυο

arc.set

arcs(e)=arc.set 
# συνεδεει με τοξα συνολο των κομβων με τον τροπο που εφτιαξα το set και με οσους κομβους
# χρησιμοποιησα απο το τον γραφο

e
plot(e)
