

install.packages("AnnotationHub")
library(AnnotationHub)
install.packages("GO.db")
library(GO.db)

ah <- AnnotationHub()

orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]  

#columns(orgdb)

df <- select(orgdb, gene_ls,keytype =   "SYMBOL",columns = "GO",multiVals = "first")


df2 <-select(GO.db,df$GO,c("TERM","DEFINITION"),"GOID")



df3 <- merge(df2,df,by.x="GOID",by.y="GO")

df4<-unique(df3)
