
library(RMongo)

mg1 <- mongoDbConnect('exac')
print(dbShowCollections(mg1))



#TTN
query <- dbGetQuery(mg1, 'genes', "{'gene_id':'ENSG00000155657' }")

print(query)


