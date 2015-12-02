
library(rmongodb)

mongo <- mongo.create()
mongo
mongo.is.connected(mongo)
mongo.get.databases(mongo)
db <- 'exac'
mongo.get.database.collections(mongo, db)
#[1] "exac.dbsnp"         "exac.base_coverage" "exac.genes"
#[4] "exac.transcripts"   "exac.exons"         "exac.variants"

genes <- mongo.find.all(mongo, 'exac.genes')

gene.ids <- lapply(genes, function(x) x$gene_id)

buf <- mongo.bson.buffer.create()
mongo.bson.buffer.start.object(buf, "genes")
mongo.bson.buffer.append(buf, "$in", gene.ids)
mongo.bson.buffer.finish.object(buf)
criteria <- mongo.bson.from.buffer(buf)

variants.per.gene <- mongo.find.all(mongo, 'exac.variants', query=criteria)


#var user = db.user.findOne({"id" : "001"}, {"friends": 1})
#db.user.find( {"id" : {$in : user.friends }}).sort("age" : 1);

#TTN
query <- dbGetQuery(mg1, 'genes', "{'gene_id':'ENSG00000155657' }")

print(query)


# all genes
genes=[g for g in db.genes.find()]
gene_ids=[g['gene_id'] for g in genes]

hom_per_gene=dict()
for code in POPS: 
    pheno=POPS[code]
    print(code,pheno)
    hom_per_gene[code]=[ len([v['pop_homs'][pheno] for v in db.variants.find({'genes': gid}, fields={'_id': False})]) for gid in gene_ids ]


