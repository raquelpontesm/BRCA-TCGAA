if(!require(org.Hs.eg.db)){install.packages("org.Hs.eg.db")}
if(!require(clusterProfiler)){install.packages("clusterProfiler")}

geneslist_AA <- readLines("genes.DEA.NT.vs.TP.lst.AA")
geneslist_EA <- readLines("genes.DEA.NT.vs.TP.lst.EA")
geneslist_TP <- readLines("genes.DEA.TP.lst")

library("org.Hs.eg.db")
library("clusterProfiler")
library("AnnotationDbi")
library("ReactomePA")

geneslist_AA_maped <- mapIds(org.Hs.eg.db, keys = geneslist_AA, keytype = "SYMBOL", column="ENTREZID")
geneslist_EA_maped <- mapIds(org.Hs.eg.db, keys = geneslist_EA, keytype = "SYMBOL", column="ENTREZID")
geneslist_TP_maped <- mapIds(org.Hs.eg.db, keys = geneslist_TP, keytype = "SYMBOL", column="ENTREZID")

geneslist_AA_maped <- geneslist_AA_maped[!is.na(geneslist_AA_maped)]
geneslist_EA_maped <- geneslist_EA_maped[!is.na(geneslist_EA_maped)]
geneslist_TP_maped <- geneslist_TP_maped[!is.na(geneslist_TP_maped)]

geneslist <- list()
geneslist$AA  <- as.vector(geneslist_AA_maped)
geneslist$EA  <- as.vector(geneslist_EA_maped)
geneslist$TP  <- as.vector(geneslist_TP_maped) 

# https://rdrr.io/bioc/clusterProfiler/man/compareCluster.html
ck <- clusterProfiler::compareCluster(geneCluster = geneslist, 
                                      fun='groupGO', OrgDb='org.Hs.eg.db')

dotplot(ck)

ck <- clusterProfiler::compareCluster(geneCluster = geneslist, 
                                      fun="enrichKEGG",
                                      organism="hsa", pvalueCutoff=0.05)

dotplot(ck)

ck <- clusterProfiler::compareCluster(geneCluster = geneslist, 
                                      fun="enrichDO")

dotplot(ck)

ck <- clusterProfiler::compareCluster(geneCluster = geneslist, 
                                      fun="enrichPathway")

dotplot(ck)

# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#visualization
if(!require(gprofiler2)){install.packages("gprofiler2")}

# enrichment analysis
gp_AA = gost(geneslist_AA, 
                user_threshold = 0.05,
                organism = "hsapiens")

gostplot(gp_AA, capped = TRUE, interactive = TRUE)

#p1 = gostplot(gp_AA, interactive = FALSE)

#publish_gostplot(p1) #highlight_terms =c("CORUM:1514", "CORUM:2160", "CORUM:6087", "GO:1990405", 
#                                        "REAC:R-HSA-500792", "REAC:R-HSA-4164", "WP:WP176"))


# enrichment analysis
gp_EA = gost(geneslist_EA, 
             user_threshold = 0.05,
             organism = "hsapiens")

gostplot(gp_EA, capped = TRUE, interactive =  TRUE)

gp_TP = gost(geneslist_TP,
             user_threshold = 0.05,
             organism = "hsapiens")

gostplot(gp_TP, capped = TRUE, interactive = TRUE)

multi_gostres2 <- gost(query = list("EA" = geneslist_EA,
                                    "AA" = geneslist_AA,
                                    "TP" = geneslist_TP), 
                       multi_query = TRUE)

gostplot(multi_gostres2, capped = TRUE, interactive = TRUE)

