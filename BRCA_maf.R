########## TCGA-BRCA Clinical Data ###########
# cleaning & exploring dataset with Tidyverse
# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/release/bioc/html/maftools.html

## Installing packages

packages_bioconductor <- c("TCGAbiolinks", "maftools", "BSgenome.Hsapiens.UCSC.hg38", "SummarizedExperiment")
packages_cran <- c("DT", "tidyverse", "data.table", "pheatmap", "NMF")

package.check <- lapply(packages_bioconductor, FUN = function(x) {
     if (!require(x, character.only = TRUE)) {
          BiocManager::install(x, dependencies = TRUE)
          library(x, character.only = TRUE)
     }
})
package.check <- lapply(packages_cran, FUN = function(x) {
     if (!require(x, character.only = TRUE)) {
          install.packages(x, dependencies = TRUE)
          library(x, character.only = TRUE)
     }
})

rm(packages_cran, packages_bioconductor, package.check)

if (!require('TCGAbiolinks')) {devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")} 

if (!require("BiocManager"))
     install.packages("BiocManager")
BiocManager::install("maftools")

## Reading Maf files 

maf <- GDCquery_Maf("BRCA", pipelines = "muse", directory = "GDCdata")
brca.maf <- read.maf(maf = maf, useAll = T)
write.mafSummary(maf = brca.maf, basename = 'brca.maf')

brca.clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical", save.csv = FALSE)
sort(colnames(brca.clinical))
colnames(brca.clinical)[1] <- "Tumor_Sample_Barcode"
# clinical$Overall_Survival_Status <- 1
# clinical$Overall_Survival_Status[which(clinical$vital_status == "Dead")] <- 0
# clinical$time <- clinical$days_to_death
# clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# create object for survival analysis 
brca.mafclin <- read.maf(maf = maf, clinicalData = brca.clinical, isTCGA = T)

## Visualization

plotmafSummary(maf = brca.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## Oncoplots

oncoplot(maf = brca.maf, top = 10)

## Transition and Transversions

brca.titv = titv(maf = brca.maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = brca.titv)

## Lollipop - ERRO

lollipopPlot(maf = brca.maf, gene = 'PIK3CA', AACol = 'AA_MAF', showDomainLabel = FALSE)

## Processing copy number data - ERRO

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

brca.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

getSampleSummary(brca.gistic) # xx samples (Tumor_Sample_Barcode) 
getGeneSummary(brca.gistic) # xx genes (hugo)
getCytobandSummary(brca.gistic) # xx cytobands

write.GisticSummary(gistic = brca.gistic, basename = 'brca.gistic')

gisticChromPlot(gistic = brca.gistic, markBands = "all")

gisticBubblePlot(gistic = brca.gistic)

gisticOncoPlot(gistic = brca.gistic, clinicalData = getClinicalData(x = brca.mafclin), clinicalFeatures = 'ajcc_pathologic_stage', sortByAnnotation = TRUE, top = 10)

## Analysis

# Somatic interations

somaticInteractions(maf = brca.maf, top = 25, pvalue = c(0.05, 0.1))

# Detecting cancer driver genes based on positional clustering

brca.sig = oncodrive(maf = brca.mafclin, AACol = 'Amino_acids', minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = brca.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# Survival analysis

mafSurvival(maf = brca.mafclin, genes = 'PIK3CA', time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)

prog_geneset = survGroup(maf = brca.mafclin, top = 20, geneSetSize = 2, time = "days_to_last_follow_up", Status = "vital_status", verbose = FALSE)

print(prog_geneset)

# Comparing two cohorts (MAFs)

primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
#Relapse APL MAF
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
relapse.apl = read.maf(maf = relapse.apl)
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt.vs.rt)

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

# Clinical enrichment analysis

fab.ce = clinicalEnrichment(maf = brca.mafclin, clinicalFeature = 'ajcc_pathologic_stage')

fab.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)