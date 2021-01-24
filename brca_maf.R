########## TCGA-BRCA Clinical Data ###########
# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files

# Instalando e lendo pacotes

install.packages("BiocManager")
BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks")

library("maftools")
library("TCGAbiolinks")

# Lendo arquivos tipo MAF

maf <- GDCquery_Maf("BRCA", pipelines = "muse", directory = "GDCdata") 

maf_brca <- read.maf(maf = maf, useAll = T)
maf_brca # maf + tumor sample barcode

# Compreendedo o arquivo tipo MAF
getSampleSummary(maf_brca)
getGeneSummary(maf_brca)
getClinicalData(maf_brca)
getFields(maf_brca)

# Escrevendo o arquivo tipo MAF
write.mafSummary(maf = maf_brca, basename = 'maf_brca')

# Separando pacientes por ancestralidade
patients_AA <- scan("patients_AA", what=list(id=""), sep="\n") 
patients_EA <- scan("patients_EA", what=list(id=""), sep="\n")

maf$Patient_id <- substr(maf$Tumor_Sample_Barcode, 1, 12) 

samples_AA <- unique(maf$Tumor_Sample_Barcode[maf$Patient_id %in% patients_AA$id])
samples_EA <- unique(maf$Tumor_Sample_Barcode[maf$Patient_id %in% patients_EA$id])

maf_AA <- maf[maf$Tumor_Sample_Barcode %in% samples_AA,]
maf_EA <- maf[maf$Tumor_Sample_Barcode %in% samples_EA,]

maf_brca_AA <- read.maf(maf = maf_AA, useAll = T)
maf_brca_EA <- read.maf(maf = maf_EA, useAll = T)

# Coletando dados clínicos do TCGA 

all.clin <- read.csv('all_clin_indexed.csv', 
                     sep = ",",
                     header = T,
                     stringsAsFactors = T) # where does it come from?

all.clin$Patient_id <- all.clin$bcr_patient_barcode
all.clin$days_to_last_followup <- all.clin$days_to_last_follow_up
all.clin$Overall_Survival_Status <- all.clin$vital_status

brca_AA.clin <- all.clin[all.clin$bcr_patient_barcode %in% patients_AA$id,]
brca_EA.clin <- all.clin[all.clin$bcr_patient_barcode %in% patients_EA$id,]

brca_AA.clin$ancestry <- "AA"
brca_EA.clin$ancestry <- "EA"

brca_AA.clin <- merge(brca_AA.clin, maf_AA[,c("Patient_id","Tumor_Sample_Barcode")], by="Patient_id")
brca_EA.clin <- merge(brca_EA.clin, maf_EA[,c("Patient_id","Tumor_Sample_Barcode")], by="Patient_id")

maf_brca_AA <- read.maf(maf = maf_AA, clinicalData = brca_AA.clin, verbose = T) 
maf_brca_EA <- read.maf(maf = maf_EA, clinicalData = brca_EA.clin, verbose = T)

# brca.clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical", save.csv = FALSE)
# sort(colnames(brca.clinical))
# colnames(brca.clinical)[1] <- "Tumor_Sample_Barcode"
# 
# brca.clinical$Patient_id <- brca.clinical$bcr_patient_barcode
# brca.clinical$days_to_last_followup <- brca.clinical$days_to_last_follow_up
# brca.clinical$Overall_Survival_Status <- brca.clinical$vital_status
# 
# brca_AA.clin <- brca.clinical[brca.clinical$bcr_patient_barcode %in% patients_AA$id,]
# brca_EA.clin <- brca.clinical[brca.clinical$bcr_patient_barcode %in% patients_EA$id,]
# 
# brca_AA.clin$ancestry <- "AA"
# brca_EA.clin$ancestry <- "EA"
# 
# brca_AA.clin <- merge(brca_AA.clin, maf_AA[,c("Patient_id","Tumor_Sample_Barcode")], by="Patient_id")
# brca_EA.clin <- merge(brca_EA.clin, maf_EA[,c("Patient_id","Tumor_Sample_Barcode")], by="Patient_id")
# 
# maf_brca_AA <- read.maf(maf = maf_AA, clinicalData = brca_AA.clin, verbose = T) 
# maf_brca_EA <- read.maf(maf = maf_EA, clinicalData = brca_EA.clin, verbose = T) 

## Plot Summary 

# AA
plotmafSummary(maf = maf_brca_AA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# EA
plotmafSummary(maf = maf_brca_EA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## Oncoplot

# AA
oncoplot(maf = maf_brca_AA, top = 10)

# EA
oncoplot(maf = maf_brca_EA, top = 10)

## Somatic Interations 

# AA
somaticInteractions(maf = maf_brca_AA, top = 10, pvalue = c(0.05, 0.1))

# EA
somaticInteractions(maf = maf_brca_EA, top = 10, pvalue = c(0.05, 0.1))

## Driver genes

# AA
brca.sig = oncodrive(maf = maf_brca_AA, AACol = 'Amino_acids', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = brca.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# EA
brca.sig = oncodrive(maf = maf_brca_EA, AACol = 'Amino_acids', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = brca.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# Comparing two cohorts (MAFs)
AA.vs.EA <- mafCompare(m1 = maf_brca_AA, m2 = maf_brca_EA, m1Name = 'AfricanAmerican', m2Name = 'EuropeanAmerican', minMut = 5)

## Forest Plot

forestPlot(mafCompareRes = AA.vs.EA, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

# Clinical enrichment analysis - ERRO

fab.ce = clinicalEnrichment(maf = maf_brca_AA, clinicalFeature = 'Variant_Classification')

fab.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

## Oncogenic Signaling Pathways

# AA
OncogenicPathways(maf = maf_brca_AA)

# EA
OncogenicPathways(maf = maf_brca_EA)
