# DEA OF BRCA

# Installing and Loading Libraries  

install.packages("dplyr")
library("dplyr")

install.packages("tidyverse")
library("tidyverse")

install.packages("ggrepel")
library("ggrepel")

if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("Glimma")
library("Glimma")

# Loading data

load("~/BRCA-TCGAA/BRCA-TCGAA/Data/brca_count.RData")

codes.type <- brca.clini %>% 
  filter(sample.type %in% c("TP", "NT")) %>%   
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  as_vector(.)

codes.typedf <- as.data.frame(codes.type)

# Running DESeq2 Differential Expression for sample.type TP (EA vs AA)

samplesTP <- brca.clini2 %>%
  filter(sample.type %in% c("TP"), eigenstrat %in% c("EA", "AA")) %>%  
  dplyr::select(c("sample.type","eigenstrat")) %>%
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  as_vector(.) %>%
  as.factor(.)

samplesTPdf <- as.data.frame(samplesTP)

ddsObj <- DESeqDataSetFromMatrix(countData = as.matrix(brca.count[,samplesTP]),
                                 colData = brca.clini2[samplesTP,],
                                 design = as.formula(~ eigenstrat))

ddsObjTP <- DESeq(ddsObjTP)

res.shrTP <- results(ddsObjTP)

summary(res.shrTP)

dea.TP <- as.data.frame(res.shrTP) %>%
  rownames_to_column("symbol") %>% 
  left_join(brca.annot, by="symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseqTP <- dea.TP  %>% filter(abs(logFC) >= 3, pvalue <= 0.05)
dim(df.deseqTP)
# [1] 1131    9
genes.DEA.TP.lst <- unique(df.deseqTP$symbol)
write(genes.DEA.TP.lst,  file = "genes.DEA.TP.lst")

save(dea.TP, file = "~/BRCA-TCGAA/BRCA-TCGAA/Data/dea.TP.rda", compress = T)

cutoffTP <- sort(dea.TP$pvalue)[10]
shrink.deseq.cutTP <- dea.TP %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffTP, symbol, ""))

ggplot(shrink.deseq.cutTP, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=0.5) +
  labs(x="mean of normalised counts", y="log Fold Change") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

cutoffTP <- sort(dea.NT.TP$pvalue)[10]
shrink.deseq.cutTP <- dea.TP %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffTP, symbol, ""))

ggplot(shrink.deseq.cutTP, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

# Running DESeq2 Differential Expression for sample.type NT (EA vs AA)

samplesNT <- brca.clini2 %>%
  filter(sample.type %in% c("NT"), eigenstrat %in% c("EA", "AA")) %>%  
  dplyr::select(c("sample.type","eigenstrat")) %>%
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  as_vector(.)

samplesNTdf <- as.data.frame(samplesNT)

ddsObjNT <- DESeqDataSetFromMatrix(countData = as.matrix(brca.count[,samplesNT]),
                                 colData = brca.clini2[samplesNT,],
                                 design = as.formula(~ eigenstrat))

ddsObjNT <- DESeq(ddsObjNT)

res.shrNT <- results(ddsObjNT)

summary(res.shrNT)

dea.NT <- as.data.frame(res.shrNT) %>%
  rownames_to_column("symbol") %>% 
  left_join(brca.annot, by="symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseqNT <- dea.NT  %>% filter(abs(logFC) >= 3, pvalue <= 0.05)
dim(df.deseq)
# [1] 1131    9
genes.DEA.NT.lst <- unique(df.deseqNT$symbol)
write(genes.DEA.NT.lst,  file = "genes.DEA.NT.lst")

save(dea.NT, file = "~/BRCA-TCGAA/BRCA-TCGAA/Data/dea.NT.rda", compress = T)

cutoffNT <- sort(dea.NT$pvalue)[10]
shrink.deseq.cutNT <- dea.NT %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, symbol, ""))

ggplot(shrink.deseq.cutNT, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=0.5) +
  labs(x="mean of normalised counts", y="log Fold Change") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

cutoffNT <- sort(dea.NT$pvalue)[10]
shrink.deseq.cutNT <- dea.NT %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, symbol, ""))

ggplot(shrink.deseq.cutNT, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

# Running DESeq2 Differential Expression for eigenstrat = EA

samplesEA <- brca.clini2 %>%
  filter(eigenstrat %in% c("EA")) %>%  
  dplyr::select(c("sample.type","eigenstrat")) %>%
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  as_vector(.)

samplesEAdf <- as.data.frame(samplesEA)

ddsObjEA <- DESeqDataSetFromMatrix(countData = as.matrix(brca.count[,samplesEA]),
                                   colData = brca.clini2[samplesEA,],
                                   design = as.formula(~ sample.type))

ddsObjEA <- DESeq(ddsObjEA)

res.shrEA <- results(ddsObjEA)

summary(res.shrEA)

dea.NT.TP.EA <- as.data.frame(res.shrEA) %>%
  rownames_to_column("symbol") %>% 
  left_join(brca.annot, by="symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseqEA <- dea.NT.TP.EA  %>% filter(abs(logFC) >= 3, pvalue <= 0.05)
dim(df.deseqEA)
# [1] 1108    9
genes.DEA.NT.vs.TP.lst.EA <- unique(df.deseqEA$symbol)
write(genes.DEA.NT.vs.TP.lst.EA,  file = "genes.DEA.NT.vs.TP.lst.EA")

save(dea.NT.TP.EA, file = "~/Projeto/BRCA-TCGAA/Data/dea.NT.TP.EA.rda", compress = T)

cutoffEA <- sort(dea.NT.TP.EA$pvalue)[10]
shrink.deseq.cutEA <- dea.NT.TP.EA %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffEA, symbol, ""))

ggplot(shrink.deseq.cutEA, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=0.5) +
  labs(x="mean of normalised counts", y="log Fold Change") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

cutoffEA <- sort(dea.NT.TP.EA$pvalue)[10]
shrink.deseq.cutEA <- dea.NT.TP.EA %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffEA, symbol, ""))

ggplot(shrink.deseq.cutEA, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

# Running DESeq2 Differential Expression for eigenstrat = aa

samplesAA <- brca.clini2 %>%
  filter(eigenstrat %in% c("AA")) %>%  
  dplyr::select(c("sample.type","eigenstrat")) %>%
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  as_vector(.)

samplesAAdf <- as.data.frame(samplesAA)

ddsObjAA <- DESeqDataSetFromMatrix(countData = as.matrix(brca.count[,samplesAA]),
                                   colData = brca.clini2[samplesAA,],
                                   design = as.formula(~ sample.type))

ddsObjAA <- DESeq(ddsObjAA)

res.shrAA <- results(ddsObjAA)

summary(res.shrAA)

dea.NT.TP.AA <- as.data.frame(res.shrAA) %>%
  rownames_to_column("symbol") %>% 
  left_join(brca.annot, by="symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseqAA <- dea.NT.TP.AA  %>% filter(abs(logFC) >= 3, pvalue <= 0.05)
dim(df.deseqAA)
# [1] 1072    9
genes.DEA.NT.vs.TP.lst.AA <- unique(df.deseqAA$symbol)
write(genes.DEA.NT.vs.TP.lst.AA,  file = "genes.DEA.NT.vs.TP.lst.AA")

save(dea.NT.TP.AA, file = "~/Projeto/BRCA-TCGAA/Data/dea.NT.TP.AA.rda", compress = T)

cutoffAA <- sort(dea.NT.TP.AA$pvalue)[10]
shrink.deseq.cutAA <- dea.NT.TP.AA %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffAA, symbol, ""))

ggplot(shrink.deseq.cutAA, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=0.5) +
  labs(x="mean of normalised counts", y="log Fold Change") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

cutoffAA <- sort(dea.NT.TP.AA$pvalue)[10]
shrink.deseq.cutAA <- dea.NT.TP.AA %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoffAA, symbol, ""))

ggplot(shrink.deseq.cutAA, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

