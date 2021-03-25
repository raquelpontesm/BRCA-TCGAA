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

# Running DESeq2 Differential Expression for sample.type

ddsObj <- DESeqDataSetFromMatrix(countData = as.matrix(brca.count[,codes.type]),
                                 colData = brca.clini[codes.type,],
                                 design = as.formula(~ sample.type))

ddsObj <- DESeq(ddsObj)

res.shr <- results(ddsObj)

summary(res.shr)

dea.NT.TP <- as.data.frame(res.shr) %>%
  rownames_to_column("symbol") %>% 
  left_join(brca.annot, by="symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseq <- dea.NT.TP  %>% filter(abs(logFC) >= 3, pvalue <= 0.05)
dim(df.deseq)
# [1] 1131    9
genes.DEA.NT.vs.TP.lst <- unique(df.deseq$symbol)
write(genes.DEA.NT.vs.TP.lst,  file = "genes.DEA.NT.vs.TP.lst")

save(dea.NT.TP, file = "~/BRCA-TCGAA/BRCA-TCGAA/Data/dea.NT.TP.rda", compress = T)

cutoff <- sort(dea.NT.TP$pvalue)[10]
shrink.deseq.cut <- dea.NT.TP %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, symbol, ""))

ggplot(shrink.deseq.cut, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=0.5) +
  labs(x="mean of normalised counts", y="log Fold Change") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

cutoff <- sort(dea.NT.TP$pvalue)[10]
shrink.deseq.cut <- dea.NT.TP %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, symbol, ""))

ggplot(shrink.deseq.cut, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

res.df <- as.data.frame(res.shr)
res.df$log10MeanNormCount <- log10(res.df$baseMean + 1)
idx <-(rowSums(counts(ddsObj)) > 0)
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1

status <- as.numeric(res.df$padj < 0.05)

glMDPlot(res.df[idx,],
         xval="baseMean",
         yval="log2FoldChange",
         counts=counts(ddsObj)[idx,],
         display.columns=c("GeneID"),
         anno=data.frame(GeneID=rownames(ddsObj)[idx]),
         groups = brca.clini[codes.type, "sample.type"],
         side.xlab = "Group",
         side.ylab = "Expression (log2)",         
         samples = brca.clini[codes.type, "codes"],
         status=status,
         folder = "MDPlot_BRCA.NT.vs.TP",
         html = "index",
         launch=TRUE)

filtTab.deseq <- dea.NT.TP %>%
  filter(!is.na(FDR)) %>%
  mutate(`-log10(FDR)` = -log10(FDR))

filtTab.deseq <- filtTab.deseq  %>%
  mutate(`-log10(FDR)`=pmin(`-log10(FDR)`))

filtTab.deseq <- filtTab.deseq[!duplicated(filtTab.deseq$symbol), ]
rownames(filtTab.deseq) <- filtTab.deseq$symbol  
de <- as.integer(abs(filtTab.deseq$logFC) >= 3 & filtTab.deseq$FDR <= 0.01)

glXYPlot(
  x = filtTab.deseq$logFC,
  y = -log10(filtTab.deseq$pvalue),
  xlab = "logFC",
  ylab = "-log10pvalue (FDR)",
  main = "NT.vs.TP",
  counts = log2( counts(ddsObj)[filtTab.deseq$symbol,] ),
  groups = brca.clini[codes.type, "sample.type"],
  side.xlab = "Group",
  side.ylab = "Expression (log2)",
  samples = brca.clini[codes.type, "codes"],
  status = de,
  side.main = "symbol",
  display.columns = c("symbol", "logFC", "FDR", "ensgene", "description"),
  anno = filtTab.deseq,
  folder = "XYPlot_BRCA.NT.vs.TP",
  html = "index",
  launch = F)

# Running DESeq2 Differential Expression for race = white

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

