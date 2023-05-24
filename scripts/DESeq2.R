rm(list=ls())
setwd('/Users/hannahcushman/Desktop/Niederhuth_Lab/Bean_Microbiome /Bean_counts')
getwd()

install.packages("BiocManager")
library(BiocManager)

BiocManager::install("DESeq2")
library("DESeq2")

df1 <- read.delim("E22_counts.tsv")
df2 <- read.delim("E23_counts.tsv")
df3 <- read.delim("E39_counts.tsv")
df4 <- read.delim("E11_counts.tsv")
df5 <- read.delim("E47_counts.tsv")
df6 <- read.delim("E64_counts.tsv")
df7 <- read.delim("E4_counts.tsv")
df8 <- read.delim("E28_counts.tsv")
df9 <- read.delim("E40_counts.tsv")
df10 <- read.delim("E15_counts.tsv")
df11 <- read.delim("E26_counts.tsv")
df12 <- read.delim("E60_counts.tsv")
df13 <- read.delim("E5_counts.tsv")
df14 <- read.delim("E20_counts.tsv")
df15 <- read.delim("E42_counts.tsv")
df16 <- read.delim("E12_counts.tsv")
df17 <- read.delim("E13_counts.tsv")
df18 <- read.delim("E19_counts.tsv")

colnames(df1) <- c("Gene","Mock_control1")
print("Dataframe with manually assigned Header")

colnames(df2) <- c("Gene","Mock_control2")
print("Dataframe with manually assigned Header")

colnames(df3) <- c("Gene","Mock_control3")
print("Dataframe with manually assigned Header")

colnames(df4) <- c("Gene","Zeb_control1")
print("Dataframe with manually assigned Header")

colnames(df5) <- c("Gene","Zeb_control2")
print("Dataframe with manually assigned Header")

colnames(df6) <- c("Gene","Zeb_control3")
print("Dataframe with manually assigned Header")

colnames(df7) <- c("Gene","Mock_drought1")
print("Dataframe with manually assigned Header")

colnames(df8) <- c("Gene","Mock_drought2")
print("Dataframe with manually assigned Header")

colnames(df9) <- c("Gene","Mock_drought3")
print("Dataframe with manually assigned Header")

colnames(df10) <- c("Gene","Zeb_drought1")
print("Dataframe with manually assigned Header")

colnames(df11) <- c("Gene","Zeb_drought2")
print("Dataframe with manually assigned Header")

colnames(df12) <- c("Gene","Zeb_drought3")
print("Dataframe with manually assigned Header")

colnames(df13) <- c("Gene","Mock_nutrient1")
print("Dataframe with manually assigned Header")

colnames(df14) <- c("Gene","Mock_nutrient2")
print("Dataframe with manually assigned Header")

colnames(df15) <- c("Gene","Mock_nutrient3")
print("Dataframe with manually assigned Header")

colnames(df16) <- c("Gene","Zeb_nutrient1")
print("Dataframe with manually assigned Header")

colnames(df17) <- c("Gene","Zeb_nutrient2")
print("Dataframe with manually assigned Header")

colnames(df18) <- c("Gene","Zeb_nutrient3")
print("Dataframe with manually assigned Header")
df18


df19 <- merge(df1, df2, by.x="Gene", by.y="Gene")
df20 <- merge(df19, df3, by.x="Gene", by.y="Gene")
df21 <- merge(df20, df4, by.x="Gene", by.y="Gene")
df22 <- merge(df21, df5, by.x="Gene", by.y="Gene")
df23 <- merge(df22, df6, by.x="Gene", by.y="Gene")
df24 <- merge(df23, df7, by.x="Gene", by.y="Gene")
df25 <- merge(df24, df8, by.x="Gene", by.y="Gene")
df26 <- merge(df25, df9, by.x="Gene", by.y="Gene")
df27 <- merge(df26, df10, by.x="Gene", by.y="Gene")
df28 <- merge(df27, df11, by.x="Gene", by.y="Gene")
df29 <- merge(df28, df12, by.x="Gene", by.y="Gene")
df30 <- merge(df29, df13, by.x="Gene", by.y="Gene")
df31 <- merge(df30, df14, by.x="Gene", by.y="Gene")
df32 <- merge(df31, df15, by.x="Gene", by.y="Gene")
df33 <- merge(df32, df16, by.x="Gene", by.y="Gene")
df34 <- merge(df33, df17, by.x="Gene", by.y="Gene")
df35 <- merge(df34, df18, by.x="Gene", by.y="Gene")



df35$Mock_control1 <- as.integer(df35$Mock_control1)
df35$Mock_drought1 <- as.integer(df35$Mock_drought1)
df35$Mock_nutrient1 <- as.integer(df35$Mock_nutrient1)
df35$Zeb_control1 <- as.integer(df35$Zeb_control1)
df35$Zeb_drought1 <- as.integer(df35$Zeb_drought1)
df35$Zeb_nutrient1 <- as.integer(df35$Zeb_nutrient1)
df35$Mock_control2 <- as.integer(df35$Mock_control2)
df35$Mock_drought2 <- as.integer(df35$Mock_drought2)
df35$Mock_nutrient2 <- as.integer(df35$Mock_nutrient2)
df35$Zeb_control2 <- as.integer(df35$Zeb_control2)
df35$Zeb_drought2 <- as.integer(df35$Zeb_drought2)
df35$Zeb_nutrient2 <- as.integer(df35$Zeb_nutrient2)
df35$Mock_control3 <- as.integer(df35$Mock_control3)
df35$Mock_drought3 <- as.integer(df35$Mock_drought3)
df35$Mock_nutrient3 <- as.integer(df35$Mock_nutrient3)
df35$Zeb_control3 <- as.integer(df35$Zeb_control3)
df35$Zeb_drought3 <- as.integer(df35$Zeb_drought3)
df35$Zeb_nutrient3 <- as.integer(df35$Zeb_nutrient3)


rownames(df35) <-df35$Gene
df36 <- df35[-c(1)]
colnames(df36) <- c("Mock_control1", "Mock_control2", "Mock_control3", "Zeb_control1", "Zeb_control2", "Zeb_control3", "Mock_drought1", "Mock_drought2", "Mock_drought3", "Zeb_drought1", "Zeb_drought2", "Zeb_drought3", "Mock_nutrient1", "Mock_nutrient2", "Mock_nutrient3", "Zeb_nutrient1", "Zeb_nutrient2", "Zeb_nutrient3")

coldata <- data.frame(
  sample = c( "Mock_control1", "Mock_control2", "Mock_control3", "Zeb_control1", "Zeb_control2", "Zeb_control3", "Mock_drought1", "Mock_drought2", "Mock_drought3", "Zeb_drought1", "Zeb_drought2", "Zeb_drought3", "Mock_nutrient1", "Mock_nutrient2", "Mock_nutrient3", "Zeb_nutrient1", "Zeb_nutrient2", "Zeb_nutrient3"),
  condition = c( "Mock_control", "Mock_drought", "Mock_nutrient", "Zeb_control", "Zeb_drought", "Zeb_nutrient" ), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)


count_matrix <- as.matrix(df36)
all(rownames(coldata) %in% colnames(count_matrix))
all(rownames(coldata) == colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, 
                              design = ~ condition)

dim(dds)

dim(dds[rowSums(counts(dds)) > 5, ])

deseq2Data <- dds[rowSums(counts(dds)) > 5, ]

BiocManager::install("BiocParallel")
library(BiocParallel)
register(MulticoreParam(4))

deseq2Data <- DESeq(deseq2Data)

deseq2Results <- results(deseq2Data, contrast=c("condition", "Mock_control", "Zeb_control"))
deseq2Results
summary(deseq2Results)

otop2Counts <- plotCounts(deseq2Data, gene="Phvul.009G134100.v2.1", intgroup=c("condition"), returnData=TRUE)

colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
ggplot(otop2Counts, aes(x=condition, y=count, colour=condition, group=condition)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Phvul.009G134100.v2.1")


colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
p<-ggplot(data=otop2Counts, aes(x=condition, y=count, colour=condition, group=condition)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=15, hjust=1)) +
  scale_colour_manual(values=colourPallette) +
  guides(colour=guide_legend(ncol=3))

p


p<-ggplot(otop2Counts, aes(x=condition, y=count, fill=condition)) +
  geom_bar(stat="identity")+theme_minimal()+ggtitle("Phvul.009G134100.v2.1")
p

plotMA(deseq2Results)

deseq2ResDF <- as.data.frame(deseq2Results)
head(deseq2ResDF)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)

deseq2VST <- vst(deseq2Data)

deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
library(reshape2)

deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

install.packages("ggdendro")
library(ggdendro)

sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

library(gtable)
library(grid)

sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

sampleDendrogramGrob$widths
heatmapGrob$widths

sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))
grid.draw(finalGrob)



