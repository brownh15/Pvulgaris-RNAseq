rm(list=ls())
setwd('/Users/hannahcushman/Desktop/Niederhuth_Lab/Bean_Microbiome /Bean_counts')
getwd()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library("edgeR")

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

sample_info <- c("Mock_control", "Mock_control", "Mock_control", "Zeb_control", "Zeb_control", "Zeb_control", "Mock_drought", "Mock_drought", "Mock_drought", "Zeb_drought", "Zeb_drought", "Zeb_drought", "Mock_nutrient", "Mock_nutrient", "Mock_nutrient", "Zeb_nutrient", "Zeb_nutrient", "Zeb_nutrient")
dge <- DGEList(counts = df36, group = factor(sample_info))
dge
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(object = dge)
dge <- estimateDisp(y = dge)
dge <- estimateTagwiseDisp(y = dge)
dge

MC_ZC <- exactTest(object = dge, pair=c("Mock_control","Zeb_control"))
MC_ZC
summary(decideTests(object = MC_ZC, lfc = 1))

MC_MD <- exactTest(object = dge, pair=c("Mock_control","Mock_drought"))
MC_MD
summary(decideTests(object = MC_MD, lfc = 1))

MC_MN <- exactTest(object = dge, pair=c("Mock_control","Mock_nutrient"))
MC_MN
summary(decideTests(object = MC_MN, lfc = 1))

ZC_ZD <- exactTest(object = dge, pair=c("Zeb_control","Zeb_drought"))
ZC_ZD
summary(decideTests(object = ZC_ZD, lfc = 1))

ZC_ZN <- exactTest(object = dge, pair=c("Zeb_control","Zeb_nutrient"))
ZC_ZN
summary(decideTests(object = ZC_ZN, lfc = 1))

MN_ZN <- exactTest(object = dge, pair=c("Mock_nutrient","Zeb_nutrient"))
MN_ZN
summary(decideTests(object = MN_ZN, lfc = 1))

MD_ZD <- exactTest(object = dge, pair=c("Mock_drought","Zeb_drought"))
MD_ZD
summary(decideTests(object = MD_ZD, lfc = 1))



de2 <- decideTestsDGE(MD_ZD, adjust.method="BH", p.value = 0.05)
summary(de2)
de2tags12 <- rownames(dge)[as.logical(de2)]
de2tags12


plotSmear(dge, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")



plotMDS(dge, method="bcv", col=as.numeric(dge$samples$group))
legend("bottomleft", as.character(unique(dge$samples$group)), col=1:3, pch=20)

names(dge)
d1 <- estimateTagwiseDisp(dge)
plotBCV(d1)






