library(dplyr)
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
library(ggpubr)
library(ggsignif)
library(rstatix)
setwd("I:ssGSEA")
############################################################################################################
########################## RIBO level
###########################################################################################################
a <- read.table('geneSym_TPM_RPF_normal_10_VS_cancer_25.txt', header = T, row.names=1)
names(a) ##### get column names
a <- a[,!names(a) %in% c("cancer_sample_21", "cancer_sample_22")] #### remove 2 adenocarcinoma

#### filtering data
Data <- a[!apply(a,1,sum)==0,]

Data <- Data[!duplicated(rownames(Data)), ]
Log_data <- log2(Data + 1)

## load data
Gene_set <- read.csv("immune_signature.csv")[,1:2]
tail(Gene_set)
List <- split(as.matrix(Gene_set)[,1], Gene_set[,2])

GSVA_matrix <- gsva(as.matrix(Data), List, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

GSVA_matrix %>% head()

GSVA_matrix %>% rownames()
#write.table(GSVA_matrix, file = "TE_immune.txt", quote = F, sep = "\t", col.names = T)
#library(pheatmap)
GSVA_matrix1 <- t(scale(t(GSVA_matrix)))
GSVA_matrix1[GSVA_matrix1 < -2] <- -2
GSVA_matrix1[GSVA_matrix1 > 2] <- 2

Anti_tumor <- c('Activated CD8 T cell', 'Type 1 T helper cell', 'Macrophage_M1', 'Natural killer cell', 'Gamma delta T cell', 'new_DC') ## 6

Anti_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Anti_tumor),]

Pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'Macrophage_M2', 'MDSC', "Type 17 T helper cell") ## 5
Pro_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Pro_tumor),]

Pro_mat %>% rownames()

Unknown <- c("B_cell", "Eosinophil", "Mast cell", "Monocyte", "Neutrophil") ### 5
Unknown_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Unknown),]


GSVA_matrix2 <- rbind(Anti_mat, Pro_mat, Unknown_mat)

normalization <- function(x){
                             return((x-min(x))/(max(x)-min(x)))}

RPF_Nor_gsva_matrix1 <- normalization(GSVA_matrix2)
RPF_Nor_gsva_matrix1 %>% rownames()

######################################################################################
######################### TE level
######################################################################################
a <- read.table('geneSym_TE_normal_10_VS_cancer_25.txt', header = T, row.names=1)
names(a) ##### get column names
a <- a[,!names(a) %in% c("cancer_sample_21", "cancer_sample_22")] #### remove 2 adenocarcinoma

#### filtering data
Data <- a[!apply(a,1,sum)==0,]

Data <- Data[!duplicated(rownames(Data)), ]
Log_data <- log2(Data + 1)

## load data
Gene_set <- read.csv("immune_signature.csv")[,1:2]
tail(Gene_set)
List <- split(as.matrix(Gene_set)[,1], Gene_set[,2])

GSVA_matrix <- gsva(as.matrix(Data), List, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

GSVA_matrix1 <- t(scale(t(GSVA_matrix)))
GSVA_matrix1[GSVA_matrix1 < -2] <- -2
GSVA_matrix1[GSVA_matrix1 > 2] <- 2

Anti_tumor <- c('Activated CD8 T cell', 'Type 1 T helper cell', 'Macrophage_M1', 'Natural killer cell', 'Gamma delta T cell', 'new_DC') ## 6

Anti_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Anti_tumor),]

Pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'Macrophage_M2', 'MDSC', "Type 17 T helper cell") ## 5
Pro_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Pro_tumor),]

Pro_mat %>% rownames()

Unknown <- c("B_cell", "Eosinophil", "Mast cell", "Monocyte", "Neutrophil") ### 5
Unknown_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Unknown),]


GSVA_matrix2 <- rbind(Anti_mat, Pro_mat, Unknown_mat)

normalization <- function(x){
  return((x-min(x))/(max(x)-min(x)))}

TE_Nor_gsva_matrix1 <- normalization(GSVA_matrix2)
TE_Nor_gsva_matrix1 %>% rownames()


##################################################################################################
########################### RNA level
#################################################################################################
a <- read.table('geneSym_RNA_normal_10_VS_cancer_25.txt', header = T, row.names=1)
names(a) ##### get column names
a <- a[,!names(a) %in% c("cancer_sample_21", "cancer_sample_22")] #### remove 2 adenocarcinoma

#### filtering data
Data <- a[!apply(a,1,sum)==0,]

Data <- Data[!duplicated(rownames(Data)), ]
Log_data <- log2(Data + 1)

## load data
Gene_set <- read.csv("immune_signature.csv")[,1:2]
tail(Gene_set)
List <- split(as.matrix(Gene_set)[,1], Gene_set[,2])

GSVA_matrix <- gsva(as.matrix(Data), List, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

GSVA_matrix1 <- t(scale(t(GSVA_matrix)))
GSVA_matrix1[GSVA_matrix1 < -2] <- -2
GSVA_matrix1[GSVA_matrix1 > 2] <- 2

Anti_tumor <- c('Activated CD8 T cell', 'Type 1 T helper cell', 'Macrophage_M1', 'Natural killer cell', 'Gamma delta T cell', 'new_DC') ## 6

Anti_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Anti_tumor),]

Pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'Macrophage_M2', 'MDSC', "Type 17 T helper cell") ## 5
Pro_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Pro_tumor),]

Pro_mat %>% rownames()

Unknown <- c("B_cell", "Eosinophil", "Mast cell", "Monocyte", "Neutrophil") ### 5
Unknown_mat <- GSVA_matrix1[which(rownames(GSVA_matrix1) %in% Unknown),]


GSVA_matrix2 <- rbind(Anti_mat, Pro_mat, Unknown_mat)

normalization <- function(x){
  return((x-min(x))/(max(x)-min(x)))}

RNA_Nor_gsva_matrix1 <- normalization(GSVA_matrix2)


############### calculate Total score
### RNA level
rm(Colsum_mat_RNA)
Colsum_mat_RNA <- RNA_Nor_gsva_matrix1[which(rownames(RNA_Nor_gsva_matrix1) %in% Unknown),] %>% colSums() %>% as.data.frame()
Type <- c(rep("normal", times = 10), rep("cancer", times = 23)) %>% as.data.frame()
Type1 <- "RNA"
Colsum_mat_RNA <- cbind(Colsum_mat_RNA, Type)
Colsum_mat_RNA <- cbind(Colsum_mat_RNA, Type1)
Colsum_mat_RNA %>% head()
 
### TE level
rm(Colsum_mat_TE)
Colsum_mat_TE <- TE_Nor_gsva_matrix1[which(rownames(TE_Nor_gsva_matrix1) %in% Unknown),] %>% colSums() %>% as.data.frame()
Type <- c(rep("normal", times = 10), rep("cancer", times = 23)) %>% as.data.frame()
Type1 <- "TE"
Colsum_mat_TE <- cbind(Colsum_mat_TE, Type)
Colsum_mat_TE <- cbind(Colsum_mat_TE, Type1)
Colsum_mat_TE %>% head()
 
### RPF level
rm(Colsum_mat_RPF)
Colsum_mat_RPF <- RPF_Nor_gsva_matrix1[which(rownames(RPF_Nor_gsva_matrix1) %in% Unknown),] %>% colSums() %>% as.data.frame()
Type <- c(rep("normal", times = 10), rep("cancer", times = 23)) %>% as.data.frame()
Type1 <- "RIBO"
Colsum_mat_RPF <- cbind(Colsum_mat_RPF, Type)
Colsum_mat_RPF <- cbind(Colsum_mat_RPF, Type1)
Colsum_mat_RPF %>% head()



