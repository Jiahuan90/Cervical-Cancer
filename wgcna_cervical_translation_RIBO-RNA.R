#remove.packages("WGCNA")
#BiocManager::install("WGCNA")
library(dplyr)
library(WGCNA)

options(stringsAsFactors = FALSE)  #### important, don't omit

setwd("I:\\博士后\\cervical_cancer_translatome\\Normal_VS_cervival_cancer\\Fig2\\6_WGCNA分析重做\\27_重做cancer部分WGCNA_top3000_RNA-TE-RIBO/")

#================================================================
# 1.Data input, cleaning, preprocesing
#================================================================
###### 1.1 RIBO level top5000 MAD input
##### 这个文件，是3个层次中每个层面的top3000
Ribo_data <- read.table(file = "RIBOlevel_MADgene.txt", header = T, row.names = 1)
Ribo_data <- log2(Ribo_data + 1)
t_Ribo_data <- t(Ribo_data) 
Ribo_data %>% head()
Ribo_data %>% nrow()

###### 1.2 RNA level correspond 5000 genes in RPF level
RNA_data <- read.table("RNAlevel_MADgene.txt", header = T, row.names = 1)
RNA_data <- log2(RNA_data + 1)
t_RNA_data <- t(RNA_data)


nSets <- 2
Set_labels <- c("RIBO", "RNA")
Short_labels <- c("RIBO", "RNA")


### Form multi-set expression data: columns starting from 9 contain actual expression data.
Multi_expr <- vector(mode = "list", length = nSets)

############# MultiExpr list 中的第一个元素也是一个name = data的list
Multi_expr[[1]] <- list(data = as.data.frame(t_Ribo_data)) 


Multi_expr[[2]] <- list(data = as.data.frame(t_RNA_data))



dim(Multi_expr[[1]]$data)
checkSets(Multi_expr)

Expr_size <- checkSets(Multi_expr)


gsg <- goodSamplesGenesMS(Multi_expr, verbose = 3);
gsg$allOK

if (!gsg$allOK) {
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(Multi_expr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:Expr_size$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(Multi_expr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    Multi_expr[[set]]$data <- Multi_expr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  } #
  #Update exprSize
  Expr_size <- checkSets(Multi_expr)
}




#### we now cluster the samples on their Euclidean distance, separately in each set
SampleTrees <- list()
for (set in 1:nSets) {
  #### 构建了一个SampleTree的list用于储存结果
  SampleTrees[[set]] <- hclust(dist(Multi_expr[[set]]$data), method = "average")
}


pdf(file = "pic_Sample_clustering_RIBO-RNA.pdf", width = 12, height = 12)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(SampleTrees[[set]], main = paste("Sample clustering on all genes in", Set_labels[set]),
       xlab="", sub="", cex = 0.7)
dev.off()


# Define data set dimensions
Expr_size$nGenes
Expr_size$nSamples
nGenes <- 6161
nSamples <- c(25, 25)

#=====================================================================================
# 2. 构建共表达网络 + module的检测
#=====================================================================================
#lnames <- load(file = "Consensus-dataInput_RIBO-RNA.RData")

##### 2.1 选择一个合适的β
nSets <- checkSets(Multi_expr)$nSets
Powers <- c(seq(4,10,by=1), seq(12,20, by=2))


##### 构建了一个包含两个元素的list
Power_tables <- vector(mode = "list", length = nSets)

# Call the network topology analysis function for each set in turn
for (set in 1:nSets) {
  Power_tables[[set]] <- list(data = pickSoftThreshold(Multi_expr[[set]]$data, 
                                                       powerVector = Powers,
                                                       networkType = "signed hybrid",
                                                       corFnc = bicor,
                                                       corOptions = list(maxPOutliers = 0.05),
                                                       verbose = 5)[[2]])
}

set <- 1
X <- Power_tables[[set]]$data[,1]
Y <- -sign(Power_tables[[set]]$data[,3])*Power_tables[[set]]$data[,2]
Z <- rep("RIBO", time = length(X))

Soft_power_ribo <- cbind(X, Y, Z)

set <- 2
X <- Power_tables[[set]]$data[,1]
Y <- -sign(Power_tables[[set]]$data[,3])*Power_tables[[set]]$data[,2]
Z <- rep("RNA", time = length(X))

X <- as.numeric(X)
Y <- as.numeric(Y)

Soft_power_rna <- cbind(X, Y, Z)
Soft_power_rna



Table <- rbind(Soft_power_ribo, Soft_power_rna)
#Table <- rbind(Table, Soft_power_te)
Table <- as.data.frame(Table)
Table$Y <- as.numeric(Table$Y)

Table
library(ggplot2)
p <- ggplot(Table, aes(x = X, y=Y, colour = Z)) + 
  geom_point(size = 0.2) + 
  scale_x_discrete(limits = c("4", "5", "6", "7", "8", "9", 
                              "10", "12", "14", "16", "18", "20")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = "solid", color = "red", size = 0.4) + 
  xlab("Soft power") + ylab("Scale free topology") +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.9, vjust = 0.9),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.line = element_line(size = 0.3, colour = "black"),
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.12, "cm"),
    legend.key.height = unit(0.14, "cm"),
    legend.key.width = unit(0.36, "cm"),
    legend.text = element_text(colour = "black", size = 8),
    legend.title = element_text(colour = "black", size = 8),
    plot.title = element_text(hjust = 0.5, size = 7),
    plot.background = element_rect(fill="transparent",colour=NA),
    panel.background = element_rect(fill="transparent",colour=NA),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.border = element_blank()
  )
  

p
ggsave("pic_soft_power_RIBO-RNA-TE.pdf", width = 3, height = 2, device = "pdf")

####### 2.2 计算Adjancies + Topological overlap
Soft_power <- c(7, 5)
# Initialize an appropriate array to hold the adjacencies
Adjacencies <- array(0, dim = c(nSets, nGenes, nGenes))

Multi_expr[[2]]$data
# Calculate adjacencies in each individual data set
for (set in 1:nSets){
  Adjacencies[set, , ] <- adjacency(Multi_expr[[set]]$data, power = Soft_power[set],
                                    type = "signed hybrid",
                                    corFnc = "bicor",
                                    corOptions = list(maxPOutliers = 0.05))
}


# Initialize an appropriate array to hold the TOMs
TOM <- array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set

for (set in 1:nSets){
  TOM[set, , ] <- TOMsimilarity(Adjacencies[set, , ], TOMType = "signed")
}

TOM[1, ,] %>% nrow()
TOM_ribo_level <- TOM[1, ,]
TOM_ribo_level[1:5, 1:5]

rownames(TOM_ribo_level) <- colnames(Multi_expr[[1]]$data)
colnames(TOM_ribo_level) <- colnames(Multi_expr[[1]]$data)


####### 2.3 对两个TOM进行scale，使得它们可比
scaleP <- 0.95
set.seed(12345)

# Sample sufficiently large number of TOM entries
nSamples <- as.integer(1/(1-scaleP) * 1000)

# Choose the sampled TOM entries
scaleSample <- sample(nGenes*(nGenes-1)/2, size = nSamples) ## 对nGens = 5000 进行随机抽样，抽样数量 = nSample


TOMScalingSamples <- list()
# These are TOM values at reference percentile
scaleQuant <- rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers <- rep(1, nSets)

# Loop over sets
for (set in 1:nSets) {
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] <- as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] <- quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1) ### ----------------- 取第二个list,这个取决于我把谁放到第二位
  {
    scalePowers[set] <- log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] <- TOM[set, ,]^scalePowers[set];
  }
}

####### 2.4 calculate consensus topological overlap
Consensus_TOM <- pmin(TOM[1, , ], TOM[2, , ])
Consensus_tree <- hclust(as.dist(1 - Consensus_TOM), method = "average")


Module_size <- 89
Module <- cutreeDynamic(dendro = Consensus_tree, 
                        distM = 1-Consensus_TOM,
                        deepSplit = 2, 
                        cutHeight = 0.995,
                        minClusterSize = Module_size,
                        pamRespectsDendro = FALSE )

Module_color <- labels2colors(Module)

####### 2.5 picture clustering (optional)
sizeGrWindow(8,6)
plotDendroAndColors(Consensus_tree, Module, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

###### 2.6 calculate module eigengenes
Eigengens <- multiSetMEs(Multi_expr, colors = NULL, universalColors = Module_color)

consMEDiss <- consensusMEDissimilarity(Eigengens)

Consensus_METree <- hclust(as.dist(consMEDiss), method = "average")

Merged_consensus_module <- mergeCloseModules(Multi_expr, Module_color, 
                                             cutHeight = 0.1, 
                                             verbose = 3, 
                                             corFnc = "bicor",
                                             corOptions = list(maxPOutliers = 0.05))

Gene_colors <- Merged_consensus_module$colors %>% as.data.frame()

rownames(Gene_colors) <- colnames(Multi_expr[[1]]$data)

# Convert labels to colors
Merged_consensus_module_colors <- labels2colors(Merged_consensus_module_labels) # ------------ is a matrix 
# Eigengenes of the new merged modules:
Consensus_MEs <- Merged_consensus_module$newMEs


save(Merged_consensus_module, Merged_consensus_module_labels, Consensus_METree, file = "Consensus-NetworkConstruction_RIBO-RNA.RData")

RIBO_colors <- Merged_consensus_module$newMEs[[1]]$validColors

####################################
#### 3. module 和临床特征关联
####################################
### 3.1 读取临床数据
PhenoData <- read.table("Table S5", header = T, row.names = 1)
PhenoData


Ribo_PhenoData <- PhenoData 

RNA_PhenoData <- PhenoData

Traits <- vector(mode="list", length = nSets)

Traits[[1]] <- list(data = Ribo_PhenoData)

Traits[[2]] <- list(data = RNA_PhenoData)

### 3.2 计算correlation
Module_traitCor <- list()
Module_trait_pvalue <- list()


for (set in 1:nSets){
  Module_traitCor[[set]] <- bicorAndPvalue(Merged_consensus_module$newMEs[[set]]$data, Traits[[set]]$data, robustY = FALSE, maxPOutliers = 0.05)$bicor
  Module_trait_pvalue[[set]] <- bicorAndPvalue(Merged_consensus_module$newMEs[[set]]$data, Traits[[set]]$data, robustY = FALSE, maxPOutliers = 0.05)$p
}


MEColors <- substring(names(Merged_consensus_module$newMEs[[1]]$data), 3)
MEColorNames <- paste("ME", MEColors, sep="")

##### 3.3 draw picture Fig 6A
# Plot the module-trait relationship table for set number 1
set <- 1
TextMatrix <- paste(signif(Module_traitCor[[set]], 2), "\n(",
                   signif(Module_trait_pvalue[[set]], 1), ")", sep = "")

dim(TextMatrix) <- dim(Module_traitCor[[set]])

pdf(file = "new_ModuleTraitRelationships-RIBO.pdf", width = 4, height = 5)

labeledHeatmap(Matrix = Module_traitCor[[set]],
               xLabelsAngle = 60,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = TextMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", Set_labels[set]))
dev.off()

# Plot the module-trait relationship table for set number 2
set <- 2
TextMatrix <- paste(signif(Module_traitCor[[set]], 2), "\n(",
                    signif(Module_trait_pvalue[[set]], 1), ")", sep = "")

dim(TextMatrix) <- dim(Module_traitCor[[set]])

pdf(file = "new_ModuleTraitRelationships-RNA.pdf", width = 4, height = 5)
labeledHeatmap(Matrix = Module_traitCor[[set]],
               xLabelsAngle = 60, 
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = TextMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", Set_labels[set]))
dev.off()



