# prepare data
data = readRDS("../Summary202301/savedRDS/intersectProt_gene.RDS")
pheno = readRDS("../Summary202301/savedRDS/pheno.RDS")
spine.dens = readRDS("../Summary202301/savedRDS/spine.RDS")
spine.dens = spine.dens[,colnames(data)]

FV_binary = pheno[colnames(data),]$Group
FV_continuous = as.numeric(spine.dens[1,])
metadata = pheno[colnames(data),c("Sex", "pH")]

# functions
source("HelperFuncs.R")


# feature selection
source("FeatureSelection.R")

# Feature selection combining limma and elastic net regression

FV = FV_continuous

# this returns a list for limma and 100 elastic_net
# feature selection with binomial FV
selectionRes = featureSelection(data = data, focalVariable = FV_binary, class = "binomial", contrast = "Sz-C", limmaCutoff = 20, glmnetCutoff = 20)
plotSelection(data, selectionRes)

# feature selection with continuous FV
selectionRes.spineL3 = featureSelection(data = data, focalVariable = FV_continuous, class = "gaussian", limmaCutoff = 17, glmnetCutoff = 23)
plotSelection(data, selectionRes.spineL3)


# simpleDcomp to get Latent variables
source("LatentVariableAnalysis.R")

prot.data = matrix(unlist(data),ncol = ncol(data))
rownames(prot.data) = rownames(data)
colnames(prot.data) = colnames(data)

decompRes = simpleDecomp(tscale(prot.data), k=20, adaptive.frac = 0.01)

# look at decomp results
# we should see some with very low variance (drop out)
boxplot(t(decompRes$B), las=2)
plot(sort(bvar<-apply(decompRes$B, 1, var)), log="y")
# 1e-2 is a good cutoff
mypheatmap(getTopEachColumn(decompRes$Z,top = 50), main="top50", show_rownames = F)
mypheatmap(getTopEachColumn(decompRes$Z[, bvar > 1e-2], top = 50), main = "Latent Variables with bvar>1e-2, top50", show_rownames = F)

# Get prediction (continuous) variable from data
predictL3 = getPrediction(data, spine.dens[1,])
predictScatter(real = spine.dens[1,], predict = predictL3$predictionRes)


# prepare data for Causal input
prot.use = selectionRes.spineL3$selected
lv.use = colnames(decompRes$Z[,bvar>1e-2])[colSums(decompRes$Z[prot.use,bvar>1e-2])>0]

metadata$Sex = metadata$Sex == "F"
metadata$pH = (metadata$pH-mean(metadata$pH))/sd(metadata$pH)
metadata = data.frame(metadata)

tscaled_prediction = (predictL3$predictionRes - mean(predictL3$predictionRes)) / sd(predictL3$predictionRes)

# bootstrapping
source("NetworkConstruction.R")
source("NotearsImplementation.R")
FocalL3 = graphSampling(data = tscale(data), LV = tscale(decompRes$B[lv.use,]),
						FV = tscaled_prediction, data.use = prot.use, group = pheno[colnames(data),]$Group=="Sz",
						name = c(lv.use, "Spine", "Sz"))

FocalL3_meta = graphSampling(data = tscale(data), LV = tscale(decompRes$B[lv.use,]),
						FV = tscaled_prediction, data.use = prot.use, group = pheno[colnames(data),]$Group=="Sz",
						metadata = metadata$pH, name = c(lv.use, "Spine", "Sz", colnames(metadata)))


# construction of summary adjacency matrix
name = c(prot.use, lv.use, "Spine", "Sz")
GraphL3 = graphConstruction(graph_list = FocalL3, name = name)

name = c(prot.use, lv.use, "Spine", "Sz", "pH")
GraphL3_meta = graphConstruction(graph_list = FocalL3_meta, name = name)

# visualization of causal network
source("NetworkVisualization.R")
# gvarShape and gvarType (type and subtype)

gvarShape = c(rep("Protein", length(prot.use)), rep("LV", length(lv.use)), rep("FV", 2))

# graph varType
Prot = prot.use
Prot[grep("p", Prot)] = "Phospho"
Prot[grep(".syn", Prot)] = "Synaptosome"
Prot[Prot != "Phospho" & Prot != "Synaptosome"] = "Homogenate"
lenLV = rep("LV",length(lv.use))
gvarType = c(Prot, lenLV, "Spine", "Diagnosis")


network_visualize(GraphL3, gvarType, gvarShape) + 
	scale_shape_manual(values=c(15, 16, 17, 18), breaks = c("Genetics", "Protein", "LV", "FV", "Metadata")) +
	scale_color_manual(breaks = c("Phospho", "Homogenate", "Synaptosome", "LV", "Spine", "Diagnosis"),
					   values = c("#a9d18e", "#fb9288", "#9dc3e6", "#ffd966", "#8ccebc", "#d4d4d4"))










as.logical.factor()
XPO7 = data[grep("XPO7",rownames(data)),]

AKT3 = data[grep("AKT3",rownames(data)),]
FAM = data[grep("FAM114A2",rownames(data)),]
DDHD2 = data[grep("DDHD2",rownames(data)),]
CLASP1 = data[grep("CLASP1",rownames(data)),]
R3HDM2 = data[grep("R3HDM2",rownames(data)),]
library(reshape2)
library(ggpubr)
df = melt(as.matrix(XPO7))
df$pheno = pheno[as.character(df$Var2),]$Group
colnames(df) = c("protein", "patient", "value", "pheno")
df$spine = rep(FV_continuous, each = 2)
ggboxplot(df, x = "pheno", y = "value", facet.by = "protein")+
	stat_compare_means(method = "t.test", label.x = 1, label.y = 2.5) +
	stat_compare_means(aes(label = ..p.signif..), method = "t.test", label.y = 1.5, label.x = 1.4)
ggscatter(df, x = "spine", y = "value", facet.by = "protein",
		  add = "reg.line", conf.int = TRUE, # Add confidence interval
		  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
		  cor.coeff.args = list(method = "pearson", label.sep = "\n", size=3),)

