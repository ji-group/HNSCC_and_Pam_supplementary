#############################Human Palmitoylation-related Gene (PRG) Set########################

rm(list=ls())
library(msigdbr)
main.dir <- getwd()
sub.dir <- "./PAM"
dir.create(file.path(main.dir, sub.dir), showWarnings = FALSE)
setwd(file.path(main.dir, sub.dir))

# Load all MSigDB gene sets for Homo sapiens
msigdb_data <- msigdbr(species = "Homo sapiens")

# Find palmitoylation-related gene sets
palmitoyl_genes <- msigdb_data[grepl("palmitoyl|lipid", msigdb_data[["gs_name"]], ignore.case = TRUE), ]

# Extract gene names
gene_list <- unique(palmitoyl_genes$gene_symbol)

# Check the number of genes
length(gene_list)

write.csv(gene_list, file="palmitoylation_related_gene.csv")


############################################################################
library(data.table)
library(dplyr)
library(tidyverse)
sub.dir <- "HNSC"
setwd(file.path(main.dir, sub.dir))

HNSC_count=fread("TCGA-HNSC.star_counts.tsv", header=T, sep='\t',data.table=F)
HNSC_pro<-fread("gencode.v36.annotation.gtf.gene.probeMap",header=T,sep='\t',data.table=F)

HNSC_pro<-HNSC_pro[,c(1,2)]

HNSC_count_pro<-merge(HNSC_pro,HNSC_count,by.y="Ensembl_ID",by.x="id")


#

HNSC_count_pro<-distinct(HNSC_count_pro,gene,.keep_all=T)

rownames(HNSC_count_pro)<-HNSC_count_pro$gene
HNSC_count_pro<-HNSC_count_pro[,-c(1,2)]

dim(HNSC_count_pro)# 59427   566

HNSC_phe<-fread("TCGA-HNSC.clinical.tsv",header=T,sep='\t',data.table=F)
HNSC_phe$sample[1:5]

#
rownames(HNSC_phe)<-HNSC_phe$sample
HNSC_cli<-HNSC_phe[,c(15,12,45, 61,62)]
colnames(HNSC_cli)<-c("age","gender","stage","N-stage","T-stage")

#If you want to remove rows based on a specific column:
HNSC_cli <- HNSC_cli[HNSC_cli$stage != "" & !is.na(HNSC_cli$stage), ]

HNSC_cli$stage <- sub("Stage", "", HNSC_cli$stage)

HNSC_cli$`N-stage` <-sub("N", "", HNSC_cli$`N-stage`)

HNSC_cli$`T-stage`<-sub("T", "", HNSC_cli$`T-stage`)

HNSC_cli <- HNSC_cli[HNSC_cli$`N-stage` != "" & !is.na(HNSC_cli$`N-stage`), ]

save(HNSC_cli, file="HNSC_cli.rdata")

library(stringr)
tumor <- colnames(HNSC_count_pro)[as.integer(substr(colnames(HNSC_count_pro),14,15)) == 01]
normal <- colnames(HNSC_count_pro)[as.integer(substr(colnames(HNSC_count_pro),14,15)) == 11]
recur <- colnames(HNSC_count_pro)[as.integer(substr(colnames(HNSC_count_pro),14,15)) == 02]
metastatic <- colnames(HNSC_count_pro)[as.integer(substr(colnames(HNSC_count_pro),14,15)) == 06]

tumor_sample <- HNSC_count_pro[,tumor]
normal_sample <- HNSC_count_pro[,normal]
HNSC_expr_by_group <- cbind(tumor_sample,normal_sample)

saveRDS(HNSC_expr_by_group,file='HNSC_expr_by_group.rds')


###############################################################################################
# BiocManager::install("edgeR")
library(edgeR)


group <- c(rep("Tumor",520), rep("Normal",44))
group
class(group)

exprSet_by_group <- readRDS("HNSC_expr_by_group.rds")
group_list <- group

########################################################################################

########################################################################################
library(limma)
## Run limma
data <- exprSet_by_group
group_list <- factor(group_list)
design <- model.matrix(~0+group_list)
rownames(design) <- colnames(data)
colnames(design) <- levels(group_list)

DGElist <- DGEList(counts = data, group = group_list)
dim(DGElist)
#keep_gene <- rowSums(cpm(DGElist) > 1) >= 2
#table(keep_gene)
#DGElist <- DGElist[keep_gene, , keep.lib.sizes = F]
DGElist <- calcNormFactors(DGElist)
dim(DGElist)
v <- voom(DGElist, design, plot = T, normalize = "quantile")
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(contrasts = c("Tumor-Normal"), levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

nrDEG_limma_voom <- topTable(fit2, coef = "Tumor-Normal", n = Inf)

nrDEG_limma_voom <- na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
dim(nrDEG_limma_voom) #59427     6

PValue <- 0.05
foldChange <- log2(1.5)
nrDEG_limma_voom_signif <- nrDEG_limma_voom[(nrDEG_limma_voom$P.Value < PValue & 
                                               (nrDEG_limma_voom$logFC > foldChange | nrDEG_limma_voom$logFC < (-foldChange))),]
dim(nrDEG_limma_voom_signif) #4138    6
nrDEG_limma_voom_signif1 <- nrDEG_limma_voom_signif[row.names(nrDEG_limma_voom_signif)!="HTR3A",]
dim(nrDEG_limma_voom_signif1) #4137    6


nrDEG_limma_voom_signif <- nrDEG_limma_voom_signif[order(nrDEG_limma_voom_signif$logFC),]
#save(nrDEG_limma_voom_signif, file = "nrDEG_limma_voom_signif.Rdata")
write.csv(nrDEG_limma_voom_signif, file = "nrDEG_limma_voom_signif.csv")


########################################################################################
#create volcano plot
DEGs1<-nrDEG_limma_voom
dir.create(file.path(main.dir, "output"), showWarnings = FALSE)


DEGs1$change <-ifelse(nrDEG_limma_voom$adj.P.Val < PValue & abs(nrDEG_limma_voom$logFC) >= foldChange, 
                      ifelse(nrDEG_limma_voom$logFC> foldChange ,'Up','Down'),
                      'Stable')
g <- ggplot(data=DEGs1,
       aes(x=logFC,
           y=-log10(adj.P.Val)))+
  geom_point(alpha=0.4,size=3.5,
             aes(color=change))+
  ylab("-log10(padj)")+
  scale_color_manual(values=c("blue4","grey","red3"))+
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=4,col="black",lwd=0.8)+
  theme_bw()

ggsave("../output/vocano.pdf", g, height=6, width=8)

##########################################################################################
#plot heatmap

library(pheatmap)

DEG=read.csv("nrDEG_limma_voom_signif.csv")
# Select top 100 DEGs by absolute logFC

top100_genes <- DEG %>% arrange(desc(abs(logFC))) %>% head(100) %>% pull(X)

DEG_matrix=exprSet_by_group[top100_genes,]

# Create metadata based on column names

Condition <- ifelse(as.integer(substr(colnames(exprSet_by_group), 14, 15)) == 01, "Tumor",
                    ifelse(as.integer(substr(colnames(exprSet_by_group), 14, 15)) == 11, "Normal", NA))

metadata <- data.frame(
  SampleID = colnames(exprSet_by_group),  # Sample names must match column names of count matrix
  Condition = Condition)  # Assign condition


# Set row names to match sample IDs
rownames(metadata) <- metadata$SampleID  
metadata=metadata[,-1,drop=FALSE]


# Create a sample annotation dataframe
annotation_col <- metadata

# Define colors for annotation
ann_colors <- list(SampleType = c("Normal" = "blue", "Tumor" = "red"))

g <- pheatmap(DEG_matrix,
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_col = annotation_col,  # Add annotation
         annotation_colors = ann_colors,  # Add color mapping
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         fontsize_row = 6,
         main = "Heatmap of Differentially Expressed Genes")

ggsave("../output/heatmap.pdf", g, height=6, width=8)

##############################282DEGS###############################################

library(VennDiagram)

# Read a CSV file

HNSC_limma<- read.csv("nrDEG_limma_voom_signif.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
HNSC_limma<-HNSC_limma$X

PAM_Genes<-read.csv("palmitoylation_related_gene.csv", header = TRUE,sep = ",", stringsAsFactors = FALSE)
PAM_Genes<-PAM_Genes$x

shareGenes <- intersect(HNSC_limma,PAM_Genes)

shareGenes<- data.frame(Gene=shareGenes)

print(shareGenes)

write.table(shareGenes,file="HNSC_limma_PAM_shared_gene.txt",quote = F,sep='\t',row.names = T)
venn.plot <- venn.diagram(
  x = list(
    'DEGs' = HNSC_limma,
    'Palmitoylation Genes' = PAM_Genes
  ),resolution = 300,
  filename = NULL,
  col = "black",
  fill = c("dodgerblue", "seagreen3"),
  alpha = 0.6,
  scaled = FALSE,
  cex = 1.4,
  cat.col = "black",
  cat.cex = 1.4,
  cat.frontface = "bold",
  cat.pos = c(-15, 20),  
  cat.dist = 0.03,
  margin = 0.05,
  main = "",
  main.cex = 0.5
)

pdf(file = "../output/HNSC_Limma_PAM_Venn02.pdf")
grid.draw(venn.plot)
dev.off()


######################### volcano& heatmap of 282 DEGs######################################

#create Volcano plot

W_DEG<-read.csv("nrDEG_limma_voom_signif.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)

rownames(W_DEG)<-W_DEG$X
W_DEG<-W_DEG[,-1]


s_DEG<-read.table("HNSC_limma_PAM_shared_gene.txt",header=TRUE,sep='\t')
s_DEG<-s_DEG$Gene

m_DEG<-W_DEG[s_DEG,]
write.csv(m_DEG, file="282DEG_limma_voom.csv")



DEGs1<-m_DEG
PValue<-0.05
foldChange<-log2(1.5)

DEGs1$change <-ifelse(DEGs1$adj.P.Val < PValue & abs(DEGs1$logFC) >= foldChange, 
                      ifelse(DEGs1$logFC> foldChange ,'Up','Down'),
                      'Stable')

graphics.off()
g <- ggplot(data=DEGs1,
       aes(x=logFC,
           y=-log10(adj.P.Val)))+
  geom_point(alpha=0.4,size=3.5,
             aes(color=change))+
  ylab("-log10(padj)")+
  scale_color_manual(values=c("blue4","grey","red3"))+
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=4,col="black",lwd=0.8)+
  theme_bw()

ggsave("../output/vocalno_282_DEGs.pdf", g, height=6, width=8)


#########################################################################################

#Plot heatmap

exprSet_by_group <- readRDS("HNSC_expr_by_group.rds")

# Select top 100 DEGs by absolute logFC
DEG=read.csv("282DEG_limma_voom.csv")

top100_genes <- DEG %>% arrange(desc(abs(logFC))) %>% head(100) %>% pull(X)

DEG_matrix=exprSet_by_group[top100_genes,]

# Create metadata based on column names

Condition <- ifelse(as.integer(substr(colnames(exprSet_by_group), 14, 15)) == 01, "Tumor",
                    ifelse(as.integer(substr(colnames(exprSet_by_group), 14, 15)) == 11, "Normal", NA))

metadata <- data.frame(
  SampleID = colnames(exprSet_by_group),  # Sample names must match column names of count matrix
  Condition = Condition)  # Assign condition


# Set row names to match sample IDs
rownames(metadata) <- metadata$SampleID  
metadata=metadata[,-1,drop=FALSE]


# Create a sample annotation dataframe
annotation_col <- metadata

# Define colors for annotation
ann_colors <- list(SampleType = c("Normal" = "blue", "Tumor" = "red"))

g <- pheatmap(DEG_matrix,
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_col = annotation_col,  # Add annotation
         annotation_colors = ann_colors,  # Add color mapping
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         fontsize_row = 6,
         main = "Heatmap of Differentially Expressed Genes")

ggsave("../output/heatmap_282_DEGs.pdf", g, height=6, width=8)


#########################Go Circle######################################################
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(GOplot)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#### GOCircle from final DEGs ####
data <- read.csv("282DEG_limma_voom.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)

ego_BP <- clusterProfiler::enrichGO(gene = data$X,
                                    ont = "BP", OrgDb = org.Hs.eg.db,
                                    keyType = "SYMBOL")

ego_BP2 <- clusterProfiler::simplify(ego_BP, cutoff=0.7, by="p.adjust", select_fun = min)

ego <- data.frame(ego_BP2) 
ego <- ego[1:10, c("ID", "Description","geneID","p.adjust")]

ego$geneID <- stringr::str_replace_all(ego$geneID, "/", ",") 
names(ego) = c("ID","Term","Genes","adj_pval")
ego$Category = "BP"

genes <- data.frame(ID = data$X, logFC = data$logFC)

circ <- GOplot::circle_dat(ego, genes)

GOplot::GOCircle(circ)

pdf("../output/GOCircle_plot1.pdf",width=14,height=7)
GOplot::GOCircle(circ)
dev.off()

###########################KEGG analysis #########################################################

library(DOSE)
library(enrichplot)

s_DEG<-read.table("HNSC_limma_PAM_shared_gene.txt",header=TRUE,sep='\t')

gene.df <- bitr(gene = s_DEG$Gene, fromType = "SYMBOL", toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
g_list<-gene.df[,2]
head(g_list)

ek <- enrichKEGG(gene = g_list, #
                 organism = "hsa",  #
                 pvalueCutoff =0.05, #
                 qvalueCutoff = 0.05) #

#write.table(ek,file="ek.txt",sep="\t",quote=F,row.names = F)   

ek2 = setReadable(ek, #
                  OrgDb = "org.Hs.eg.db", 
                  keyType = "ENTREZID") #
head(ek2@result$geneID)
#write.table(ek2,file="ek2.txt",sep="\t",quote=F,row.names = F)  

pdf(file="../output/ek_barplot.pdf",width = 7,height = 5) 
barplot(ek, x = "GeneRatio", color = "p.adjust", #
        showCategory =10) 
dev.off()


#######Univariate Cox analysis###########################################################


HNSC_exp=readRDS("HNSC_expr_by_group.rds")
PAM=read.table("HNSC_limma_PAM_shared_gene.txt", header=T, sep="\t")

PAM=PAM$Gene
HNSC_PAM_exp=HNSC_exp[PAM,]


HNSC_surv=read.table("TCGA-HNSC.survival.tsv", header=T,sep="\t")


rownames(HNSC_surv)=HNSC_surv$sample
HNSC_surv=HNSC_surv[,-1]
HNSC_surv=HNSC_surv[,-3]

## 
Exp_surv=intersect(rownames(HNSC_surv),colnames(HNSC_PAM_exp))
HNSC_surv=HNSC_surv[Exp_surv,]
HNSC_PAM_exp=HNSC_PAM_exp[,Exp_surv]

HNSC_exp=HNSC_exp[,Exp_surv]

save(HNSC_surv, file="meta.rdata")
save(HNSC_exp, file="HNSC_exp.rdata")

exp=t(HNSC_PAM_exp)
COX_data=cbind(exp,HNSC_surv)

#

#
library(survival)
library(survminer)


rt <-COX_data

problematic_values <- rt$OS.time[is.na(as.numeric(rt$OS.time))]
print(problematic_values)
rt_clean <- rt[complete.cases(rt), ]
dim(rt_clean)
head(rt_clean)
rt=rt_clean
rt$OS.time=rt$OS.time/365

#
coxPfilter=0.05 
sigGenes <- c()  # Initialize sigGenes
outTab <- data.frame()  # Initialize output table

for (i in colnames(rt[,3:ncol(rt)])) {
  # Perform Cox regression
  cox <- coxph(Surv(OS.time, OS) ~ rt[, i], data = rt)
  coxSummary <- summary(cox)
  
  # Extract p-value
  coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
  
  # Check for NA before filtering
  if (!is.na(coxP) && coxP < coxPfilter) {
    sigGenes <- c(sigGenes, i)
    outTab <- rbind(outTab,
                    cbind(id = i,
                          HR = coxSummary$conf.int[, "exp(coef)"],
                          HR.95L = coxSummary$conf.int[, "lower .95"],
                          HR.95H = coxSummary$conf.int[, "upper .95"],
                          pvalue = coxP)
    )
  }
}


#
write.table(outTab,file="TCGA_uniCox.txt",sep="\t",row.names=F,quote=F)

#
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="TCGA_uniSigExp.txt",sep="\t",row.names=F,quote=F)


#################Lasso Regression ##############################################

# Lasso regression
proj = "TCGA-HNSC"

load("meta.rdata")
meta=HNSC_surv
expr=read.table("TCGA_uniSigExp.txt", header=T, sep="\t")
rownames(expr)=expr$id
expr=expr[,-1]
expr=expr[,-50]

x = expr  # 
x <- as.matrix(x)
mode(x) <- "numeric"  # Ensure numeric type
y = meta$OS

library(glmnet)
set.seed(1006)
cv_fit <- cv.glmnet(x=x, y=y)
plot(cv_fit)

lambda_min <- cv_fit$lambda.min

# 
coef_path <- coef(cv_fit, s = "lambda.1se")  # 
coef_path <- as.matrix(coef_path)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE)

# 
model_lasso_min <- glmnet(x=x, y=y, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, lambda=cv_fit$lambda.1se)

coef(model_lasso_min)

choose_gene_min = rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta) != 0]
choose_gene_1se = rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta) != 0]

# 
save(choose_gene_min, file = paste0(proj, "_lasso_choose_gene_min.Rdata"))
save(choose_gene_1se, file = paste0(proj, "_lasso_choose_gene_1se.Rdata"))

write.csv(choose_gene_1se,file="HNSC_lasso_choose_gene_1se.csv")

####################Multivariate Cox Analysis##################################
# Cox regression
library(survival)
library(survminer)

# 
cox_data <- expr[,choose_gene_1se]
cox_data$time <- meta$OS.time
cox_data$status <- meta$OS
save(cox_data, file="TCGA_cox_data.rdata")


cox_model <- coxph(Surv(time, status) ~ CYP3A4+UNC13C+BRINP1+TIMP4+HCAR1+RSPO1+CYP4A11+HTR2C, data = cox_data)
summary(cox_model)

cox_summary <- summary(cox_model)
cox_table <- data.frame(
  Gene = rownames(cox_summary$coefficients),
  Coef = cox_summary$coefficients[, "coef"],
  HR = exp(cox_summary$coefficients[, "coef"]),
  HR_95L = exp(cox_summary$conf.int[, "lower .95"]),
  HR_95H = exp(cox_summary$conf.int[, "upper .95"]),
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

# Round values for better readability
cox_table <- cox_table %>%
  mutate(
    Coef = round(Coef, 3),
    HR = round(HR, 3),
    HR_95L = round(HR_95L, 3),
    HR_95H = round(HR_95H, 3),
    p_value = format.pval(p_value, digits = 3, eps = 0.001)
  )

cox_table

library(flextable)

cox_table %>%
  flextable() %>%
  theme_vanilla() %>%
  autofit()

write.csv(cox_table, "multivariate_cox_results.csv", row.names = FALSE)


g <- ggforest(model = cox_model,  # Cox model
         data = cox_data,  # Explicitly provide dataset
         main = "",
         cpositions = c(0.05, 0.15, 0.35),  
         fontsize = 1,  
         refLabel = 'reference',  
         noDigits = 3)

ggsave("../output/cox.pdf", g, height=6, width=8)

# 
cox_pred <- predict(cox_model, type = "risk")
cox_data$risk_group <- ifelse(cox_pred > median(cox_pred), "High", "Low")

save(cox_data,file="Cox_risk_group.rdata")

# Add risk scores to the dataset
cox_data$risk_score <- cox_pred
save(cox_data, file="Cox_riskscore.rdata")



cox_data <- cox_data[order(cox_data$risk_score), ]  # Sort by risk score
cox_data$PatientIndex <- 1:nrow(cox_data)  # Assign patient index

g <- ggplot(cox_data, aes(x = PatientIndex, y = risk_score, color = risk_group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "", x = "Patients (Sorted)", y = "Risk Score") +
  theme_minimal()
ggsave("../output/cox_risk.pdf", g, height=6, width=8)

g <- ggplot(cox_data, aes(x = PatientIndex, y = time, color = as.factor(status))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "red"), labels = c("Alive", "Dead")) +
  labs(title = "",
       x = "Patients (Increasing Risk Score)",
       y = "Survival Time (Days)",
       color = "Patient Status") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../output/cox_time.pdf", g, height=6, width=8)

#Extract the signature genes for heatmap 
cox_data_new<-cox_data[,c("CYP3A4","UNC13C","BRINP1","TIMP4","HCAR1","RSPO1","CYP4A11","HTR2C")]

# Transpose the dataset so rows = genes and columns = patients

gene_expression_matrix <- as.matrix(t(cox_data_new))  

# Assuming `RiskGroup` column is in `cox_data`
annotation_col <- data.frame(RiskGroup = cox_data$risk_group)
rownames(annotation_col) <- rownames(cox_data)  # Set patient IDs as row names

library(pheatmap)

# Define color scheme for the risk groups
ann_colors <- list(risk_group = c("High Risk" = "red", "Low Risk" = "blue"))

# Create the heatmap
g <- pheatmap(gene_expression_matrix, 
         scale = "row",  # Normalize by gene
         cluster_rows = TRUE,  # Cluster genes
         cluster_cols = FALSE,  # Keep patients ordered
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         main = "")

ggsave("../output/cox_heat.pdf", g, height=6, width=8)


fit <- survfit(Surv(time, status) ~ risk_group, data = cox_data)
g <- ggsurvplot(fit, data = cox_data, pval = TRUE, conf.int = TRUE,risk.table = TRUE,risk.table.height=0.25,
           title = "")
pdf("../output/cox_survival.pdf")
print(g, newpage = FALSE)
dev.off()



#####################Expression and Survival###################################
library(pROC)
library(survival)
library(survminer)

load("TCGA_cox_data.rdata")

# Create a binary prognosis outcome (e.g., short vs long survival)
cox_data$Prognosis <- ifelse(cox_data$time <= median(cox_data$time), "Poor", "Good")

cox_data$Prognosis <- factor(cox_data$Prognosis, levels = c("Good", "Poor"))

# ROC analysis using the pROC package
roc_curve <- roc(response = cox_data$Prognosis, predictor = cox_data$TIMP4, levels = c("Good", "Poor"))

# Print AUC
print(paste("Area Under the Curve (AUC):", auc(roc_curve)))

# Get optimal threshold using Youden's index
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
print(optimal_coords)

# Add a column for high/low gene expression based on the threshold
cox_data$GeneGroup <- ifelse(cox_data$TIMP4 > optimal_coords$threshold, "High", "Low")

# Fit Kaplan-Meier survival curves
surv_obj <- Surv(time = cox_data$time, event = cox_data$status)
km_fit <- survfit(surv_obj ~ GeneGroup, data = cox_data)

# Perform log-rank test to get the p-value
log_rank_test <- survdiff(surv_obj ~ GeneGroup, data = cox_data)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

# Plot Kaplan-Meier survival curve with p-value
g <- ggsurvplot(
  km_fit,
  data = cox_data,
  xlab = "Time (days)",
  ylab = "Survival Probability",
  title = "",
  legend.title = "TIMP4  Gene Expression",
  legend.labs = c("Low", "High"),
  pval = TRUE,              # Automatically displays p-value
  pval.method = TRUE,       # Display method (e.g., Log-rank test)
  risk.table = TRUE,        # Add risk table
  risk.table.height = 0.2   # Adjust risk table height
)

pdf("../output/Kaplan_Meier.pdf")
print(g, newpage = FALSE)
dev.off()



######################Time-dependent ROC analysis##############################
library(timeROC)
library(survival)
library(ggplot2)

load("Cox_riskscore.rdata")
# Time-dependent ROC analysis

roc <- timeROC(T = cox_data$time,            # Survival time
               delta = cox_data$status,      # Event status (1 = event, 0 = censored)
               marker = cox_data$risk_score, # Predicted risk score
               cause = 1,                    # Event of interest
               times = c(1, 3, 5) * 365)     # Time points in days (1 year, 3 years, 5 years)

# Plot ROC curves
pdf("../output/Roc.pdf")
plot(roc, time = 365, col = "red", lwd = 2, title=FALSE)      # 1 year ROC curve
plot(roc, time = 1095, col = "blue", add = TRUE, title = FALSE)  # 3 year ROC curve
plot(roc, time = 1825, col = "green", add = TRUE, title=FALSE) # 5 year ROC curve

# Add legend
legend("bottomright", 
       legend = c(paste0("1-year AUC: ", round(roc$AUC[1], 3)),
                  paste0("3-year AUC: ", round(roc$AUC[2], 3)),
                  paste0("5-year AUC: ", round(roc$AUC[3], 3))),
       col = c("red", "blue", "green"), lwd = 2)

dev.off()



#########################GSE41613 data prepossessing##################################################
# BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
library(dplyr)
# BiocManager::install(c("WGCNA", "devtools", "impute", "preprocessCore", "GO.db"))
# library(devtools)
# install_github('jdrudolph/PerseusR')
library(WGCNA)
library(stringr)
library(flashClust)
library(iterators)
# data acquisition and preprocessing
# download GSE16561 dataset
gset <- getGEO('GSE41613',
               AnnotGPL = TRUE,   # download annotGPL
               getGPL = TRUE)     # download getGPL 
gset[[1]]
exp<-exprs(gset[[1]]) #  to get expression data
cli<-pData(gset[[1]]) # to get clinical data
GPL<-fData(gset[[1]]) # to get platform data
gpl<-GPL[,c(3,1)] #only need column 3 and column 1
colnames(gpl)<-c("Symbol","Probe_Id") # assign the column names
#gpl$Symbol<-data.frame(sapply(gpl,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1] #
exp<-as.data.frame(exp) # convert the object exp into a data frame
exp$Probe_Id<-rownames(exp) #add a column into the data frame in order for the merge function later.

exp_symbol<-merge(gpl,exp,by="Probe_Id")

# to remove if one probe corrsponds to multiple genes,
exp_symbol <- exp_symbol %>%
  distinct(Probe_Id, .keep_all = TRUE) #unique can do,but unique is the base R, distinct used in dplyr
# No one has been removed

# to remove empty symbol
exp_symbol <- filter(exp_symbol, Symbol!= "") 

#to get the average value if multiple probes correspond to one gene
expr<-exp_symbol[,-1]  # to remove the first column in order for the group_by function later

# to get the average number if the same Symbol has multiple expression data
exp_symbol_mean <- expr %>%
  group_by(Symbol) %>%
  summarise(across(everything(), mean, .names = "{col}")) %>%
  ungroup()

Normalize_EXP<-exp_symbol_mean

Normalize_EXP1<-Normalize_EXP[,-1] #to make the data frame perfect
rownames(Normalize_EXP1)<- Normalize_EXP$Symbol 

#save(Normalize_EXP1, file="GSE41613.rda")
#save(Normalize_EXP1, file="GSE41613.rdata")

write.csv(Normalize_EXP1,file = "GSE41613_Normalize_EXP.csv")# save the data for further usage

# Extract multiple columns by name
cli_data <- cli[, c(14, 15)]
cli_data$characteristics_ch1.5 <- sub("fu time:", "", cli_data$characteristics_ch1.5)

#cli_data$characteristics_ch1.16 <- sub("^os_event:", "", cli_data$characteristics_ch1.16)

# Change multiple column names
colnames(cli_data)[colnames(cli_data) %in% c("characteristics_ch1.4", "characteristics_ch1.5")] <- c("OS", "OS.time")

cli_data$OS <- ifelse(cli_data$OS == "vital: Alive", 0, 1)

str(cli_data$OS)
table(cli_data$OS)
class(cli_data$OS)

save(cli_data, file="GSE41613_cli_da.rdata")


###############################################################################
exp<-read.csv("GSE41613_Normalize_EXP.csv")
rownames(exp)<-exp$X
exp<-exp[,-1]

choose_gene_1se<-read.csv("HNSC_lasso_choose_gene_1se.csv")
choose_gene_1se<-choose_gene_1se$x

val_expr<-exp[choose_gene_1se,]
#Normalize the microarray data
###########################################
library(preprocessCore)

# Quantile normalization
normlized_val_expr <- normalize.quantiles(as.matrix(val_expr))
# Z-score normalization
normlized_val_expr <- scale(normlized_val_expr)

# Restore row and column names
rownames(normlized_val_expr) <- rownames(val_expr)
colnames(normlized_val_expr) <- colnames(val_expr)


load("GSE41613_cli_da.rdata")
normlized_val_expr<-t(normlized_val_expr)

val_dat<-cbind(normlized_val_expr,cli_data)

colnames(val_dat)[colnames(val_dat) == "OS.time"] <- "time"
colnames(val_dat)[colnames(val_dat) == "OS"] <- "status"

colnames(val_dat)[6] <- "UGT1A10"
colnames(val_dat)[15] <- "BGLAP"
save(val_dat, file="GSE41613_val.data.rdata")

# Cox regression
library(survival)
library(survminer)

load("TCGA_cox_data.rdata")

cox_model <- coxph(Surv(time, status) ~ CYP3A4+UNC13C+BRINP1+TIMP4+HCAR1+RSPO1+CYP4A11+HTR2C, data = cox_data)
summary(cox_model)

mean_train <- colMeans(cox_data, na.rm = TRUE)
val_dat$UGT1A10[is.na(val_dat$UGT1A10)] <- mean_train
val_dat$BGLAP[is.na(val_dat$BGLAP)] <- mean_train

# Predict risk scores on validation data
cox_pred <- predict(cox_model, newdata = val_dat, type = "risk")

# Add predicted risk to the validation data
val_dat$risk_score <- cox_pred
val_dat$risk_group <- ifelse(cox_pred > median(cox_pred), "High", "Low")

#################################################################################
str(val_dat$time)
val_dat$time <- as.numeric(val_dat$time)



# Create KM curve
fit <- survfit(Surv(time, status) ~ risk_group, data =val_dat)

# Plot KM curve
g <- ggsurvplot(fit, data = val_dat,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("High Risk", "Low Risk"),
           palette = c( "#f87669", "#2fa1dd"))  

pdf("../output/GSE41613_Kaplan_Meier.pdf")
print(g, newpage = FALSE)
dev.off()


########################################################################################
####################################################################################
val_dat$time<-val_dat$time*30

library(timeROC)
library(survival)
library(ggplot2)

pdf("../output/GSE41613_Roc.pdf")

# Fit timeROC model
roc <- timeROC(T = val_dat$time,            # Survival time
               delta = val_dat$status,      # Event status (1 = event, 0 = censored)
               marker = val_dat$risk_score, # Predicted risk score
               cause = 1,                   # Event of interest
               times = c(2, 3, 5) * 365)    # Time points in days (1 year, 3 years, 5 years)

# Plot the first ROC curve (1-year)
plot(roc, time = 730, col = "red", lwd = 2, title=FALSE)   

# Add the 3-year and 5-year ROC curves
plot(roc, time = 1095, col = "blue", add = TRUE, lwd = 2) 
plot(roc, time = 1825, col = "green", add = TRUE, lwd = 2)

# Add legend
legend("bottomright", 
       legend = c(paste0("1-year AUC: ", round(roc$AUC[1], 3)),
                  paste0("3-year AUC: ", round(roc$AUC[2], 3)),
                  paste0("5-year AUC: ", round(roc$AUC[3], 3))),
       col = c("red", "blue", "green"), lwd = 2)


dev.off()

############################################################################################

# Add risk scores to the dataset

val_dat$risk_score <- cox_pred

val_dat <- val_dat[order(val_dat$risk_score), ]  # Sort by risk score
val_dat$PatientIndex <- 1:nrow(val_dat)  # Assign patient index

g <- ggplot(val_dat, aes(x = PatientIndex, y = risk_score, color = risk_group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "", x = "Patients (Sorted)", y = "Risk Score") +
  theme_minimal()
ggsave("../output/GSE41613_risk.pdf", g, height=6, width=8)

ggplot(val_dat, aes(x = PatientIndex, y = time, color = as.factor(status))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "red"), labels = c("Alive", "Dead")) +
  labs(title = "",
       x = "Patients (Increasing Risk Score)",
       y = "Survival Time (Days)",
       color = "Patient Status") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../output/GSE41613_time.pdf", g, height=6, width=8)

#Extract the signature genes for heatmap 
val_dat_new<-val_dat[,c("CYP3A4","UNC13C","BRINP1","TIMP4","HCAR1","RSPO1","CYP4A11","HTR2C")]

# Transpose the dataset so rows = genes and columns = patients

gene_expression_matrix <- as.matrix(t(val_dat_new))  

# Assuming `RiskGroup` column is in `cox_data`
annotation_col <- data.frame(RiskGroup = val_dat$risk_group)
rownames(annotation_col) <- rownames(val_dat)  # Set patient IDs as row names

library(pheatmap)

# Define color scheme for the risk groups
ann_colors <- list(risk_group = c("High Risk" = "red", "Low Risk" = "blue"))

# Create the heatmap
g <- pheatmap(gene_expression_matrix, 
         scale = "row",  # Normalize by gene
         cluster_rows = TRUE,  # Cluster genes
         cluster_cols = FALSE,  # Keep patients ordered
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         main = "")
ggsave("../output/GSE41613_heatmap.pdf", g, height=6, width=8)



###########################Nomogram and ROC analysis###########################
library(rms)

load("Cox_riskscore.rdata")
load("HNSC_cli.rdata")
index<-intersect(rownames(cox_data), rownames(HNSC_cli))
cox_data<-cox_data[index,]
HNSC_cli<-HNSC_cli[index,]

cox_data$time<-cox_data$time/365


HNSC_cli$stage <- sub("[A,B,C]$", "", HNSC_cli$stage)
HNSC_cli$gender <- ifelse(HNSC_cli$gender == "female", "F", "M")

newdata<-cbind(cox_data,HNSC_cli)


dd <- datadist(newdata)  # Ensure `data` contains SYT4
options(datadist = "dd")

# Fit Cox model using cph()
cox_model <- cph(Surv(time, status) ~ age+gender+stage+risk_score, data =newdata, x = TRUE, y = TRUE, surv = TRUE)

# Create a nomogram for 1-, 3-, and 5-year survival
nom <- nomogram(cox_model, fun = list(function(x) Survival(cox_model)(1, x),
                                      function(x) Survival(cox_model)(3, x),
                                      function(x) Survival(cox_model)(5, x)),
                funlabel = c("1-year survival", "3-year survival", "5-year survival"))

pdf("../output/nomogram.pdf")

plot(nom)

dev.off()

###################
load("HNSC_cli.rdata")
load("Cox_riskscore.rdata")
newdata<-cox_data[,c(17,18,20)]

index<-intersect(rownames(newdata), rownames(HNSC_cli))
newdata<-newdata[index,]
HNSC_cli<-HNSC_cli[index,]
HNSC_cli_riskscore<-cbind(newdata, HNSC_cli)
HNSC_cli_riskscore$gender <- ifelse(HNSC_cli_riskscore$gender == "female", "0", "1")
HNSC_cli_riskscore$stage <- sub("[A,B,C]$", "", HNSC_cli_riskscore$stage)

HNSC_cli_riskscore$stage <- trimws(toupper(HNSC_cli_riskscore$stage))
HNSC_cli_riskscore$stage <- suppressWarnings(as.numeric(as.roman(HNSC_cli_riskscore$stage)))

HNSC_cli_riskscore$`N-stage` <- sub("[a,b,c]$", "", HNSC_cli_riskscore$`N-stage`)
HNSC_cli_riskscore$`T-stage` <- sub("[a,b,c]$", "", HNSC_cli_riskscore$`T-stage`)

# Replace all 'X' values with 0 in the whole dataset
HNSC_cli_riskscore[HNSC_cli_riskscore == "X"] <- NA
HNSC_cli_riskscore <- na.omit(HNSC_cli_riskscore)
HNSC_cli_rs<-HNSC_cli_riskscore[,-1]

write.table(HNSC_cli_rs, file="HNSC_cli_sc.txt", sep="\t", row.names = TRUE, col.names = TRUE,quote = FALSE)

#############################################################################

library(pROC)

inputFile = "HNSC_cli_sc.txt"

rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
# 
y = colnames(rt)[1]


bioCol = c("red", "blue", "green", "yellow","turquoise","brown")


pdf(file = "../output/ROC2.pdf", width = 5, height = 5)


roc1 = roc(rt[, y], as.vector(rt[, 2]))
plot(1 - roc1$specificities, roc1$sensitivities, type = "l", 
     col = bioCol[1], xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity", 
     lwd = 2, cex.axis = 1, cex.lab = 1, main = "")


aucText = c(paste0(colnames(rt)[2], ", AUC = ", sprintf("%0.3f", auc(roc1))))


for (i in 3:ncol(rt)) {
  roc1 = roc(rt[, y], as.vector(rt[, i]))
  lines(1 - roc1$specificities, roc1$sensitivities, col = bioCol[i-1], lwd = 2) 
  aucText = c(aucText, paste0(colnames(rt)[i], ", AUC = ", sprintf("%0.3f", auc(roc1)))) 
}


legend("bottomright", aucText[order(-as.numeric(gsub(".*, AUC = ", "", aucText)))],
       lwd = 2, bty = "n", col = bioCol[1:(ncol(rt)-1)], cex = 0.8)

dev.off()


########################Immune Infiltration Analysis######################################################

options(stringsAsFactors = FALSE)

library(limma) 
library(dplyr)
library(tidyverse)

dir.create('immune_results')
expFile = "TCGA_HNSC.txt"
dir.create('immune_results/estimate')

library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)

tcga_exp=read.delim(expFile,sep='\t',header = T,check.names = F)
#
filterCommonGenes(input.f=expFile , output.f="immune_results/estimate/estimate_input.gct", id="GeneSymbol")


## ###
estimateScore("immune_results/estimate/estimate_input.gct", "immune_results/estimate/estimate_score.gct", platform="affymetrix")
ESTIMATE_score = read.table("immune_results/estimate/estimate_score.gct",
                            skip = 2,
                            header = T,
                            row.names = 1,check.names = F)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[,2:ncol(ESTIMATE_score)]))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]
row.names(ESTIMATE_score)=colnames(tcga_exp)
ESTIMATE_score=ESTIMATE_score[,-1]
write.csv(ESTIMATE_score, "immune_results/estimate/ESTIMATE_score.CSV")


load("Cox_risk_group.rdata")
dat_group <- cox_data[, "risk_group", drop = FALSE]

ESTIMATE_score$risk_group=dat_group$risk_group

library(ggpubr)
library(ggsci)
head(ESTIMATE_score)
PurityScore=ggplot(ESTIMATE_score,aes(x = risk_group, y = TumorPurity, fill = risk_group)) +
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values = ggsci::pal_aaas()(10))+theme_bw()+
  geom_signif(
    comparisons = list(
      c("High","Low")), 
    map_signif_level = T, 
    test = "t.test", 
    vjust=0.1, 
    tip_length = 0.05)+
  theme_classic(base_size = 16)+theme(legend.position = "top")

ImmuneScore=ggplot(ESTIMATE_score,aes(x = risk_group, y = ImmuneScore, fill = risk_group)) +
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values = ggsci::pal_aaas()(10))+theme_bw()+
  geom_signif(
    comparisons = list(
      c("High","Low")), 
    map_signif_level = T, 
    test = "t.test", 
    vjust=0.1, 
    tip_length = 0.05)+
  theme_classic(base_size = 16)+theme(legend.position = "top")

StromalScore=ggplot(ESTIMATE_score,aes(x = risk_group, y = StromalScore, fill = risk_group)) +
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values = ggsci::pal_aaas()(10))+theme_bw()+
  geom_signif(
    comparisons = list(
      c("High","Low")), 
    map_signif_level = T, 
    test = "t.test", 
    vjust=0.1, 
    tip_length = 0.05)+
  theme_classic(base_size = 16)+theme(legend.position = "top")

ESTIMATEScore=ggplot(ESTIMATE_score,aes(x = risk_group, y = ESTIMATEScore, fill = risk_group)) +
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values = ggsci::pal_aaas()(10))+theme_bw()+
  geom_signif(
    comparisons = list(
      c("High","Low")), 
    map_signif_level = T, 
    test = "t.test", 
    vjust=0.1, 
    tip_length = 0.05)+
  theme_classic(base_size = 16)+theme(legend.position = "top")


ESTIMATE=ggarrange(ImmuneScore,StromalScore,ESTIMATEScore,PurityScore,ncol = 2,nrow = 2,common.legend = T)
ESTIMATE
ggsave(ESTIMATE,filename = "../output/ESTIMATE.pdf",wi=8,he=8)


##

checkpoint=read.delim('checkpoint.txt',sep='\t',header = F,check.names = F)
tcga=read.table(expFile, header=T, sep="\t", check.names=F)
checkpoint_exp=tcga[checkpoint$V1,]
#
checkpoint_exp=na.omit(checkpoint_exp)
checkpoint_exp=as.data.frame(t(checkpoint_exp))
#
load("Cox_risk_group.rdata")
dat_group <- cox_data[, "risk_group", drop = FALSE]

identical(rownames(checkpoint_exp),rownames(dat_group))
checkpoint_exp$risk_group=dat_group$risk_group

#
Muti_Boxplot<-function(dat,group,group_cols,leg,
                       test_method = 'wilcox.test',ylabs){
  library(ggpubr)
  library(reshape2)
  dat1=reshape2::melt(cbind.data.frame(dat,group))
  p=ggboxplot(dat1, x='variable', y='value', 
              fill = "group", color = "black",
              palette = group_cols, 
              ylab=ylabs,xlab='',
              add = "boxplot")+ 
    stat_compare_means(aes(group=group),method = test_method,
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")+
    #theme(axis.text.x = element_text(angle = 0,hjust = 1))+
    labs(fill=leg)
  return(p)
}

checkpoint=Muti_Boxplot(dat =checkpoint_exp,
                        group = checkpoint_exp$risk_group,
                        group_cols = ggsci::pal_aaas()(9),
                        test_method = 'wilcox.test',
                        leg = 'Checkpoint',ylab = 'Gene expression')+ coord_flip()
checkpoint
ggsave("../output/checkpoint.pdf", checkpoint,wi=6,he=10)


##########
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
source('../Cibersort.R')


load("HNSC_exp.rdata")
exprSet<-HNSC_exp

LM22.file <- "LM22.txt"
TCGA_HNSC.file <- "TCGA_HNSC.txt"
# 1. Cibersort
write.table(exprSet, file = TCGA_HNSC.file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


results <- CIBERSORT(LM22.file, TCGA_HNSC.file, perm = 1000, QN = F)  

load("Cox_risk_group.rdata")

risk_group <- ifelse(cox_data$risk_group == "High", 1, 0)
TME_data <- as.data.frame(results[, 1:22])  
TME_data$group <- risk_group  
TME_data$sample <- row.names(TME_data)  

# 
TME_New <- melt(TME_data, id.vars = c("group", "sample"), variable.name = "Cell_Type", value.name = "Fraction")

TME_New$group <- ifelse(TME_New$group == 0, "low risk", "high risk")


# 
g <- ggplot(TME_New, aes(x = sample, y = Fraction, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ group, scales = "free_x") +  
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Immune Infiltration by Risk Group",
       x = "Sample",
       y = "Fraction",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c", "#ff7f00", 
                               "#984ea3", "#a65628", "#f781bf", "#999999",
                               "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
                               "#1b9e77", "#d95f02", "#7570b3", "#e7298a", 
                               "#66a61e", "#e6ab02"))

ggsave("../output/TME.pdf", g, wi=10,he=6)




###############################################################################
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
results <- read.table("CIBERSORT-Results.txt", sep = '\t', header = TRUE)
rownames(results)=results$Mixture
results=results[,-1]
theme_set(theme_bw())
results_NEW<-results
results_NEW<-as.data.frame(results_NEW)
results_NEW$risk<-as.factor(cox_data$risk_group)

results_NEW<-melt(results_NEW,
                  id.vars=c("P.value", "Correlation","RMSE", "risk"),
                  variable.name="immune_cell",
                  value.name = "proportion")
results_NEW$risk<-factor(results_NEW$risk)
results_NEW$immune_cell<-factor(results_NEW$immune_cell)
str(results_NEW)

pdf("../output/immune_infiltration.pdf")
ggplot(results_NEW, aes(x = immune_cell, y = proportion, fill = risk)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.colour = "black") +
  stat_compare_means(aes(group = risk), label = "p.signif", method = "wilcox.test", 
                     label.x = 1.5, size = 5, color = "black") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "Immune Cell Type",
    y = "Proportion",
    fill = "Risk Group",
    title = ""
  )
dev.off()


#########################GSEA analysis#####################################################
load("HNSC_exp.rdata")
load("Cox_risk_group.rdata")
library(limma)

# Create the design matrix
design <- model.matrix(~ 0 + cox_data$risk_group)
colnames(design) <- c("HighRisk", "LowRisk")

fit <- lmFit(HNSC_exp, design)
contrast.matrix <- makeContrasts(HighRisk - LowRisk, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

deg_results <- topTable(fit, coef = 1, number = Inf, adjust.method = "BH")
geneList <- deg_results$logFC
names(geneList) <- rownames(deg_results)
geneList <- sort(geneList, decreasing = TRUE)

library(clusterProfiler)
library(msigdbr)

gene_set_Hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_results_hallmark <- GSEA(geneList, 
                              TERM2GENE = gene_set_Hallmark, 
                              pvalueCutoff = 0.05,
                              minGSSize = 10,
                              maxGSSize = 500)

result_hallmark <- gsea_results_hallmark@result 

pdf("../output/gsea_hallmark_myogenesis.pdf")
gseaplot2(gsea_results_hallmark, 
          geneSetID = "HALLMARK_MYOGENESIS",
          title = "",
          color = "green",
)

dev.off()

pdf("../output/gsea_hallmark_dotplot.pdf")
dotplot(gsea_results_hallmark)
dev.off()

pdf("../output/gsea_hallmark_head_enrichment.pdf")
gseaplot2(gsea_results_hallmark, 
          geneSetID =,rownames(gsea_results_hallmark@result)[head(order(gsea_results_hallmark@result$enrichmentScore))]
)
dev.off()

pdf("../output/gsea_hallmark_tail_enrichment.pdf")
gseaplot2(gsea_results_hallmark, 
          geneSetID =,rownames(gsea_results_hallmark@result)[tail(order(gsea_results_hallmark@result$enrichmentScore))]
)
dev.off()



#######################drug sensitivity####################################################
#oncoPredict

options(stringsAsFactors = F)
# BiocManager::install(c("sva", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TCGAbiolinks"))
# install.packages("oncoPredict")
library(oncoPredict)

get_oncoPredict_res <- function(data = data,
                                traData = c('GDSC2', 'GDSC1', 'CTRP2')[1],
                                minNumSamples = 10) {
  library(oncoPredict)
  
  data <- as.matrix(data)
  if (traData == 'GDSC2') {
    traDataExp <- readRDS("../DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
    traDataRes <- readRDS(file = "../DataFiles/Training Data/GDSC2_Res.rds")
  } else if (traData == 'GDSC1') {
    traDataExp <- readRDS("../DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
    traDataRes <- readRDS(file = "../DataFiles/Training Data/GDSC1_Res.rds")
  } else if (traData == 'CTRP2') {
    traDataExp <- readRDS("../DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
    traDataExp <- log2(traDataExp + 1)
    traDataRes <- readRDS(file = "../DataFiles/Training Data/CTRP2_Res.rds")
  } else {
    stop('Please check if traData is correct, only support GDSC2, GDSC1, CTRP2 databases')
  }
  
  calcPhenotype(trainingExprData = traDataExp,
                trainingPtype = traDataRes,
                testExprData = data,
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = minNumSamples, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData')
}

# 
load("HNSC_exp.rdata")
dat_exp <- HNSC_exp

boxplot(dat_exp[, 1:5])
dat_exp[1:5, 1:5]

load("Cox_risk_group.rdata")

dat_group <- cox_data[, "risk_group", drop = FALSE]

head(dat_group)
table(dat_group$risk_group)

#  ###############
library(stringr)
if (file.exists('calcPhenotype_Output/DrugPredictions.csv')) {

  GDSC2_res <- read.csv('calcPhenotype_Output/DrugPredictions.csv',
                        row.names = 1, check.names = F, stringsAsFactors = F)
} else {
  GDSC2_res <- get_oncoPredict_res(data = dat_exp,
                                   traData = 'GDSC2')
}

GDSC2_res[1:5, 1:5]

# 
dat_GDSC2_res <- cbind(dat_group,
                       GDSC2_res)
dat_GDSC2_res[1:5, 1:6]
class(dat_GDSC2_res)
str(dat_GDSC2_res)

#  
library(ggsci)
library(ggplot2)
library(ggpubr)

#
Cisplatin <- ggplot(dat_GDSC2_res, aes(x=risk_group, y=Cisplatin_1005)) +
  geom_boxplot(aes(fill = risk_group), position=position_dodge(0.9), outlier.colour = NA, notch = T) +
  stat_compare_means(aes(group=risk_group)) +
  scale_fill_d3()+
  theme_bw()+ 
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5,  hjust = 0.5,colour="black"),
        legend.position = 'top') + xlab('risk_group') + 
  ylab('Estimated IC50 (Cisplatin)')#
Cisplatin
ggsave(Cisplatin,filename = "../output/Cisplatin.pdf",he=6,wi=6)


Fulvestrant <- ggplot(dat_GDSC2_res, aes(x=risk_group, y=Fulvestrant_1200)) +
  geom_boxplot(aes(fill = risk_group), position=position_dodge(0.9), outlier.colour = NA, notch = T) +
  stat_compare_means(aes(group=risk_group)) +
  scale_fill_d3()+
  theme_bw()+ 
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5,  hjust = 0.5,colour="black"),
        legend.position = 'top') + xlab('risk_group') + 
  ylab('Estimated IC50 (Fulvestrant)')#

Fulvestrant

ggsave(Fulvestrant,filename = "../output/Fulvestrant.pdf",he=6,wi=6)



