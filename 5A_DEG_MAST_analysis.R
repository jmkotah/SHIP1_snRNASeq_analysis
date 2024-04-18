# Author: Janssen Kotah
# snRNAseq analysis for WT/SHIP1 KO mice as part of Matera et al. project
# MAST analysis code adapted Astrid Alsema, this was analyzed on a different computing cluster, therefore not in jupyter notebook
# current code was for non-neuronal celltypes, similarly analyzed for neuronal subtypes

rm(list = ls())

library(lme4)
library(MAST)
library(Seurat)
library(reshape2)
library(data.table)
library(SingleCellExperiment)
library(textshape)
library(ggplot2)
library(dplyr)

Input.dir <- paste0("/scratch/p310674/20240305_SHIP1_Reseq_nonNeurons_MAST//")
batch_number = "01_"

# import the processed seurat object 
datasets <- readRDS(paste0(Input.dir, "/003B5_nonNeurons_Rd5_harmony_integration.rds"))

print("number of cells in this dataset: ")
print(nrow(datasets@meta.data))

#merge ependymal due to lack of secreted cells in WT
datasets@meta.data$non_neur_annot = case_when(datasets@meta.data$non_neur_annot %in% c("Ependymal_Ciliated", "Ependymal_Secretory") ~ "Ependymal", T ~ datasets@meta.data$non_neur_annot)

Idents(datasets) <- "non_neur_annot"

subtypes <- unique(datasets$non_neur_annot) %>% sort()
print(paste0("Number of celltypes to analyze: ", length(subtypes)))

for (each in 1:length(subtypes)){
print(paste0("Analyzing celltype: ", subtypes[each]))
print(Sys.time())

# subset out celltype of interest
data_new <- subset(datasets, idents = subtypes[each])

# remove cells with less than 5 cells from one sample in this subcluster
table(data_new$Genotype, data_new$sample) %>% t() %>% data.frame() %>%
    filter(Freq != 0) %>%
    filter(Freq > 5) %>%
    select(Var1) %>% unlist() -> keep.samples

data_new <- subset(data_new, sample %in% keep.samples)

##First part is done with seurat: identify features that are expressed by a min percentage of cells(min.pct) and with a min FC (logfc.threshold) between the comparison groups
temp <- data_new
DefaultAssay(temp) <- 'RNA'   

# feature selection thresholds
min.pct=0.05 #only test genes that are detected in at least 5% of cells in either of the two groups
min.diff.pct = -Inf  #  minimum difference in expression between the two groups, set to -Inf for our purpose
logfc.threshold=0.15 # select genes to be modeled, for computational efficiency select based on at least logFC difference of +/-  0.15 between the two groups of cells

# set reference group
temp$Genotype <- factor(temp$Genotype, levels = c("WT", "KO"))

# Calculate FC between the two groups  
fc.results<-Seurat::FoldChange(temp,  group.by="Genotype", ident.1="KO", ident.2 = "WT", # 
                  slot = 'data', pseudocount.use = 0.5)

print('dimensions before filter:')
print(dim(temp))

alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
names(x = alpha.min) <- rownames(x = fc.results)
features <- names(x = which(x = alpha.min > min.pct))
alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
features <- names(x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct))

# feature selection based on logFC
total.diff <- fc.results[, 1] 
names(total.diff) <- rownames(fc.results)
features.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
features <- intersect(x = features, y = features.diff)

# subset for the selected features 
seu_subset<-subset(x = temp, features = features)
DefaultAssay(seu_subset) <- "RNA"
print('dimensions after prefilter:')
print(dim(seu_subset))

logcounts <- as.matrix(GetAssayData(seu_subset, assay = "RNA", slot = 'data'))
meta <- as.data.frame(seu_subset@meta.data)
fdata <- data.frame(primerid=rownames(seu_subset), stringsAsFactors = F)

print('making mast object..')
sca <- MAST::FromMatrix(exprsArray = list(logNorm = logcounts), cData = meta, fData = fdata)
rm(logcounts, meta, fdata)

print('saving mast object')

saveRDS(sca, file = paste0(subtypes[each], '_-sca-filt.rds'))

# 1. prepare the donor-related variables inside the sca
sca$group <- factor(sca$Genotype, levels = c("WT", "KO"))
sca$sample <- factor(sca$sample)

#adding covariate to account for sequencing depth per cell
cdr <-colSums(assay(sca)>0)
sca$cngeneson <- scale(cdr)

## 2. fit the model. orig.ident is the donor ID 
zlmCond <-MAST::zlm(~ group + cngeneson + (1|sample),
                    sca, 
                    exprs_value='logNorm', 
                    method = 'glmer', 
                    ebayes = FALSE)

saveRDS(zlmCond, paste0(batch_number, "_majorCelltypes_", subtypes[each], '_zlmCond.rds'))

## 3. create contrast
summaryCond <- summary(zlmCond, doLRT="groupKO") 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='groupKO' & component=='H',.(primerid, `Pr(>Chisq)`)],
#hurdle P values
summaryDt[contrast=='groupKO' & component=='logFC', .(primerid, coef, ci.hi,ci.lo)], by='primerid') 

saveRDS(summaryDt, paste0(batch_number, subtypes[each], '_SummaryDT.rds'))

## 4. investigate and merge output
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='groupKO' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='groupKO' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficient, log10
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
setorder(fcHurdle, fdr)

# optional: print some sanity checks, some gene models might not converge
print(paste0('we modelled ', nrow(sca), 'genes for ', subtypes[each]))
print('number of converged models for continuous component:')
print(sum(zlmCond@converged[,1] == "TRUE")) # continuous component 

print('number of converged models for discrete component:')
print(sum(zlmCond@converged[,2] == "TRUE")) # discrete model 

# Write CSV out in a fast way
fwrite(fcHurdle, paste0(subtypes[each], "-MAST-mixed.csv")) # this is all the genes modelled
fwrite(fcHurdleSig, paste0(subtypes[each], "-MAST-mixed-sign.csv")) # this is only the significant genes

a<-c(dim(fcHurdleSig))
print(paste('number of DE genes: ', a[1]))

print(Sys.time())
print(warnings())

}

