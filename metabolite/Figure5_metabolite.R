setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)#‘3.5.1’
library(ggrepel)
library(ggpubr)
library(plyr)
library(tidyverse)
library(stringr)
library(patchwork)#‘1.3.0’
library(ggsci)
library("DEswan")
library("UpSetR")
library(VennDiagram)
library(openxlsx)
library(clusterProfiler)
library(RColorBrewer)
library(colorspace)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure5/'

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig5 Data prepare =====
load(file = paste0(dataDir,'/Metabolite_anno_ForKEGG.Rdata'),verbose = T)
class_color = unique(class$class_color)
names(class_color) = levels(class$Class)

load(file = paste0(dataDir,'/Metabolite_expFor_DEM.Rdata'),verbose = T)
metab.normal = as.data.frame(t(life.metab.normal))
metab.normal.log = log10(metab.normal)

metab_anno = data.frame(
  Metabolites_names = rownames(metab.normal.log),
  Metabolite_IDs = paste0('Metab_', seq(nrow(metab.normal.log)))
)
class$metabolite_id = metab_anno$Metabolite_IDs[match(class$metabolite, metab_anno$Metabolites_names)]

data = as.data.frame(t(metab.normal.log))
colnames(data) = metab_anno$Metabolite_IDs[match(colnames(data), metab_anno$Metabolites_names)]
data$group = str_sub(rownames(data), 1, 1)

phenotype = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Sample_phenotype.csv',header=TRUE,row.names = 1)
phenotype$sample = rownames(phenotype)
phenotype$Group = str_sub(rownames(phenotype),1,1)

data$sex <- phenotype$sex[match(rownames(data),rownames(phenotype))]
data$Age <- phenotype$Age[match(rownames(data),rownames(phenotype))]

data[1:4,1:4] 

## Fig5a b metabolomics loess =====
set.seed(2023)
data_scaled = data

non_metab_cols = c('Age','sex','group')
data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)] <- scale(data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)])
data_scaled$Age <- as.numeric(as.character(data_scaled$Age))
data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)] <- data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)][ , colSums(is.na(data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)])) == 0]
data_scaled[1:4,1:4]

# Reassign metabolites
metabs <- names(data_scaled)
metabs <- metabs[ metabs %!in% non_metab_cols ]

age_min = min(data_scaled$Age)
age_max = max(data_scaled$Age)
lo_predict_list <- lapply(metabs, run_loess, data = data_scaled)

# Merge lo_predict
lo_predict <- as.data.frame(lo_predict_list)
lo_predict$age <- seq(age_min, age_max, 1)
lo_predict <- dplyr::relocate(lo_predict, age) ##Change column order

#### Convert wide to long add Age Group ####
lo_predict_long <- gather(lo_predict, metab, expression, 2:ncol(lo_predict)) ## reshape data
lo_predict_long_group = mutate(lo_predict_long, group = ifelse(lo_predict_long$age == 0, 'A', 
                                                               ifelse(0< lo_predict_long$age & lo_predict_long$age < 7, 'B', #1-6
                                                                      ifelse(6< lo_predict_long$age & lo_predict_long$age < 13, 'C',#7-12
                                                                             ifelse(12< lo_predict_long$age & lo_predict_long$age < 19, 'D',#13-18
                                                                                    ifelse(18< lo_predict_long$age & lo_predict_long$age < 41, 'E',#19-40
                                                                                           ifelse(40< lo_predict_long$age & lo_predict_long$age < 66, 'F', 'G' #41-65,65-84
                                                                                           )
                                                                                    )
                                                                             )
                                                                      )
                                                               )
)
)
table(lo_predict_long_group$age, lo_predict_long_group$group)

# Order metabolites by hclustering
clust_order <- hclust(dist(t(lo_predict[,2:ncol(lo_predict)])))$order
lo_predict_long_group$metab <- factor(lo_predict_long_group$metab, levels = unique(lo_predict_long_group$metab)[clust_order])
levels(lo_predict_long_group$metab)
# Add limits for heatmap plotting
zscore_limit <- 1
lo_predict_long_plotting <- lo_predict_long_group
lo_predict_long_plotting[ which(lo_predict_long_plotting$expression < -zscore_limit), "expression"] <- -zscore_limit
lo_predict_long_plotting[ which(lo_predict_long_plotting$expression > zscore_limit), "expression"] <- zscore_limit

# Generate heatmap
p <- ggplot(lo_predict_long_plotting , aes(x = age, y = metab, fill = expression)) +
  theme_bw() +
  geom_raster(interpolate=TRUE) +
  theme(axis.text.y=element_blank()) +
  scale_fill_gradient2(low = "#384f91", mid = "white", high = "#9E4434",
                       breaks = c(-zscore_limit, zscore_limit),
                       limits = c(-zscore_limit, zscore_limit))

p

ggsave(file = paste0(outDir,"metabolite_loess_predicted_heatmap.pdf"),p, width = 5, height = 4)

cols <- c(colnames(data_scaled[ ,(names(data_scaled) %!in% non_metab_cols)]), "Age")
data_scaled_subset <- data_scaled[ ,cols ]

data_long <- gather(data_scaled_subset, metab, expression, -Age)

line_plot_of_loess(data_long, alpha = 0.8, color = '#2D2D2D')

ggsave(file = paste0(outDir,"metabolite_loess_trajectories_norm.pdf"), width = 5, height = 4)

data_long$metabolite = metab_anno$Metabolites_names[match(data_long$metab, metab_anno$Metabolite_IDs)]
openxlsx::write.xlsx(data_long,file = paste0(outDir, '/Metabolite_scaledForLoess.xlsx'),rowNames = TRUE)

lo_predict_long_plotting$metabolite = metab_anno$Metabolites_names[match(lo_predict_long_plotting$metab, metab_anno$Metabolite_IDs)]
openxlsx::write.xlsx(lo_predict_long_plotting,file = paste0(outDir, '/Metabolite_Loess_predict_z-score.xlsx'),rowNames = TRUE)

## Fig5c metabolomics loess clusters =====

dist_mat <- dist(t(lo_predict[,2:ncol(lo_predict)]))

# Plot hierarchal clusters
pdf(paste0(outDir, "/Metabolite_dendogram_loess_predicted.pdf"))
clust <- hclust(dist_mat, method = "complete")
plot(clust, labels = FALSE)
dev.off()

save(dist_mat,lo_predict,clust,file = paste0(outDir, '/Metabolite_loess_cluster.Rdata'))
load(file = paste0(outDir, '/Metabolite_loess_cluster.Rdata'),verbose = T)

k = 8
hclust_cut <- data.frame(cutree(clust, k = k))
color_vector <- colorRampPalette(pal_npg("nrc", alpha = 0.9)(5))(k)
hclust_plot_list = list()
for (i in seq(k)) {
  cols <- c(rownames(hclust_cut[ which(hclust_cut[,1] == i), ,drop = FALSE]), "Age")
  data_scaled_subset <- data_scaled[ ,cols ]
  
  data_long <- gather(data_scaled_subset, metab, expression, -Age)
  
  plot <- line_plot_of_loess(data_long, color = color_vector[i])
  
  # Append plot to list
  hclust_plot_list[[i]] <- plot
}

# Generate figure of clustered trajectories
patchwork <- wrap_plots(hclust_plot_list, guides = "collect")
patchwork <- patchwork + plot_annotation(subtitle = "AllSample")

pdf(paste0(outDir, "/hclust_k", k, ".pdf"))
print(patchwork)
dev.off()


hclust_cut$metabolite_id = rownames(hclust_cut)
colnames(hclust_cut)[1] = 'cutree_hclust_k8'
hclust_cut_k8 = merge(hclust_cut, class, by = 'metabolite_id')
save(hclust_cut_k8, class, file = paste0(outDir, '/Hclust8_anno.Rdata'))
write.xlsx(hclust_cut_k8, file = paste0(outDir, '/Loess_hclustcomplete_cutk8.xlsx'))

data_keg = list()
for (i in 1:8) {
  data_keg[[i]] =  MarkerKEGG(hclust_cut_k8$metabolite[ which(hclust_cut_k8$cutree_hclust_k8 == i)], path_metab = path_metab, path_intro = path_intro,
                              minGSSize = 2, pAdjustMethod = "fdr",
                              pvalueCutoff = 0.99, outDir = outDir, Clust = paste0('clust_',i),
                              return_data = TRUE
  )
}


names(data_keg) = paste0("clust_", 1:8)

save(data_keg, file = paste0(outDir, '/Hclust_k8KEGG.Rdata'))

load(file = paste0(outDir, '/Hclust_k8KEGG.Rdata'))

path_num = 4
dataMerge = as.data.frame(matrix(nrow=length(data_keg)*path_num,ncol=length(data_keg)+1)) ## top10 pathway 2 crests
colnames(dataMerge) = c('Description',names(data_keg))

dataMerge$Description = as.vector(sapply(data_keg, function(x) x$Description[1:path_num]))
dataMerge = dataMerge[!duplicated(dataMerge$Description),]
rownames(dataMerge) = dataMerge$Description

for (i in names(data_keg)) {
  data_keg[[i]][,'-LogP'] = -log10(data_keg[[i]][,'p.adjust'])
  dataMerge[,i] = data_keg[[i]][,'-LogP'][match(rownames(dataMerge),data_keg[[i]]$Description)]
}

dataMerge[is.na(dataMerge)] <- 0

##### Dot plot ========
dataMerge.long <-gather(dataMerge,key=Change,value=LogPvalue,-Description)
for (i in names(data_keg)) {
  dataMerge.long[dataMerge.long$Change == i,'FoldEnrich'] = data_keg[[i]][,'FoldEnrich'][match(dataMerge.long$Description[dataMerge.long$Change == i],data_keg[[i]]$Description)]
}
dataMerge.long$Description <- factor(dataMerge.long$Description, levels = rev(rownames(dataMerge)))
dataMerge.long$Change <- factor(dataMerge.long$Change, levels = names(data_keg))
dataMerge.long = dataMerge.long[!is.na(dataMerge.long$FoldEnrich),]

cat('Max size: ', max(dataMerge.long$FoldEnrich),'; Min size: ', min(dataMerge.long$FoldEnrich), '\n',
    'Max color: ', max(dataMerge.long$LogPvalue),'; Min color: ', min(dataMerge.long$LogPvalue), '\n'
)

ggplot(dataMerge.long, aes(Change, Description)) +
  geom_point(aes(fill = LogPvalue, size = FoldEnrich), color = "black", shape = 21) +
  scale_size(range = c(1, 8), breaks = c(5,25,50,75)) +
  scale_fill_gradientn(colours = colorRampPalette(rev(sequential_hcl(5, palette = "Reds")))(50),breaks = c(3, 6, 9,12)) +
  theme_bw() +
  coord_fixed(ratio = 1.2) +
  labs(x = "", y = "", title = "", 
       fill = "-Log10 (p.adjust)", size = "Fold Enrichment") +
  theme(#aspect.ratio = 1.6,
    axis.text.x = element_text(size = 10,  angle = 90, hjust = 1, vjust = 0.5,
                               face = "plain", colour = "black"), 
    axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
  ) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black")
  )

ggsave(file = paste0(outDir,"/hclust8_KEGG_top4Dotplot.pdf"),width = 6,height = 7)

## Fig5h metabolomics DEswan remove Group A =====
data = data[data$group != 'A', ]
padj_deswan.thresh <- 0.05
beta_deswan.thresh <- 0.05 ##coefficients

res.DEswan=DEswan(data.df = data[,!(colnames(data) %in% c('sex','Age','group'))],
                  qt = data[,'Age'],
                  window.center = seq(2,80,5),
                  covariates = data[,colnames(data) %in% c('sex')],
                  buckets.size = 20
)

save(res.DEswan,file = paste0(outDir,'/DEswan_res_win5_bucket20.Rdata'))
load(file = paste0(outDir,'/DEswan_res_win5_bucket20.Rdata'))
head(res.DEswan$p)
head(res.DEswan$coeff)

# Convert to wide and calculate BH
res_wide_p <- reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")
res_wide_q <- q.DEswan(res_wide_p, method = "BH")

# Find coefficients
res_wide_b <- reshape.DEswan(res.DEswan, parameter = 2, factor = "qt")

# Subset BH and coefficients
q <- lapply(res_wide_q, function(x) {
  res_wide_q$variable[which(x <= padj_deswan.thresh)]
})
b <- lapply(res_wide_b, function(x) {
  res_wide_b$variable[which((x >= beta_deswan.thresh) | (x <= -beta_deswan.thresh))]
})
q[["variable"]] <- NULL
b[["variable"]] <- NULL

# Subset sig metabolites
sig_genes_list <- list()
for (age in names(b)) {
  sig_genes_list[[age]] <- intersect(q[[age]], b[[age]])
}
sig_genes <- sapply(sig_genes_list, length)

# Generate dataframe of number of sig metabolites
sig_gene_df <- data.frame(sig_genes)
names(sig_gene_df) <- 'Metabolites'
sig_gene_df$Age <- gsub("X", "", rownames(sig_gene_df))
colnames(sig_gene_df)[1] = "value"
sig_gene_df$Age = factor(sig_gene_df$Age,levels = sig_gene_df$Age)

sig_gene_df$group = rep('Metabolites',nrow(sig_gene_df))

# Generate plot
p <- ggplot(sig_gene_df, aes(x = Age, y = value,color = group,group = group)) +
  geom_line(size = 2) +
  scale_color_manual(values = '#DD5129FF') +
  theme_bw() +
  theme(legend.title = element_blank(),
        #plot.margin = margin(1, 1, 1, 1, "cm"),
        aspect.ratio = 0.7,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+  
  theme(legend.position = "right") +
  labs(title = "Num dems with q < 0.05")

ggsave(p, file = paste0(outDir,"/M_q0.05_b0.05_num_5centers_bucket20.pdf"), width = 12,height = 6)

save(sig_gene_df, sig_genes_list, file = paste0(outDir,'/DEswan_res_plot_win5_bucket20.Rdata'))
sig_gene_anno = list()

for (i in names(sig_genes_list)) {
  print(i)
  if (length(sig_genes_list[[i]])!=0) {
    crest_ann = class[class$metabolite_id %in% sig_genes_list[[i]],]
    crest_ann$lmBeta = res_wide_b[res_wide_b$variable %in% sig_genes_list[[i]],i]
    crest_ann$p.value = res_wide_p[res_wide_p$variable %in% sig_genes_list[[i]],i]
    crest_ann$p.adj = res_wide_q[res_wide_q$variable %in% sig_genes_list[[i]],i]
    crest_ann$status = ifelse(crest_ann$lmBeta>0.05,'Up','Down')
    sig_gene_anno[[i]] = crest_ann
  }
}
save(sig_gene_anno, file = paste0(outDir,'/Deswan_annoGeneList.Rdata'))
load(file = paste0(outDir,'/DEswan_res_plot_win5_bucket20.Rdata'),verbose = T)
load(file = paste0(outDir,'/Deswan_annoGeneList.Rdata'),verbose = T)

## FigS5d metabolomics DEswan DEM ====
temdir = outDir
Cre_names = c('X7','X67')
col_vector = structure(names = Cre_names,c('#CE9344FF', '#17486FFF'))

sel_crest = sig_gene_anno[Cre_names]

sel_crest_new = lapply(seq_along(sel_crest), function(x) {sel_crest[[x]]$crest_name = names(sel_crest)[[x]];
return(sel_crest[[x]])}
)

resultdata_all_long = do.call("rbind", sel_crest_new)
resultdata_all_long$crest_name = factor(resultdata_all_long$crest_name, levels = Cre_names)

# Generate plot
resultdata_beta_long = arrange(resultdata_all_long, crest_name, lmBeta)
resultdata_beta_long$label <- rep(FALSE, nrow(resultdata_beta_long ))
for (i in unique(resultdata_beta_long $crest_name)) {
  df <- resultdata_beta_long[which(resultdata_beta_long$crest_name == i),]
  indeces = c(as.vector(top_n(df,10,lmBeta)[['metabolite']]), as.vector(top_n(df,-10,lmBeta)[['metabolite']]))
  resultdata_beta_long [which(resultdata_beta_long$crest_name==i & resultdata_beta_long$metabolite%in%indeces),'label'] = TRUE
}
table(resultdata_beta_long$crest_name, resultdata_beta_long$label)
table(resultdata_beta_long$crest_name, resultdata_beta_long$status)
df_t = resultdata_beta_long[which(resultdata_all_long$label == TRUE), ]

pos <- position_jitter(width = 0.3, seed = 2)
p <- ggplot(resultdata_beta_long, aes(x = crest_name, y = lmBeta, color = crest_name,size = label)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-3, 2)) +
  geom_jitter(alpha = 0.8,position = pos) +
  scale_size_manual(values = c(1,3)) +
  theme_bw() +
  #coord_fixed(ratio = 0.8) +
  labs(title = "DE-SWAN at Peaks") +
  scale_color_manual(values = col_vector) +
  geom_text_repel(label = ifelse(resultdata_beta_long$label,resultdata_beta_long$metabolite,''),
                  inherit.aes = T,max.overlaps = 30,position = pos,
                  color = 'black', size = 3) +
  theme(legend.position = "right",aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank())
p
ggsave( paste0(temdir, '/Manhattan_beta_DESWAN.pdf'),p,width = 8,height = 7)


## FigS5d metabolomics DEswan Upsets ====
sets1 = lapply(sig_gene_anno, function(x) x$metabolite[which(x$status=='Up')])
sets2 = lapply(sig_gene_anno, function(x) x$metabolite[which(x$status=='Down')])
names(sets1) = paste0(names(sets1),'_Up')
names(sets2) = paste0(names(sets2),'_Down')

sets = c(sets1[paste0(Cre_names,'_Up')], sets2[paste0(Cre_names,'_Down')])
sets_le <- sets[ lapply(sets, length) > 0 ]
names(sets)
num_degs <- sapply(sets_le, length)
sets_le <- sets_le[order(-num_degs)]

# Define colors
color_vector <- mapvalues(
  names(sets_le),
  from = names(sets),
  to = c(col_vector, col_vector)
)

pdf(file = paste0(outDir, "/CombineUpDown_upset.pdf"))
print(
  upset(
    fromList(sets_le),
    nsets = length(sets_le),
    order.by = "freq",
    sets.bar.color = color_vector,
    text.scale = 2
  ),newpage = FALSE
)
dev.off()

# Write out gene lists
inter <- get.venn.partitions(sets_le)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[,colnames(inter) %!in% c("..set..", "..values..")], file = paste0(outDir,"/CombineUpDown_intersection.txt"), sep = '\t', quote = FALSE)

