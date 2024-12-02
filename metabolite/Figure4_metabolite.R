setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)
library(ggrepel)
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(paletteer)
library(openxlsx)
library(car)
library(lsr)
library(clusterProfiler)
library(ComplexHeatmap)
library(grid)
library(circlize)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure4/'
mycol = paletteer_d("MetBrewer::Tiepolo")
mycol2 = as.vector(mycol[1:7])

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig4a metabolomics linear model =====
load(file = paste0(dataDir,'/Metabolite_anno_ForKEGG.Rdata'),verbose = T)
class_color = unique(class$class_color)
names(class_color) = levels(class$Class)


load(file = paste0(dataDir,'/Metabolite_expFor_DEM.Rdata'),verbose = T)
metab.normal = as.data.frame(t(life.metab.normal))
metab.normal.log = log10(metab.normal)
metab.normal.log$Metabolite = rownames(metab.normal.log)
metab.normal.log$Metab_ID = paste0('Metab_', as.numeric(as.factor(metab.normal.log$Metabolite)) )
data = as.data.frame(t(metab.normal.log[,colnames(metab.normal.log)%!in%c("Metab_ID","Metabolite")]))
colnames(data) = metab.normal.log$Metab_ID
data$group = str_sub(rownames(data), 1, 1)

phenotype = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Sample_phenotype.csv',header=TRUE,row.names = 1)
phenotype$sample = rownames(phenotype)
phenotype$Group = str_sub(rownames(phenotype),1,1)

data$Gender <- phenotype$sex[match(rownames(data),rownames(phenotype))]
data$Age <- phenotype$Age[match(rownames(data),rownames(phenotype))]

metabdata = data[which(data$group != 'A'),]
table(metabdata$group)

age_stats_df <- data.frame(matrix(ncol = 3, nrow = 0))
sex_stats_df <- data.frame(matrix(ncol = 3, nrow = 0))

for (metab in colnames(metabdata)) {
  if (metab == "Age" |
      metab == "Gender" |
      metab == "group") {
    next
  }
  
  tmp <- metabdata[, c(metab, "Age", "Gender", "group")]
  
  lm <- lm(reformulate(c("Age", "Gender"), metab), data = tmp)
  
  # Run ANOVA
  anova <- Anova(lm, type = "II")
  
  pval <- anova["Age", "Pr(>F)"]
  beta <- coef(summary(lm))["Age", "Estimate"]
  df <- data.frame(metab, beta, pval)
  
  age_stats_df <- rbind(age_stats_df, df)
  
  pval <- anova["Gender", "Pr(>F)"]
  beta <- coef(summary(lm))["Gendermale", "Estimate"]
  df <- data.frame(metab, beta, pval)
  
  sex_stats_df <- rbind(sex_stats_df, df)
  
}

p_adj = 0.05
LOG2FC = 0

colnames(age_stats_df) <- c("Metab_ID", "avg_log2FC", "pval")
rownames(age_stats_df) <- age_stats_df$Metab_ID

# BH correct p values
age_stats_df$p.adj <- p.adjust(age_stats_df$pval, method = "BH")
age_stats_df$Metabolite = metab.normal.log$Metabolite[match(rownames(age_stats_df),metab.normal.log$Metab_ID)]
age_stats_df$state <- ifelse(age_stats_df$p.adj < p_adj & abs(age_stats_df$avg_log2FC) >= LOG2FC, 
                             ifelse(age_stats_df$avg_log2FC > LOG2FC ,'Up','Down'), "Not changed")
table(age_stats_df$state)

ggplot(age_stats_df, 
       aes(x = avg_log2FC, y = -log10(p.adj), colour=state)) +
  geom_point(size=3,alpha = 0.8, shape = 19,stroke = 0) +
  scale_color_manual(values=c('Down' = '#3681B7', 'Up' = '#C65248', 'Not changed' = "#d2dae2"))+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8) +
  labs(x="LM beta",y="-log10(p.adj)")+
  theme_bw(  
  )+
  xlim(-0.04,0.04) +
  theme(aspect.ratio = 1,
        axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=15),
        axis.title.x  = element_text(size=15),
        axis.ticks=element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
ggsave(paste0(outDir, '/Age_pvalue_volcano.pdf'),width = 7,height = 6)

colnames(sex_stats_df) <- c("Metab_ID", "avg_log2FC", "pval")
rownames(sex_stats_df) <- sex_stats_df$Metab_ID

# BH correct p values
sex_stats_df$p.adj <- p.adjust(sex_stats_df$pval, method = "BH")
sex_stats_df$Metabolite = metab.normal.log$Metabolite[match(rownames(sex_stats_df),metab.normal.log$Metab_ID)]
sex_stats_df$state <- ifelse(sex_stats_df$p.adj < p_adj & abs(sex_stats_df$avg_log2FC) >= LOG2FC, 
                             ifelse(sex_stats_df$avg_log2FC > LOG2FC ,'Up','Down'), "Not changed")

age_sex_comm = intersect(sex_stats_df$Metabolite[sex_stats_df$state != 'Not changed'], age_stats_df$Metabolite[age_stats_df$state != 'Not changed'])
length(age_sex_comm)

age_stats_df$Class = class$Class[match(age_stats_df$Metabolite,class$metabolite)]
age_stats_df$cpd = class$cpd[match(age_stats_df$Metabolite,class$metabolite)]

##### remove age-sex metabolites ===========
age_only_dem = age_stats_df[which(age_stats_df$Metabolite %!in% age_sex_comm & age_stats_df$state != 'Not changed'), ]

save(age_only_dem, file = paste0(outDir, '/DEM_signif_nocomm.Rdata'))
#load(paste0(outDir, '/DEM_signif_nocomm.Rdata'), verbose = TRUE)

#### write out DEM table ======
write.xlsx(age_only_dem, file = paste0(outDir, '/LM_Age_DEM_anno.xlsx'))

### linear DEMs line loess Age ####
##### Convert wide to long add Age Group ####
metab_age = metab.normal.log[age_only_dem$Metabolite,rownames(phenotype)]
metab_mean.scal = apply(metab_age, 1, AutoNorm)
metab_mean.scal = as.data.frame(metab_mean.scal)

metab_mean.scal$group = str_sub(rownames(metab_mean.scal), 1, 1)
table(metab_mean.scal$group)

metab_age_long <- gather(metab_mean.scal[which(metab_mean.scal$group != 'A'),], metab, expression, -group) ## reshape data
table(metab_age_long$group)

color_vet = c('Down' = '#3681B7', 'Up' = '#C65248')
for (i in c('Up', 'Down')) {
  
  rows <- age_only_dem$Metabolite[age_only_dem$state == i]
  
  data_sub = metab_age_long[metab_age_long$metab %in% rows, ]
  
  line_plot_of_group(data_sub, alpha = 0.3, color = color_vet[[i]],line_Tpos = NULL)
  
  ggsave(file = paste0(outDir,"metabolite_", i, "_trajectories_scaled.pdf"), width = 5, height = 4)
  
}

save(metab_mean.scal,metab_age_long, age_only_dem, file = paste0(outDir, '/linear_lineByGroup.Rdata'))


## Fig4c metabolomics linear model KEGG enrichment =====
data_keg = list()
for (i in c('Up', 'Down')) {
  data_keg[[i]] =  MarkerKEGG(age_only_dem$Metabolite[age_only_dem$state == i], path_metab = path_metab, path_intro = path_intro,
                              minGSSize = 2, pAdjustMethod = "fdr",
                              pvalueCutoff = 0.99, outDir = outDir, Clust = i,
                              return_data = TRUE
  )
}

names(data_keg) = paste0('Age_', names(data_keg))

save(data_keg, file = paste0(outDir, '/Age_metabKEGG.Rdata'))

path_num = 10
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
  scale_size(range = c(1, 8), breaks = c(5,15,30,45)) +
  scale_fill_gradientn(colours = colorRampPalette(rev(sequential_hcl(5, palette = "Reds")))(50),breaks = c(2, 4, 6)) +
  theme_bw() +
  coord_fixed(ratio = 0.4) +
  labs(x = "", y = "", title = "", 
       fill = "-Log10 (p.adjust)", size = "Fold Enrichment") +
  theme(#aspect.ratio = 1.6,
    axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
    axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
  ) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black")
  )

ggsave(file = paste0(outDir,"/KEGG_Age_top10Dotplot.pdf"),width = 6,height = 7)

