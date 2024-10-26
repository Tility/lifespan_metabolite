setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(plyr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(paletteer)
library(openxlsx)
library(clusterProfiler)
library("UpSetR")
library(VennDiagram)
library(ComplexHeatmap)
library(grid)
library(circlize)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure1/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig2a metabolite class mean expression Heatmap =====
class = read.xlsx('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/Matabolite_Class_heatmap/metab_class.xlsx',
                  startRow = 1,sheet = 1)
rownames(class) = class$metabolite
class = class[,-1]
colnames(class) = c('metabolite','Class','cpd','class_color') ##compound	cpd
class$Class = factor(class$Class, levels = c('Fatty Acyls', 'Steroids and steroid derivatives', 'Glycerophospholipids', 'Prenol lipids', 
                                             'Organic acids and derivatives', 'Unclassified', 'Organoheterocyclic compounds', 'Benzenoids', 
                                             'Organic oxygen compounds', 'Phenylpropanoids and polyketides', 'Organic nitrogen compounds', 
                                             'Nucleosides, nucleotides, and analogues', 'Others'
))

class_color = unique(class$class_color)
names(class_color) = levels(class$Class)

metab.normal = read.csv(paste0(dataDir,"Metabolite_norm.csv"),header = TRUE,row.names = 1)
metab.normal.log = log10(metab.normal)

metab_mean = as.data.frame(matrix(nrow=length(levels(class$Class)),ncol=ncol(metab.normal.log)))
rownames(metab_mean) = levels(class$Class)
colnames(metab_mean) = colnames(metab.normal.log)

for (i in levels(class$Class)) {
  metab_mean[i,] <- apply(metab.normal.log[class$metabolite[class$Class==i],], 2, mean)
}

phenotype = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Sample_phenotype.csv',header=TRUE,row.names = 1)
phenotype$sample = rownames(phenotype)
phenotype = phenotype[colnames(metab_mean),]
phenotype$Group = str_sub(rownames(phenotype),1,1)

##color
mycol = paletteer_d("MetBrewer::Tiepolo")
mycol2 = as.vector(mycol[1:7])

group_file = data.frame(sample = rownames(phenotype),
                        Group2 = phenotype$Group
)
rownames(group_file) = group_file$sample
group_file$col = group_file$Group2
Group = unique(group_file$Group2)

for (i in 1:length(unique(group_file$Group2))) {
  group_file$col[which(str_detect(group_file$Group2,Group[i]))] = mycol2[i]
}

##metab_mean cols are samples, rows are metabolics
metab_mean.scal = apply(metab_mean, 1, AutoNorm)
metab_mean2 = t(metab_mean.scal)
openxlsx::write.xlsx(as.data.frame(metab_mean2), file = paste0(outDir, '/Metabolite_Class_MeanExp.xlsx'),rowNames = TRUE)

#color
col_fun = colorRamp2(c(-2,-1, 0,1, 2), c("#3B66A2","#B9BFCC",'white',"#D9B4AE","#9E3D3D"))
column_ha = HeatmapAnnotation(
  Group = group_file$Group2,
  col = list(
    Group = structure(names = unique(group_file$Group2),
                      mycol2)
  ))

pdf(paste0(outDir, 'Metabolite_SuperClass_heatmap_rm.pdf'),width = unit(15, "cm"), height = unit(5, "cm"))
Heatmap(metab_mean2,
        top_annotation = column_ha,
        col = col_fun,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 2),
        column_names_rot = 0,
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = 12),
          title = " "
        ),
        width = 100, height = 50 
)
dev.off()

## Fig2b metabolite A vs Other groups common =====
data = read.csv(file = "/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Metabolite_PCA_input.csv",header = FALSE,row.names = 1)
datat = as.data.frame(t(data))
colnames(datat)[1] = 'Metabolites'
datat1=as.data.frame(apply(datat[,-1],2,function(x) as.numeric(x))) 
rownames(datat1) = datat$Metabolites
#life.metab = datat1
life.metab.group = data.frame(t(datat1),check.names = F)
life.metab.group$group = str_sub(rownames(life.metab.group), 1, 1)
table(life.metab.group$group)

life.metab.normal <- as.data.frame(t(metab.normal))
dim(life.metab.normal)
#[1]  136 1931
life.metab.normal.log <- log10(life.metab.normal)
ncol(life.metab.normal)

save(life.metab.group,life.metab.normal,file = paste0(outDir,'/Metabolite_expFor_DEM.Rdata'))

AgeG = unique(life.metab.group$group)
base_idents = AgeG
logFC = 0 
temDir = paste0(outDir,'/FC_',logFC,'_used/')

if (!dir.exists(temDir)){
  dir.create(temDir)}

for(base_ident in base_idents){
  print(base_ident)
  comp_idents = AgeG[AgeG %!in% base_ident]
  tem_dir = paste0(temDir, 'baseIdents_', base_ident)
  
  if (!dir.exists(tem_dir)){
    dir.create(tem_dir)}
  
  if (base_ident == 'A') {
    marker = TRUE ## A as numerator, others as denominator (AvsOthers)
  }
  else {
    marker = FALSE ## for example (Others vs B)
  }
  
  marker_list = list()
  DEG_Num = data.frame()
  set.seed(2023)
  
  for (i in 1:length(comp_idents)) {
    cat(base_ident, 'vs')
    print(comp_idents[i])
    marker_list[[i]] = FindMarker_metabolite(life.metab.group, life.metab.normal, life.metab.normal.log, 
                                             class = class,
                                             idents.1 = base_ident, idents.2 = comp_idents[i], marker = marker,
                                             LOG2FC = logFC, p_adj = 0.05,
                                             pAdjustMethod = "fdr", outDir = tem_dir, 
                                             return_data = TRUE)
    deg_up = marker_list[[i]][which(marker_list[[i]]$change_t.test == 'Up'),]
    deg_down = marker_list[[i]][which(marker_list[[i]]$change_t.test == 'Down'),]
    DEG_Num[i, 'Group'] = comp_idents[i]
    DEG_Num[i, 'Up'] = nrow(deg_up)
    DEG_Num[i, 'Down'] = nrow(deg_down)
    DEG_Num[i, 'ALL'] = nrow(deg_down)+nrow(deg_up)
  }
  
  names(marker_list) = comp_idents
  save(marker_list,DEG_Num, file = paste0(tem_dir, '/Pairwise_DEG.Rdata'))
  write.xlsx(DEG_Num,paste0(tem_dir, "/",base_ident,"vs_DEMnum.xlsx"))
  
  marker_num.long <-gather(DEG_Num,key=Change,value=Count,-Group)
  marker_num.long$Group<-factor(marker_num.long$Group)
  marker_num.long$Change<-factor(marker_num.long$Change, levels = c('ALL', 'Up', 'Down'))
  marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] <- marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] * -1 
  
  # Generate plot
  p <- ggplot(marker_num.long, aes(x = Group, y = Count,color = Change,group = Change)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c('#EA967C', '#DA4B35', '#3B5182')) +
    theme_bw() +
    theme(legend.title = element_blank(),
          aspect.ratio = 0.5,
          legend.text = element_text(size=14),
          axis.text = element_text(color="black",size=10),
          axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
          axis.title.y = element_text(color = "black",size=14),
          axis.title.x = element_text(color = "black",size=14))+  
    theme(legend.position = "right") +
    labs(title = "Num dems with q < 0.05")
  
  ggsave(file = paste0(tem_dir,"/DEG_Num_Metabolites_FC", logFC,".pdf"),width = 6,height = 6)
}


base_ident = 'A' ## A B 
tem_dir = paste0(temDir, 'baseIdents_', base_ident,'/')
comDir = paste0(tem_dir, '/KEGG_enrichment_common/')

if (!dir.exists(comDir)){
  dir.create(comDir)}
load(file = paste0(tem_dir, '/Pairwise_DEG.Rdata'), verbose = TRUE)

### kegg pathway annotion file =====
pathway_cpd_anno = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/KEGG_human/method3_dataset/KEGG_human_291_metab_pathway.csv',
                            header = TRUE)
pathway_cpd_anno$Description = substr(pathway_cpd_anno$Description, start = 1, stop = nchar(pathway_cpd_anno$Description)-23)
pathway_cpd_anno$CPD_ID = str_split(pathway_cpd_anno$CPD_ID,':',simplify = T)[,2]
path_metab = pathway_cpd_anno[,c(1,2)]
path_intro = pathway_cpd_anno[,c(1,4)]

save(class, pathway_cpd_anno, path_intro, path_metab, file = paste0(outDir,'/Metabolite_anno_ForKEGG.Rdata'))

### Read upset files =====
count_list = list()
data_keg = list()
DEM_list = lapply(marker_list, function(x) x[which(x$change_t.test != 'Not changed'),])

col_vector = as.vector(paletteer_d("MetBrewer::Tiepolo")[2:7])
sets1 = lapply(DEM_list, function(x) x$Metabolites[which(x$change_t.test=='Up')])
sets2 = lapply(DEM_list, function(x) x$Metabolites[which(x$change_t.test=='Down')])
names(sets1) = paste0(base_ident,'vs',names(marker_list),'_Up')
names(sets2) = paste0(base_ident,'vs',names(marker_list),'_Down')
sets = c(sets1, sets2)
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

# Create upset plot and export
pdf(file = paste0(comDir, "/All_upsetTop5.pdf"), width = 6, height = 14)
print(
  upset(
    fromList(sets_le),
    nsets = length(sets_le),
    mb.ratio = c(0.3, 0.7),
    nintersects = 5,
    order.by = "freq",
    sets.bar.color = color_vector,
    text.scale = 2
  )
)
dev.off()
print("Plot Done")

### Get intersection =======
for (status in c('Up','Down')) {
  sets = lapply(DEM_list, function(x) x$Metabolites[which(x$change_t.test==status)])
  sets_le <- sets[ lapply(sets, length) > 0 ]
  
  num_degs <- sapply(sets_le, length)
  sets_le <- sets_le[order(-num_degs)]
  
  inter <- get.venn.partitions(sets_le)
  
  count_list[[status]] = inter[1, '..values..'][[1]]
  
  data_keg[[status]] =  MarkerKEGG(count_list[[status]], path_metab = path_metab, path_intro = path_intro,
                                   minGSSize = 2, pAdjustMethod = "fdr",
                                   pvalueCutoff = 0.99, outDir = comDir, Clust = status,
                                   return_data = TRUE
  )
}

path_num = 10
dataMerge = as.data.frame(matrix(nrow=length(data_keg)*path_num,ncol=length(data_keg)+1)) 
colnames(dataMerge) = c('Description',names(data_keg))

dataMerge$Description = as.vector(sapply(data_keg, function(x) x$Description[1:path_num]))
dataMerge = dataMerge[!duplicated(dataMerge$Description),]
rownames(dataMerge) = dataMerge$Description
rownames(dataMerge) = factor(rownames(dataMerge))

for (i in names(data_keg)) {
  data_keg[[i]][,'-LogP'] = -log10(data_keg[[i]][,'p.adjust'])
  dataMerge[,i] = data_keg[[i]][,'-LogP'][match(rownames(dataMerge),data_keg[[i]]$Description)]
}
dataMerge[is.na(dataMerge)] <- 0

## Fig2c metabolite A vs Other groups common DEM KEGG ========
dataMerge.long <-gather(dataMerge,key=Change,value=LogPvalue,-Description)
for (i in names(data_keg)) {
  dataMerge.long[dataMerge.long$Change == i,'FoldEnrich'] = data_keg[[i]][,'FoldEnrich'][match(dataMerge.long$Description[dataMerge.long$Change == i],data_keg[[i]]$Description)]
}
dataMerge.long$Description <- factor(dataMerge.long$Description, levels = rev(rownames(dataMerge)))
dataMerge.long$Change <- factor(dataMerge.long$Change, levels = c('Up', 'Down'))
dataMerge.long[is.na(dataMerge.long)] <- 0

ggplot(dataMerge.long[dataMerge.long$LogPvalue!=0,], aes(Change, Description)) +
  geom_point(aes(fill = LogPvalue, size = FoldEnrich), color = "black", shape = 21) +
  scale_size(range = c(1, 6), breaks = c(6,10,14,18)) +
  scale_fill_gradientn(colours = colorRampPalette(rev(sequential_hcl(5, palette = "Reds")))(50),breaks = c(2, 4, 6,8)) +
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

ggsave(file = paste0(comDir,"/MarkerA_KEGG_top10_commonDotplot.pdf"),width = 6,height = 7)

save(data_keg,file = paste0(comDir,'/AvsOthersCommDEM_KEGG.Rdata'))

## Fig2d metabolite KEGG Pathway metabolite line Plot ========
metab_comm = metab.normal.log
metab_mean.scal = apply(metab_comm, 1, AutoNorm)
metab_mean.scal = as.data.frame(metab_mean.scal)
metab_mean.scal$group = str_sub(rownames(metab_mean.scal), 1, 1)
table(metab_mean.scal$group)
metab_age_long <- gather(metab_mean.scal, metab, expression, -group) ## reshape data

save(metab_mean.scal, metab_age_long, file = paste0(outDir,'/AllMetabolite_scaled.Rdata'))

### Up =======
keg_up = data_keg$Up
colman = c('#608B99','#807FA2','#B2D28D','#B19BBC','#7F5067',
           '#54629B','#E19B5E','#D06E5D','#B4DCD5','#F7DDB9',
           '#D97484','#E88A5E','#BCBCBA','#B36F9F'
)
pathys = keg_up$Description[1:10]

for (pathy in pathys) {
  keg_up_metab = str_split(keg_up$geneID[which(keg_up$Description == pathy)], pattern = '\\/', simplify = TRUE)[1,]
  
  metb_need = DEM_list$B[which(DEM_list$B$cpd %in% keg_up_metab & DEM_list$B$change_t.test == 'Up'),]
  metb_need = top_n(metb_need,n = 10,LOG2FC)
  
  data_sub = metab_age_long[metab_age_long$metab %in% metb_need$Metabolite, ]
  data_sub$metab = factor(data_sub$metab, levels = unique(metb_need$Metabolite))
  
  color_v = colorRampPalette(colman)(length(metb_need$Metabolite))
  
  ggplot(data_sub, aes(x=group, y=expression, fill=metab))+
    scale_fill_manual(values = color_v) +  
    stat_summary(fun=mean, 
                 geom="point", 
                 aes(color = metab))+
    stat_summary(fun=mean, 
                 aes(color = metab,group = metab),
                 geom="line")+
    scale_color_manual(values = color_v)+
    labs(x= "", y= "Mean Expression",  title = pathy)+
    theme(panel.grid = element_blank(), 
          legend.position = 'right',
          panel.background =  element_rect(fill = "white", colour = 'black', linewidth = 0.7), 
          plot.title = element_text(size = 15, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
          axis.text = element_text(size = 15, color = 'black'), aspect.ratio = 1, 
          axis.text.x = element_text(angle=0), axis.title = element_text(size = 15, color = 'black')) 
  ggsave(file = paste0(comDir,"/Kegg_up_",pathy,"_scaled_top10.pdf"), width = 6, height = 5)
  
}



### Down =======
keg_down = data_keg$Down
colman = c('#608B99','#807FA2','#B2D28D','#B19BBC','#7F5067',
           '#54629B','#E19B5E','#D06E5D','#B4DCD5','#F7DDB9',
           '#D97484','#E88A5E','#BCBCBA','#B36F9F'
)
pathys = keg_down$Description[1:10]

for (pathy in pathys) {
  keg_up_metab = str_split(keg_down$geneID[which(keg_down$Description == pathy)], pattern = '\\/', simplify = TRUE)[1,]
  
  if (pathy == 'Glycolysis / Gluconeogenesis') {
    pathy = 'Glycolysis_Gluconeogenesis'
  }
  
  metb_need = DEM_list$B[which(DEM_list$B$cpd %in% keg_up_metab & DEM_list$B$change_t.test == 'Down'),]
  metb_need = top_n(metb_need,n = -10,LOG2FC)
  
  data_sub = metab_age_long[metab_age_long$metab %in% metb_need$Metabolite, ]
  data_sub$metab = factor(data_sub$metab, levels = rev(unique(metb_need$Metabolite)))
  
  color_v = colorRampPalette(colman)(length(metb_need$Metabolite))
  
  ggplot(data_sub, aes(x=group, y=expression, fill=metab))+
    scale_fill_manual(values = color_v) +  
    stat_summary(fun=mean, 
                 geom="point", 
                 aes(color = metab))+
    stat_summary(fun=mean, 
                 aes(color = metab,group = metab),
                 geom="line")+
    scale_color_manual(values = color_v)+
    labs(x= "", y= "Mean Expression",  title = pathy)+
    theme(panel.grid = element_blank(), 
          legend.position = 'right',
          panel.background =  element_rect(fill = "white", colour = 'black', linewidth = 0.7), 
          plot.title = element_text(size = 15, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
          axis.text = element_text(size = 15, color = 'black'), aspect.ratio = 1, 
          axis.text.x = element_text(angle=0), axis.title = element_text(size = 15, color = 'black')) 
  ggsave(file = paste0(comDir,"Kegg_down_",pathy,"_scaled_top10.pdf"), width = 6, height = 5)
  
}
