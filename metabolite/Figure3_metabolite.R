setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(plyr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(colorspace)
library(paletteer)
library(openxlsx)
library(clusterProfiler)
library(ComplexHeatmap)
library(grid)
library(circlize)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure3/'
mycol = paletteer_d("MetBrewer::Tiepolo")
mycol2 = as.vector(mycol[1:7])

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig3a Others Group vs B  rmA =====
base_ident='B'
tem_dir = paste0(dataDir, 'FC_0_used/baseIdents_', base_ident)
load(file = paste0(tem_dir, '/Pairwise_DEG.Rdata'),verbose = T)

marker_num.long <-gather(DEG_Num,key=Change,value=Count,-Group)
marker_num.long$Group<-factor(marker_num.long$Group)
marker_num.long$Change<-factor(marker_num.long$Change, levels = c('ALL', 'Up', 'Down'))
marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] <- marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] * -1 

# Generate plot
ggplot(marker_num.long[marker_num.long$Group!='A',], aes(x = Group, y = Count,color = Change,group = Change)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('#EA967C', '#DA4B35', '#3B5182')) +
  theme_bw() +
  theme(legend.title = element_blank(),
        #plot.margin = margin(1, 1, 1, 1, "cm"),
        aspect.ratio = 0.5,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+  
  theme(legend.position = "right") +
  labs(title = "Num dems with q < 0.05")

ggsave(file = paste0(tem_dir,"/DEG_Num_Metabolites_FC0_rmA.pdf"),width = 6,height = 6)

## Fig3b Compare adjacent groups =====
load(file = paste0(dataDir,'/Metabolite_expFor_DEM.Rdata'),verbose = T)
load(file = paste0(dataDir,'/Metabolite_anno_ForKEGG.Rdata'),verbose = T)
life.metab.normal.log <- log10(life.metab.normal)
AgeG = unique(life.metab.group$group)

marker_list = list()
DEG_Num = data.frame()
set.seed(2023)

base_idents = AgeG
logFC = 0 
tem_dir = paste0(outDir, '/logFC_', logFC,'_adjacentgroups')

if (!dir.exists(tem_dir)){
  dir.create(tem_dir)}

for (i in 2:length(base_idents)) {
  
  print(base_idents[i])
  
  marker_list[[paste0(base_idents[i], 'vs',base_idents[i-1])]] = FindMarker_metabolite(life.metab.group, life.metab.normal, life.metab.normal.log, 
                                                                                       class = class,
                                                                                       idents.1 = base_idents[i], idents.2 = base_idents[i-1],
                                                                                       marker = TRUE,
                                                                                       LOG2FC = logFC, p_adj = 0.05,
                                                                                       pAdjustMethod = "fdr", outDir = tem_dir, 
                                                                                       return_data = TRUE)
  deg_up = marker_list[[paste0(base_idents[i], 'vs',base_idents[i-1])]][which(marker_list[[paste0(base_idents[i], 'vs',base_idents[i-1])]]$change_t.test == 'Up'),]
  deg_down = marker_list[[paste0(base_idents[i], 'vs',base_idents[i-1])]][which(marker_list[[paste0(base_idents[i], 'vs',base_idents[i-1])]]$change_t.test == 'Down'),]
  DEG_Num[i-1, 'Group'] = paste0(base_idents[i], 'vs',base_idents[i-1])
  DEG_Num[i-1, 'Up'] = nrow(deg_up)
  DEG_Num[i-1, 'Down'] = nrow(deg_down)
  DEG_Num[i-1, 'ALL'] = nrow(deg_down)+nrow(deg_up)
}
names(marker_list)
save(marker_list,DEG_Num, file = paste0(tem_dir, '/PairwiseAdjacent_DEM.Rdata'))
write.xlsx(DEG_Num,paste0(tem_dir, "/PairwiseAdjacent_DEMnum.xlsx"))

load(file = paste0(tem_dir, '/PairwiseAdjacent_DEM.Rdata'))

marker_num.long <-gather(DEG_Num,key=Change,value=Count,-Group)
marker_num.long$Group<-factor(marker_num.long$Group)
marker_num.long$Change<-factor(marker_num.long$Change, levels = c('ALL', 'Up', 'Down'))
marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] <- marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] * -1 

# Generate plot
ggplot(marker_num.long[marker_num.long$Group!='BvsA',], aes(x = Group, y = Count,color = Change,group = Change)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('#EA967C', '#DA4B35', '#3B5182')) +
  theme_bw() +
  theme(legend.title = element_blank(),
        #plot.margin = margin(1, 1, 1, 1, "cm"),
        aspect.ratio = 0.5,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+  
  theme(legend.position = "right") +
  labs(title = "Num dems with q < 0.05")

ggsave(file = paste0(tem_dir,"/DEG_Num_Metabolites_logFC", logFC,"_rmA.pdf"),width = 6,height = 6)

## Fig3c Marker metabolites of Age Group =====
marker_list = list()
DEG_Num = data.frame()
set.seed(2023)

for (i in 1:length(AgeG)) {
  print(AgeG[i])
  marker_list[[AgeG[i]]] = FindMarker_metabolite(life.metab.group, life.metab.normal, life.metab.normal.log, 
                                           class = class,
                                           idents.1 = AgeG[i],
                                           LOG2FC = 1, p_adj = 0.05,
                                           pAdjustMethod = "fdr", outDir = outDir, 
                                           return_data = TRUE)
  deg_up = marker_list[[AgeG[i]]][which(marker_list[[AgeG[i]]]$change_t.test == 'Up'),]
  deg_down = marker_list[[AgeG[i]]][which(marker_list[[AgeG[i]]]$change_t.test == 'Down'),]
  DEG_Num[i, 'Group'] = AgeG[i]
  DEG_Num[i, 'Up'] = nrow(deg_up)
  DEG_Num[i, 'Down'] = nrow(deg_down)
  DEG_Num[i, 'ALL'] = nrow(deg_down)+nrow(deg_up)
}

save(marker_list,class, DEG_Num,file = paste0(outDir, 'Metabolite_Group_Marker.Rdata'))

marker_num.long <-gather(DEG_Num,key=Change,value=Count,-Group)
marker_num.long$Group<-factor(marker_num.long$Group)
marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] <- marker_num.long[which(marker_num.long$Change == 'Down'), c('Count')] * -1 

marker_num.long_sub = marker_num.long[marker_num.long$Change =='Up'&marker_num.long$Group!='A',]

ggplot(marker_num.long[marker_num.long$Change =='Up'&marker_num.long$Group!='A',],aes(Group,Count,fill=Group))+
  geom_col(width=0.5)+
  scale_y_continuous(breaks = seq(0, 300,50), 
                     labels = as.character(abs(seq(0, 300,50))), #设置y轴合适范围，并将负数转化为正数
                     limits = c(0, 300))+  
  scale_fill_manual(values =  structure(names = AgeG,mycol2))+
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        aspect.ratio = 1,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+
  ylab("Number of Marker Metabolites(FC)>2")

ggsave(file = paste0(outDir,"/Num_Marker_Metabolites_FC2_rmA.pdf"),width = 6,height = 6)

## Fig3e Marker gene class proportion ===== 
marker.fc.lst = lapply(seq_along(marker_list), function(x){
  marker_list[[x]]$Group = names(marker_list)[x]
  marker_list[[x]][marker_list[[x]]$change_t.test != 'Not changed',c('Metabolites','LOG2FC','t.test_BH','Group'), ] 
})

marker.fc = do.call('rbind',marker.fc.lst)
marker.fc$Group = paste0('MarkerOf',marker.fc$Group)
table(is.na(marker.fc))

marker.fcdf = marker.fc
colnames(marker.fcdf)[1] = 'metabolite'
marker.fcdf = merge(marker.fcdf, class, by = 'metabolite', all.x = T)
marker.fcdf$change = ifelse(marker.fcdf$LOG2FC > 0, 'Up', 'Down')
marker.fcdf$Group_change = paste0(marker.fcdf$Group, '_',marker.fcdf$change)
save(marker.fcdf, class, file = paste0(outDir, '/Marker_signif_FC.Rdata'))

class_color = unique(class$class_color)
names(class_color) = levels(class$Class)

class.ratio <- as.data.frame(prop.table(table(marker.fcdf$Group_change, marker.fcdf$Class), margin = 1))
colnames(class.ratio) = c('Group_change','Class', 'prop')
class.ratio$Class = factor(class.ratio$Class, levels = levels(class$Class))
class.ratio$Group_change = factor(class.ratio$Group_change, levels = c(paste0('MarkerOf',AgeG, '_', 'Up'),paste0('MarkerOf',AgeG, '_', 'Down'))
)
ggplot(class.ratio[which(str_detect(class.ratio$Group_change,'_Up')&str_detect(class.ratio$Group_change,'A_',negate = T)),]) + 
  geom_bar(aes(x =Group_change, y= prop, fill = Class),stat = "identity",width = 0.9,size = 0.5,colour = 'white')+ 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 14, hjust = 0.5), 
        axis.text = element_text(size = 14, color = 'black'), aspect.ratio = 2, 
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5), axis.title = element_text(size = 14, color = 'black')) +
  labs(x='Group_change',y = 'proportion')+
  scale_fill_manual(values = class_color)+
  guides(fill=guide_legend(title=NULL))
ggsave(paste0(outDir,"/Marker_stackedbar_Class_up_rmA.pdf"), width = 7, height = 7)


## Fig3d Marker gene line plot of each Age Group ===== 
marker.fcdf.up = marker.fcdf[marker.fcdf$change=='Up',]
metab_comm_all = life.metab.normal.log[, marker.fcdf.up$metabolite]
metab_mean.scal_all = apply(metab_comm_all, 2, AutoNorm)
metab_mean.scal_all = as.data.frame(metab_mean.scal_all)
metab_mean.scal_all$group = str_sub(rownames(metab_mean.scal_all), 1, 1)
table(metab_mean.scal_all$group)
metab_age_long_all <- gather(metab_mean.scal_all, metab, expression, -group) ## reshape data
color_vet = structure(names = unique(marker.fcdf.up$Group), mycol2)

for (i in unique(marker.fcdf.up$Group)) {
  
  rows <- marker.fcdf.up$metabolite[marker.fcdf.up$Group == i]
  
  data_sub = metab_age_long_all[metab_age_long_all$metab %in% rows, ]
  
  line_plot_of_group(data_sub, alpha = 0.8, color = color_vet[[i]],line_Tpos = NULL)
  
  ggsave(file = paste0(outDir,"/", i, "_trajectories_scaled.pdf"), width = 5, height = 4)
  
}

## Fig3f Marker gene KEGG pathway ===== 
data_keg = list()
for (i in 1:length(AgeG)) {
  data_keg[[AgeG[i]]] =  MarkerKEGG(marker.fcdf.up$metabolite[marker.fcdf.up$Group ==paste0('MarkerOf', AgeG[i])], path_metab = path_metab, path_intro = path_intro,
                              minGSSize = 2, pAdjustMethod = "fdr",
                              pvalueCutoff = 0.99, outDir = outDir, Clust = AgeG[i],
                              return_data = TRUE
  )
}
save(data_keg,file = paste0(outDir,'/Marker_KEGGlist.Rdata'))
data_keg_r = data_keg[2:7]
path_num = 3
dataMerge = as.data.frame(matrix(nrow=length(data_keg_r)*path_num,ncol=length(data_keg_r)+1)) ## 
colnames(dataMerge) = c('Description',names(data_keg_r))

dataMerge$Description = as.vector(sapply(data_keg_r, function(x) x$Description[1:path_num]))
dataMerge = dataMerge[!duplicated(dataMerge$Description),]
dataMerge = dataMerge[!is.na(dataMerge$Description),]
rownames(dataMerge) = dataMerge$Description

for (i in names(data_keg_r)) {
  data_keg_r[[i]][,'-LogP'] = -log10(data_keg_r[[i]][,'p.adjust'])
  dataMerge[,i] = data_keg_r[[i]][,'-LogP'][match(rownames(dataMerge),data_keg_r[[i]]$Description)]
}

dataMerge[is.na(dataMerge)] <- 0

#### Dot plot ========
dataMerge.long <-gather(dataMerge,key=Change,value=LogPvalue,-Description)
for (i in names(data_keg_r)) {
  dataMerge.long[dataMerge.long$Change == i,'FoldEnrich'] = data_keg_r[[i]][,'FoldEnrich'][match(dataMerge.long$Description[dataMerge.long$Change == i],data_keg_r[[i]]$Description)]
}
dataMerge.long$Description <- factor(dataMerge.long$Description, levels = rev(rownames(dataMerge)))
dataMerge.long$Change <- factor(dataMerge.long$Change, levels = names(data_keg_r))

ggplot(dataMerge.long, aes(Change, Description)) +
  geom_point(aes(fill = LogPvalue, size = FoldEnrich), color = "black", shape = 21) +
  scale_size(range = c(2, 8), breaks = c(5,15,50,75,200)) +
  scale_fill_gradientn(colours = colorRampPalette(rev(sequential_hcl(5, palette = "Reds")))(50),breaks = c(2, 4, 6)) +
  theme_bw() +
  coord_fixed(ratio = 1.8) +
  labs(x = "", y = "", title = "", 
       fill = "-Log10 (p.adjust)", size = "Fold Enrichment") +
  theme(#aspect.ratio = 1.6,
    axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
    axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
  ) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black")
  )

ggsave(file = paste0(outDir,"/KEGG_Marker_top3Dotplot_All_rmA.pdf"),width = 6,height = 7)


## Fig3g select Marker gene Heatmap ===== 
marker_top = marker.fcdf %>%
  group_by(Group) %>%
  arrange(Group, t.test_BH,desc(LOG2FC))%>%
  dplyr::filter(LOG2FC > 1)

write.xlsx(as.data.frame(marker_top), file = paste0(outDir, '/Maker_p.adjFC2_byGroup_signif.xlsx'))

marker.count = read.xlsx('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure3/Maker_select_1.xlsx',
                         startRow = 1,sheet = 1)
marker.count = marker.count[,-1]

phenotype = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Sample_phenotype.csv',header=TRUE,row.names = 1)
phenotype$sample = rownames(phenotype)
phenotype$Group = str_sub(rownames(phenotype),1,1)

group_file = data.frame(sample = rownames(phenotype),
                        Age = phenotype$Age,
                        Group2 = phenotype$Group
)
rownames(group_file) = group_file$sample
group_file$col = group_file$Group2
Group = unique(group_file$Group2)

for (i in 1:length(unique(group_file$Group2))) {
  group_file$col[which(str_detect(group_file$Group2,Group[i]))] = mycol2[i]
}
### rmA marker Heatmap ====
group_file = group_file[group_file$Group2!='A',]
metab_comm = life.metab.normal.log[group_file$sample, marker.count$metabolite]
metab_mean.scal = apply(metab_comm, 2, AutoNorm)
metab_comm2 = t(metab_mean.scal)

openxlsx::write.xlsx(as.data.frame(metab_comm2),file = paste0(outDir, '/Marker_select_metabolite_scaled.xlsx'),rowNames = TRUE)

#color
col_fun = colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#3B66A2","#B9BFCC",'white',"#D9B4AE","#9E3D3D"))

column_ha = HeatmapAnnotation(
  Age = anno_barplot(group_file$Age, 
                     gp = gpar(col = structure(names = group_file$Age,
                                               group_file$col))
  ),
  Group = group_file$Group2,
  col = list(
    Group = structure(names = unique(group_file$Group2),
                      mycol2[2:7])
  ))


pdf(paste0(outDir, '/Marker_common_rmA.pdf'),width = unit(10, "cm"), height = unit(7, "cm"))
Heatmap(metab_comm2,
        top_annotation = column_ha,#annotation
        col = col_fun,#color
        cluster_rows = FALSE,
        cluster_columns = FALSE,# turn off row clustering
        row_names_gp = gpar(fontsize = 8),#rowname size
        column_names_gp = gpar(fontsize = 2),#columnname size
        column_names_rot = 0,#columnname angle
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = 12),
          title = " "
        ),
        width = 50, height = 50 #heatmap size
)
dev.off()

