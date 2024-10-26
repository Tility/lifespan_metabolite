setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(factoextra)
library(pvca)
library(Biobase)
library(stringr)
library(RColorBrewer)
library(paletteer)
library(openxlsx)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure1/'

if(!dir.exists(outDir)){
  dir.create(outDir)
}

## Fig1d metabolite cv QC =====
metab.class = read.csv('matabolite_CV/metabolite_CV_QC.csv',header=TRUE,row.names = 1)
metab.class$SuperClass[which(metab.class$SuperClass == "")] = 'Unclassified'
metab.class$Metabolite = rownames(metab.class)

class = read.xlsx('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/Matabolite_Class_heatmap/metab_class.xlsx',
                  startRow = 1,sheet = 1
                  )
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
#class_color = colorRampPalette(brewer.pal(12,"Paired")[c(1:4, 7:12)])(nlevels(class$Class))
##"#A6CEE3" "#408DBF" "#68AB9F" "#92CF72" "#33A02C" "#CAB75E" "#FE9F37" "#F18B35" "#CAB2D6" "#815AA8" "#B49E99" "#EBD57C" "#B15928"
#class$class_color = sapply(class$Class, function(x) class_color[[x]])

meta_cv = metab.class[,colnames(metab.class)%!in%c("SuperClass","Metabolite")]
CV = data.frame(Metabolite = rownames(metab.class),
                CV = apply(meta_cv, 1, cv_metab),
                Class = class$Class[match(rownames(metab.class), rownames(class))])
CV<- CV %>% 
  mutate(CV = CV*100)
CV$Class = factor(CV$Class, levels = levels(class$Class))

ggplot(CV, aes(x=Class, y=CV, fill=Class))+
  stat_boxplot(geom = "errorbar",
               width=0.3,
               size=0.4,
               position = position_dodge(0.8))+
  geom_boxplot( width = 0.7, outlier.shape = NA, 
                fatten = 1, ## middle line size
                position = position_dodge(0.8)) +  
  scale_fill_manual(values =class_color) +  
  geom_hline(yintercept = 30, linetype="dashed")+
  theme(panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white", colour = 'black', size = 0.7), 
        plot.title = element_text(size = 8, hjust = 0.5), plot.subtitle = element_text(size = 8, hjust = 0.5), 
        axis.text = element_text(size = 8, color = 'black'), aspect.ratio = 0.5, 
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0), axis.title = element_text(size = 8, color = 'black')) 
ggsave(paste0(outDir,"Metablite_CV_NewClass.pdf"), width = 10, height = 5)

table(CV$CV<30)
##FALSE  TRUE 
##349  1582 
1582/1931*100
##[1] 81.92646
save(metab.class,meta_cv,CV, class_color, file = paste0(outDir,'/metabolite_CV_RSD_newClass.Rdata'))
table(CV$CV<20)
##FALSE  TRUE 
##970   961 
write.xlsx(CV, file = paste0(outDir, '/Metabolite_CV.xlsx'))

## Fig1b metabolite class proportion =====
class = arrange(class, Class)
unique(class$Class)
class.ratio <- as.data.frame(table(class$Class))
colnames(class.ratio) = c('Class','value')
class.ratio<- class.ratio %>% 
  mutate(prop = value / sum(class.ratio$value) *100)
levels(class.ratio$Class)

class.ratio<- class.ratio %>% 
  arrange(desc(Class))%>% 
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(class.ratio, aes(x="", y=prop, fill=Class)) +
  geom_bar(stat="identity", width=1, color=NA) + 
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  geom_text(aes(y = ypos, 
                x = sum(prop)/60),  
            label = paste0(round(class.ratio$prop,2),'%'), 
            color = 'gray20', size=6)+
  scale_fill_manual(values = class_color )
ggsave(paste0(outDir,"Metablite_pie_SuperClass.pdf"), width = 10, height = 5)

## Fig1e Figs1c metabolite PCA =====
data = read.csv(file = "/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Metabolite_PCA_input.csv",header = FALSE,row.names = 1)
datat = as.data.frame(t(data))
colnames(datat)[1] = 'Metabolites'
datat1=as.data.frame(apply(datat[,-1],2,function(x) as.numeric(x))) 
rownames(datat1) = datat$Metabolites
phenotype = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Sample_phenotype.csv',header=TRUE,row.names = 1)
phenotype$sample = rownames(phenotype)
phenotype = phenotype[colnames(datat1),]
phenotype$Group = str_sub(rownames(phenotype),1,1)

## normalization by Sum 
## metab.t -- cols are sample, rows are metabolite
## Areai/Sum(Areai)
metab.t = datat1
metab.normal <- sweep(metab.t, 2, colSums(metab.t, na.rm=TRUE)/100, FUN="/")
dim(metab.normal)
write.csv(metab.normal,paste0(outDir,"Metabolite_norm.csv"))
metab.normal = read.csv(paste0(outDir,"Metabolite_norm.csv"),header = TRUE,row.names = 1)
metab.normal.log = log2(metab.normal)
metab.normal.log.t = t(metab.normal.log)
metab.normal.log.scal = apply(metab.normal.log.t, 2, AutoNorm)

## compute pca
# normalize to zero mean and unit variance
pca_res <- prcomp(metab.normal.log.scal, scale. = FALSE)
##point shape
mycol = paletteer_d("MetBrewer::Tiepolo")
shape = rep(19,7)
batch_col = c('#5773CCFF', '#FFB900FF')
sex_col = c('female' = '#DE7862FF', 'male' = '#5A6F80FF','unknown' = '#A9B4BCFF')

fviz_pca_ind(pca_res, col.ind = phenotype$Group,pointsize = 4,invisible="quali",geom = c("point"), 
             palette = mycol, addEllipses = TRUE,axes.linetype=NA,ellipse.level=0.99,
             ellipse.alpha = 0)+ 
  scale_shape_manual(values=shape)+ 
  #lims(x = c(-1000,1000),y = c(-1000,1000))+
  theme_bw()+
  theme(
    axis.text.y   = element_text(size=12),
    axis.text.x   = element_text(size=12),
    axis.title.y  = element_text(size=15),
    axis.title.x  = element_text(size=15),
    axis.ticks=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
  )

ggsave(paste0(outDir,"PCA_Group_addEllipses.pdf"),width = 7,height = 6)

fviz_pca_ind(pca_res, col.ind = phenotype$sex,pointsize = 4,invisible="quali",geom = c("point"),
             palette = sex_col, addEllipses = TRUE,axes.linetype=NA,ellipse.level=0.99,
             ellipse.alpha = 0)+ 
  scale_shape_manual(values=shape)+ ###控制点的形状
  theme_bw()+
  theme(
    axis.text.y   = element_text(size=12),
    axis.text.x   = element_text(size=12),
    axis.title.y  = element_text(size=15),
    axis.title.x  = element_text(size=15),
    axis.ticks=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
  )

ggsave(paste0(outDir,"PCA_sex_addEllipses.pdf"),width = 7,height = 6)

save(datat1,metab.normal,metab.normal.log,metab.normal.log.scal,pca_res,phenotype,
     file = paste0(outDir,'/PCA_result.Rdata'))

metab.normal.out = as.data.frame(metab.normal.log.scal)
metab.normal.out$Group = str_sub(rownames(metab.normal.out), 1, 1)
#write.xlsx(metab.normal.out,file = paste0(outDir, '/Metabolite_PCA_plot.xlsx'),rowNames = TRUE)

## Figs1b Gender distribution ratio =====
gender.data = as.data.frame(prop.table(table(phenotype$Group, phenotype$sex), margin = 1))
gender.data

colnames(gender.data)[1:3] = c('Group', 'Gender', 'prop')

ggplot(gender.data) + 
  geom_bar(aes(x =Group, y= prop, fill = Gender),stat = "identity",width = 1,size = 0.5,colour = 'white')+ 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 14, hjust = 0.5), 
        axis.text = element_text(size = 14, color = 'black'), aspect.ratio = 2, 
        axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5), axis.title = element_text(size = 14, color = 'black')) +
  labs(x='Age Stage',y = 'proportion')+
  scale_fill_manual(values = sex_col)+
  guides(fill=guide_legend(title=NULL))

ggsave(paste0(outDir,"GenderProp.pdf"), width = 7, height = 7)

## Figs1d metabolite pvca =====
set.seed(2024)
pvca_data = log10(metab.normal)
phenoData <- new("AnnotatedDataFrame",data=phenotype)
fd=as.data.frame(row.names(pvca_data))
colnames(fd)[1] <- "metabolite_name"
rownames(fd)<-row.names(pvca_data)

featureData = new("AnnotatedDataFrame",data=fd)
expr = as.matrix(pvca_data)
metab.set<- ExpressionSet(assayData = expr, phenoData = phenoData, featureData = featureData) #assayData :columns are samples and rows are variables
metab.set2 <- metab.set[ , metab.set$sex != "unknown"]
batch.factors <- c("Age", "sex")

pct_threshold <- 0.6
pvcaObj <- pvcaBatchAssess(metab.set2, batch.factors, pct_threshold)

plot.dat <- data.frame(eff=pvcaObj$label, prop=pvcaObj$dat[1,])
head(plot.dat)
plot.dat = plot.dat %>% arrange(desc(prop)) %>% mutate(prop = prop*100)
plot.dat$eff = factor(plot.dat$eff, levels = unique(plot.dat$eff))

ggdotchart(plot.dat, x = "eff", y = "prop",
           color = "orange",
           sorting = "descending",
           add = "segments",
           add.params = list(color = "lightgray", size = 1), 
           dot.size = 6,
           label = paste0(round(plot.dat$prop,2), '%'),                       
           font.label = list(color = "black", size = 14,
                             hjust = 0.5,
                             vjust = 0),              
           ggtheme = theme_pubr() 
           )+
  ylim(0,100)+
  ylab("Weighted average proportion variance")+ggtitle("Metabolism")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1), aspect.ratio = 1)

ggsave(filename = paste0(outDir,"pvac_metabolite_norm_lollipop.pdf"),width = 6,height = 6)

openxlsx::write.xlsx(plot.dat, file = paste0(outDir, '/PVCA_out.xlsx'),rowNames = TRUE)






