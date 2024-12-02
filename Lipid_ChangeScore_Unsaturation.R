##Most of the lipid analyses were consistent with the metabolome, and only metabolome data needed to be replaced with lipidome data
setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/lipid/')
library(plyr)
library(tidyverse)
library(ggthemes)
library(openxlsx)
library(stringr)
library(colorspace)
library(paletteer)

source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/lipid/1_lipid_figure/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/lipid/1_lipid_figure/LipidFigure/'
mycol = paletteer_d("MetBrewer::Tiepolo")
mycol2 = as.vector(mycol[1:7])

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Change Score =====
load(file = paste0(dataDir, '/DEM_signif_anno.Rdata'), verbose = T)
class.level <- as.data.frame(prop.table(table(age_only_dem$state, age_only_dem$Class), margin = 1))
colnames(class.level) = c('Change','Class', 'prop')
class.level$Change = factor(class.level$Change, levels = c('Up','Down') )

class.level$count <- as.data.frame(table(age_only_dem$state, age_only_dem$Class))$Freq
class.level.wide = spread(class.level[,colnames(class.level)!='prop'], key = Change, value = count )

class.level.wide$change.radio = apply(class.level.wide[,2:3], 1,
                                      function(x) (x[1]-x[2])/sum(x))
class.level.wide = arrange(class.level.wide, desc(change.radio))
class.level.wide$Class = factor(class.level.wide$Class, levels = unique(class.level.wide$Class))
class.level.wide$ALL = apply(class.level.wide[,2:3], 1,
                             function(x) sum(x))

ggplot(class.level.wide, aes(x=change.radio, y=ALL)) +
  geom_point(aes(color=Class,size=ALL),shape = 19) +
  scale_color_manual(values = class_color)+
  scale_size(range = c(2, 6), breaks = c(2, 75, 150, 300),name = 'Size') + 
  geom_text(
    label =class.level.wide$Class,
    size = 14/.pt,
    hjust = 0.5,
    vjust = 0,
    color = "black")+
  geom_vline(xintercept = 0,linetype = 2)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 14,angle = 90,hjust = 1,vjust =0.5, color = 'black'), 
        text=element_text(size = 14, color = 'black'),
        aspect.ratio = 1)+
  xlim(-1,1) +
  labs(x= "Change Score", y= "Size", title = 'Lipidomics')

ggsave(paste0(outDir,"/Lipid_Class_ChangeScore.pdf"), width = 9, height = 7)

## Total expression of metabolites at different saturation levels =====
class_color
class$Class = factor(class$Class,levels = names(class_color))
class$unsaturation_3Class = ifelse(class$unsaturation_new == 0, 'Saturated', ifelse(class$unsaturation_new >1, 'Polyunsaturated', 'Monounsaturated'))

lipid.exp = read.csv(paste0(dataDir,"Lipid_norm.csv"),header = TRUE,row.names = 1)
lipid.exp$lipid = rownames(lipid.exp)
lipid.exp = merge(class, lipid.exp, by = 'lipid')

lipid.exp$length_unsaturation_PUMU_Class  = paste0(lipid.exp$length, '_',lipid.exp$unsaturation_3Class, '_',lipid.exp$Class)
table(lipid.exp$length_unsaturation_PUMU_Class)
not_expsscol = c(colnames(class),'length_unsaturation_PUMU_Class')

lipid.exp.sum.PUMU <- aggregate(apply(lipid.exp[,colnames(lipid.exp) %!in% not_expsscol],2,function(x) as.numeric(x)),by=list(lipid.exp$length_unsaturation_PUMU_Class),sum,na.rm=TRUE)
row.names(lipid.exp.sum.PUMU)<-lipid.exp.sum.PUMU$Group.1

leng_unsat = as.data.frame(table(lipid.exp$length_unsaturation_PUMU_Class))
lipid.exp.sum.PUMU$Type_count <- leng_unsat$Freq[match(lipid.exp.sum.PUMU$Group.1, leng_unsat$Var1)]
lipid.exp.sum.PUMU$length <- as.numeric(str_split(lipid.exp.sum.PUMU$Group.1, '_', simplify = T)[,1])
lipid.exp.sum.PUMU$unsaturation_PUMU <- str_split(lipid.exp.sum.PUMU$Group.1, '_', simplify = T)[,2]
lipid.exp.sum.PUMU$Class <- str_split(lipid.exp.sum.PUMU$Group.1, '_', simplify = T)[,3]

AgeG = c('All','A','B','C','D','E','F','G')
tem_dir = paste0(outDir, '/Class_length_unsaturation3PUMU/')
if (!dir.exists(tem_dir)){
  dir.create(tem_dir)}

lipid.exp.sum.PUMU.plot = lipid.exp.sum.PUMU

lipid.exp.sum.PUMU.plot = arrange(lipid.exp.sum.PUMU.plot ,length)
lipid.exp.sum.PUMU.plot$length_chr = factor(lipid.exp.sum.PUMU.plot$length, levels = unique(lipid.exp.sum.PUMU.plot$length))
lipid.exp.sum.PUMU.plot$unsaturation_PUMU = factor(lipid.exp.sum.PUMU.plot$unsaturation_PUMU, levels = c("Polyunsaturated", "Monounsaturated", "Saturated"))
lipid.exp.sum.PUMU.plot$Class = factor(lipid.exp.sum.PUMU.plot$Class, levels = names(class_color))


for (class_lipid in names(class_color)) {
  
  lipid.exp.sum.PUMU.plot_sub = lipid.exp.sum.PUMU.plot[lipid.exp.sum.PUMU.plot$Class == class_lipid,]
  
  for (i in AgeG) {
    if (i == 'All') {
      lipid.exp.sum.PUMU.plot_sub[,paste0(i,'_sum')] = apply(lipid.exp.sum.PUMU.plot_sub[,grep('^[A-G]\\d', colnames(lipid.exp.sum.PUMU.plot_sub))], 1, sum)
      
    }else {
      lipid.exp.sum.PUMU.plot_sub[,paste0(i,'_sum')] = apply(lipid.exp.sum.PUMU.plot_sub[,grep(paste0("^",i,'\\d'), colnames(lipid.exp.sum.PUMU.plot_sub))], 1, sum)
      
    }
  }
  
  lipid.class.sum.wide.plot = lipid.exp.sum.PUMU.plot_sub[,c(which(colnames(lipid.exp.sum.PUMU.plot_sub)%in%c('Group.1',"Type_count","unsaturation_PUMU","length_chr")), 
                                                             grep("^\\w_sum", colnames(lipid.exp.sum.PUMU.plot_sub)) )]
  
  lipid.class.sum.long.plot = gather(lipid.class.sum.wide.plot, Group_sum,expression, -c('Group.1',"Type_count","unsaturation_PUMU","length_chr"))
  
  lipid.class.sum.long.plot$expression_log = log10(lipid.class.sum.long.plot$expression)
  
  ggplot(lipid.class.sum.long.plot,aes(unsaturation_PUMU, length_chr))+
    geom_tile(aes(fill=expression),color="white",linewidth=1)+ 
    scale_x_discrete(expand = expansion(add = 0))+ 
    scale_y_discrete(expand = expansion(mult = 0))+
    scale_fill_gradientn(colours = colorRampPalette(rev(sequential_hcl(5, palette = "Reds")))(50)) +
    coord_fixed(ratio = 0.3) +
    facet_wrap(vars(Group_sum),ncol = 7)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(size=8,angle = 90,hjust = 1,vjust =0.5),
      axis.text.y.left = element_text(size = 8), 
      legend.title = element_blank() 
    )
  
  ggsave(file = paste0(tem_dir,"/Combine_",class_lipid,"_length_unsaturation_Expsize_Heatmap.pdf"),width = 10,height = 10)
  
}


