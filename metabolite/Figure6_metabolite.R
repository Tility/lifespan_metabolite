setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/')
library(ggplot2)#‘3.5.1’
library(tidyverse)
library(stringr)
library(grid)
library(circlize)
library(patchwork)#‘1.3.0’
library(openxlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorspace)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure6/'

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig6 Data prepare =====
load(file = paste0(dataDir,'/Metabolite_anno_ForKEGG.Rdata'),verbose = T)
class_color = unique(class$class_color)
names(class_color) = levels(class$Class)

load(file = paste0(dataDir,'/Metabolite_expFor_DEM.Rdata'),verbose = T)

data = log10(life.metab.normal)

pt_clinical = read.csv('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/patient_clinical_FG_2.csv',header=TRUE,row.names = 1)
pt_clinical1 = as.data.frame(t(pt_clinical))

data = data[colnames(pt_clinical), ]
dim(data)
##[1]   33 1931
data[1:4,1:4] 

## Marker Metabolite of F or G
load(file = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure3/Marker_signif_FC.Rdata', verbose = T)
marker.fcdf.up = marker.fcdf[marker.fcdf$change=='Up',]
marker.fcdf.up = arrange(marker.fcdf.up,Group,t.test_BH,desc(LOG2FC))


### F G Class with clinical corralation ======
data_f2 <- data[1:11, marker.fcdf.up$metabolite[marker.fcdf.up$Group_change == 'MarkerOfF_Up']]
marker_g = marker.fcdf.up[marker.fcdf.up$Group_change == 'MarkerOfG_Up',]
marker_gTop50 = marker_g[1:50,]
data_g2 <- data[12:33, marker_gTop50$metabolite]

head(data_f2);head(data_g2)


Class_clc_f.r <- round(cor(data_f2, pt_clinical1[1:11,],method = "spearman"), 3)

result = cor.test.mod(x = data_f2,y = pt_clinical1[1:11,], method = "spearman")
head(result$p)
min(result$p)


Class_clc_f.p <- result$p
cat(paste0('Min p value:',min(Class_clc_f.p)))
print('-----')
min(Class_clc_f.r)

min(Class_clc_f.p)
save(Class_clc_f.r,Class_clc_f.p, file = paste0(outDir,'/Marker_clinical_spearman_f.Rdata'))


Class_clc_g.r <- round(cor(data_g2, pt_clinical1[12:33,],method = "spearman"), 3)
result = cor.test.mod(x = data_g2,y = pt_clinical1[12:33,], method = "spearman")
head(result$p)
min(result$p)


Class_clc_g.p<- result$p
cat(paste0('Min p value:',min(Class_clc_g.p)))
min(Class_clc_g.r)

min(Class_clc_g.p)
save(Class_clc_g.r,Class_clc_g.p, file = paste0(outDir,'/Marker_clinical_spearman_g.Rdata'))

### Heat map ======
write.xlsx(Class_clc_f.r, file = paste0(outDir, '/F_class_Marker.xlsx'))
write.xlsx(Class_clc_f.p, file = paste0(outDir, '/F_class_Marker_pvalue.xlsx'))

write.xlsx(Class_clc_g.r, file = paste0(outDir, '/G_class_Marker.xlsx'))
write.xlsx(Class_clc_g.p, file = paste0(outDir, '/G_class_Marker_pvalue.xlsx'))

load(file = paste0(outDir,'/Marker_clinical_spearman_f.Rdata'))
load(file = paste0(outDir,'/Marker_clinical_spearman_g.Rdata'))

#color
col_fun = colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#3B66A2","#B9BFCC",'white',"#D9B4AE","#9E3D3D"))

Class_clc_g.r.sig = Class_clc_g.r
Class_clc_g.r.sig[Class_clc_g.p>0.05]=0
write.xlsx(Class_clc_g.r.sig, file = paste0(outDir, '/G_class_Marker_r_signifi.xlsx'))

pdf(paste0(outDir,'Marker_clinical_PCC_metab_anno_g.pdf'),width = unit(10, "cm"), height = unit(10, "cm"))
Heatmap(Class_clc_g.r.sig,
        col = col_fun,
        border = 'grey30',
        rect_gp = gpar(col = "grey30", lwd = 1),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 90,
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = 12),
          title = " "
        ),
        width = 200, height = 100 
)
dev.off()


Class_clc_f.r.sig = Class_clc_f.r
Class_clc_f.r.sig[Class_clc_f.p>0.05]=0
write.xlsx(Class_clc_f.r.sig, file = paste0(outDir, '/F_class_Marker_r_signifi.xlsx'))
pdf(paste0(outDir,'Marker_clinical_PCC_metab_anno_f.pdf'),width = unit(10, "cm"), height = unit(6, "cm"))
Heatmap(Class_clc_f.r.sig,
        col = col_fun,
        border = 'grey30',
        rect_gp = gpar(col = "grey30", lwd = 1),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 90,
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = 12),
          title = " "
        ),
        width = 200, height = 100 
)
dev.off()


