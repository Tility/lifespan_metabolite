setwd('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite')
library(glmnet)
library(Metrics) ## MAE mean absolute error
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(stringr)
library(openxlsx)
library(RColorBrewer)
source('/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/function_use.R')

dataDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure2/'
outDir = '/Users/liangtingting/hp/lab_file/lifespan_metabolism/metabolite/1_metabolite_figure/Figure7/'

if(!dir.exists(outDir)){
  dir.create(outDir)
}
## Fig7 Data prepare =====
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

## splite sample for trainning dataset ====

set.seed(20240)
Index_p <- sample(2, nrow(data), replace = TRUE, prob = c(0.8, 0.2)) 

data_train <- data[Index_p==1, ] #the training data set
data_test <- data[Index_p==2, ] #the test data set
save(data_train, data_test, file = paste0(outDir, '/Metab_trainData_seed20240_0.8.Rdata'))

Gender_info  = c(TRUE, FALSE)

for (Gender in Gender_info) {
  
  if (Gender) {
    gender_dir = paste0(outDir, '/Gender_add_result/')
  }else {
    data_train = data_train[, colnames(data_train)%!in%c('Gender_num') ]
    data_test = data_test[, colnames(data_test)%!in%c('Gender_num') ]
    gender_dir = paste0(outDir, '/Gender_rm_result/')
  }
  if (!dir.exists(gender_dir)){
    dir.create(gender_dir)}
  
  MAE.list = list()
  train_sizes = seq(0.1, 0.8, 0.1)
  for (train_size in train_sizes) {
    
    tem_dir = paste0(gender_dir, '/train_size', train_size, '_result/')
    
    if (!dir.exists(tem_dir)){
      dir.create(tem_dir)}
    
    set.seed(202401)
    
    prop = 1.25*train_size
    
    Index_p <- sample(2, nrow(data_train), replace = TRUE, prob = c(prop, 1-prop)) 
    
    data_train_sub <- data_train[Index_p==1, ] #the training data set
    print(dim(data_train_sub))
    data_test_sub <- data_train[Index_p==2, ] #the test data set
    print(dim(data_test_sub))
    
    save(data_train_sub, data_test_sub,file = paste0(tem_dir, '/Metab_trainData_seed202401_train_size', train_size, '.Rdata'))
    
    ## glmnet elastic Net ==== n << d 
    data_train_x = data_train_sub[, colnames(data_train_sub)%!in%c("Age","group", 'Gender') ]
    
    x <- as.matrix(data_train_x) ## only log(expression_norm) + Gender
    
    y <- as.matrix(data_train_sub$Age)
    
    alphas = seq(0.1, 0.9, 0.1)
    
    for (alpha in alphas) {
      cv.fit <- cv.glmnet(x,y,alpha =alpha,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 10)
      plot(cv.fit)
      
      cv.fit$lambda.min
      
      predc.trainDF = data.frame(real_Age = data_train_sub[, 'Age']
      )
      predc.testDF = data.frame(real_Age = data_test[, 'Age']
      )
      
      predCV.train <- predict(cv.fit, newx = x[,],
                              s = "lambda.min"
      )
      predc.trainDF$predic_Age = predCV.train
      predCV.test <- predict(cv.fit, newx = as.matrix(data_test[, colnames(data_test)%!in%c("Age","group", 'Gender') ]),
                             s = "lambda.min"
      )
      predc.testDF$predic_Age = predCV.test
      MAE.list[[paste0('MAE.test_', 'trainsize',train_size,'_alpha', alpha)]] = mae(data_test[, 'Age'], predCV.test)
      MAE.list[[paste0('MAE.train_', 'trainsize',train_size,'_alpha', alpha)]] = mae(data_train_sub[, 'Age'], predCV.train)
      
      ggplot(data = predc.trainDF,aes(x = real_Age,y = predic_Age))+
        geom_point(size = 3,shape = 19,colour = "black")+
        labs(x = "Chronological age (years)",y = "Predicted age (years)", title = 'Train')+
        stat_smooth(method = lm,se = F,colour = "dodgerblue3")+
        stat_cor(method = "pearson",label.x.npc="right", label.y.npc="top",color = '#C65248',hjust=1)+
        stat_cor(method = "spearman",color = '#3681B7')+
        theme_bw(  
        )+
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 14, hjust = 0.5),
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 14),
              panel.grid = element_blank())
      
      ggsave(paste0(tem_dir, '/Plot_alpha',alpha,'train_dot.pdf'), width = 6, height = 6)
      
      ggplot(data = predc.testDF,aes(x = real_Age,y = predic_Age))+
        geom_point(size = 3,shape = 19,colour = "black")+
        labs(x = "Chronological age (years)",y = "Predicted age (years)", title = 'Test')+
        stat_smooth(method = lm,se = F,colour = "dodgerblue3")+
        stat_cor(method = "pearson",label.x.npc="right", label.y.npc="top",color = '#C65248',hjust=1)+
        stat_cor(method = "spearman",color = '#3681B7')+
        theme_bw(  
        )+
        theme(aspect.ratio = 1,
              plot.title = element_text(size = 14, hjust = 0.5),
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 14),
              panel.grid = element_blank())
      ggsave(paste0(tem_dir, '/Plot_alpha',alpha,'test_dot.pdf'), width = 6, height = 6)
      
      if (train_size!=0.8) {
        predc.testsubDF = data.frame(real_Age = data_test_sub[, 'Age']
        )
        predCV.testsub <- predict(cv.fit, newx = as.matrix(data_test_sub[, colnames(data_test_sub)%!in%c("Age","group", 'Gender') ]),
                                  s = "lambda.min"
        )
        predc.testsubDF$predic_Age = predCV.testsub
        MAE.list[[paste0('MAE.testsub_', 'trainsize',train_size,'_alpha', alpha)]] = mae(data_test_sub[, 'Age'], predCV.testsub)
        
        ggplot(data = predc.testsubDF,aes(x = real_Age,y = predic_Age))+
          geom_point(size = 3,shape = 19,colour = "black")+
          labs(x = "Chronological age (years)",y = "Predicted age (years)", title = 'Test')+
          stat_smooth(method = lm,se = F,colour = "dodgerblue3")+
          stat_cor(method = "pearson",label.x.npc="right", label.y.npc="top",color = '#C65248',hjust=1)+
          stat_cor(method = "spearman",color = '#3681B7')+
          theme_bw(  
          )+
          theme(aspect.ratio = 1,
                plot.title = element_text(size = 14, hjust = 0.5),
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 14),
                panel.grid = element_blank())
        ggsave(paste0(tem_dir, '/Plot_alpha',alpha,'testsub_dot.pdf'), width = 6, height = 6)
        
        save(predCV.testsub, predc.testsubDF, file = paste0(tem_dir, '/Glmnet_testsub_','trainsize',train_size,'_alpha',alpha,'.Rdata'))
      }
      
      save(cv.fit, predCV.train, predCV.test, predc.trainDF, predc.testDF, file = paste0(tem_dir, '/Glmnet_','trainsize',train_size,'_alpha',alpha,'.Rdata'))
      
    }
    
  }
  
  save(MAE.list, file = paste0(gender_dir, '/Glmnet_MAE_all.Rdata'))
  
  load(file = paste0(gender_dir, '/Glmnet_MAE_all.Rdata'))
  
  MAE.df = data.frame(MAE.list)
  
  MAE.df = as.data.frame(t(MAE.df))
  colnames(MAE.df) = 'MAE'
  MAE.df[,c('anno' ,'trainsize','alpha_need')] = str_split(rownames(MAE.df), '_', simplify = TRUE)[,1:3]
  MAE.df$Type = str_split(MAE.df$anno, '\\.', simplify = TRUE)[,2]
  
  colman = c('#608B99','#807FA2','#B2D28D','#B19BBC','#7F5067',
             '#54629B','#E19B5E','#D06E5D','#B4DCD5','#F7DDB9',
             '#D97484','#E88A5E','#BCBCBA','#B36F9F'
  )
  color_v = colorRampPalette(colman)(length(unique(MAE.df$alpha_need)))
  
  ggplot(data = MAE.df[MAE.df$Type == 'test',],aes(x = trainsize,y = MAE, color = alpha_need,group = alpha_need))+
    geom_line(
      alpha = 1, linewidth = 1) +
    geom_point(size = 3,shape = 19)+
    scale_color_manual(values = color_v) +
    labs(x = "transet size",y = "MAE", title = 'Test')+
    theme_bw()+
    theme(aspect.ratio = 0.5,
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 14),
          panel.grid = element_blank())
  
  ggsave(paste0(gender_dir, '/Plot_MAEtest_dot_test.pdf'), width = 8, height = 6)
  
  ggplot(data = MAE.df[MAE.df$Type == 'train',],aes(x = trainsize,y = MAE, color = alpha_need,group = alpha_need))+
    geom_line(
      alpha = 1, linewidth = 1) +
    geom_point(size = 3,shape = 19)+
    scale_color_manual(values = color_v) +
    labs(x = "transet size",y = "MAE", title = 'Train')+
    theme_bw()+
    theme(aspect.ratio = 0.5,
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 14),
          panel.grid = element_blank())
  
  ggsave(paste0(gender_dir, '/Plot_MAEtest_dot_train.pdf'), width = 8, height = 6)
  
  ggplot(data = MAE.df[MAE.df$Type == 'testsub',],aes(x = trainsize,y = MAE, color = alpha_need,group = alpha_need))+
    geom_line(
      alpha = 1, linewidth = 1) +
    geom_point(size = 3,shape = 19)+
    scale_color_manual(values = color_v) +
    labs(x = "transet size",y = "MAE", title = 'Test sub')+
    theme_bw()+
    theme(aspect.ratio = 0.5,
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 14),
          panel.grid = element_blank())
  
  ggsave(paste0(gender_dir, '/Plot_MAEtest_dot_testsub.pdf'), width = 8, height = 6)
  
}

## min MAE result ====

train_size = 0.8
alpha = 0.1

gender_dir = paste0(outDir, '/Gender_add_result/')
load(file = paste0(gender_dir, '/Glmnet_MAE_all.Rdata'))
MAE_test = MAE.df[MAE.df$Type == 'test',]
top_n(MAE_test, n = -1, MAE)
openxlsx::write.xlsx(MAE_test, file = paste0(gender_dir, '/MAE_multiMedel.xlsx'),rowNames = T)


tem_dir = paste0(gender_dir, '/train_size', train_size, '_result/')
load(file = paste0(tem_dir, '/Glmnet_','trainsize',train_size,'_alpha',alpha,'.Rdata'), verbose = TRUE)
load(file = paste0(tem_dir, '/Metab_trainData_seed202401_train_size', train_size, '.Rdata'), verbose = TRUE)

data_test_out = as.data.frame(t(data_test))
data_test_out$metabolite = metab.normal.log$Metabolite[match(rownames(data_test_out),metab.normal.log$Metab_ID)]

data_train_sub_out = as.data.frame(t(data_train_sub))
data_train_sub_out$metabolite = metab.normal.log$Metabolite[match(rownames(data_train_sub_out),metab.normal.log$Metab_ID)]

openxlsx::write.xlsx(data_test_out, file = paste0(tem_dir, '/sample_test_N33.xlsx'),rowNames = TRUE)
openxlsx::write.xlsx(data_train_sub_out, file = paste0(tem_dir, '/sample_train_N103.xlsx'),rowNames = TRUE)

predc.testDF$sample = rownames(data_test)

openxlsx::write.xlsx(predc.testDF, file = paste0(tem_dir, '/predc.testAge.xlsx'),rowNames = TRUE)
openxlsx::write.xlsx(predc.trainDF, file = paste0(tem_dir, '/predc.trainAge.xlsx'),rowNames = TRUE)

coef_df = as.data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
colnames(coef_df) = 'coefficients'
coef_df$Metabolite = metab.normal.log$Metabolite[match(rownames(coef_df),metab.normal.log$Metab_ID)]

coefm_df = coef_df[!is.na(coef_df$Metabolite),]
colnames(class)[1] = 'Metabolite'
coefm_df = merge(coefm_df, class,by = 'Metabolite')
coefm_df = coefm_df %>%
  arrange(desc(coefficients)) %>%
  mutate("Rank" = row_number())

write.xlsx(coefm_df, file = paste0(tem_dir, '/Metab_coeff_all.xlsx'))

save(coef_df, coefm_df, class, class_color,file = paste0(tem_dir, '/Feature_coefficients.Rdata'))

ggplot(data = coefm_df,aes(x = Rank,y = coefficients, color = Class))+
  geom_point(size = 3,shape = 19)+
  scale_color_manual(values = class_color) +
  geom_text_repel(
    data = coefm_df[c(1:2,1928:1931),],
    aes(x = Rank,y = coefficients,label = Metabolite),
    size = 4,
    min.segment.length = 0,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  ylim(-21,10) +
  labs(x = "Rank",y = "coefficients", title = 'Test')+
  theme_bw()+
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14),
        panel.grid = element_blank())

ggsave(paste0(tem_dir, '/Plot_coeff_Metab.pdf'), width = 8, height = 6)


## For Combine metabolite and lipid ====
## We used log-transformed data; the above steps were repeated in combination with metabolite and lipid data.

