# Negation of %in%
'%!in%' <- Negate('%in%')

cv_metab=function(x){  #  sd/mean
  i=na.omit(x)
  return(sd(i)/mean(i))
}

AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

FindMarker_metabolite = function(metab.group, metab.normal, metab.normal.log, 
                                 class = class,
                                 idents.1 = 'A',idents.2 = NULL,  marker = FALSE,
                                 LOG2FC = 0.584963, p_adj = 0.05,
                                 pAdjustMethod = "fdr", outDir = outDir, 
                                 return_data = FALSE){
  if (is.null(idents.2)) {
    marker = TRUE
    metab.group = mutate(metab.group, group = ifelse(group==idents.1, idents.1, 'Others'))
    idents.2 = 'Others'
  }
  
  test1.metab.group <-metab.group[metab.group$group%in%c(idents.1,idents.2),]
  test1.metab.normal.log <-metab.normal.log[rownames(test1.metab.group),]
  
  test1.test_results <- data.frame(Metabolites = colnames(metab.normal)[1:ncol(metab.normal)])
  
  # t.test
  test1.test_results$t.test <- apply(test1.metab.normal.log[,1:ncol(metab.normal)],2,
                                     function(x) unlist(t.test(as.numeric(x) ~ test1.metab.group$group, var.equal=TRUE,#student t
                                                               data = test1.metab.normal.log)[3]))
  # adjust p value using BH method
  test1.test_results$t.test_BH <- p.adjust(test1.test_results$t.test, method = pAdjustMethod)
  
  if (marker) {
    #fold-change use raw peak area
    test1.test_results$FC <- apply(test1.metab.group[,1:ncol(metab.normal)], 2, 
                                   function(x) mean(x[which(test1.metab.group$group == idents.1)])/mean(x[which(test1.metab.group$group == idents.2)]))
    test1.test_results$LOG2FC <- log2(test1.test_results$FC)
    # differential altered metabolites FC 
  }else {
    #fold-change use raw peak area
    test1.test_results$FC <- apply(test1.metab.group[,1:ncol(metab.normal)], 2, 
                                   function(x) mean(x[which(test1.metab.group$group == idents.2)])/mean(x[which(test1.metab.group$group == idents.1)]))
    test1.test_results$LOG2FC <- log2(test1.test_results$FC)
    # differential altered metabolites FC 
  }
  
  test1.test_results$change_t.test <- ifelse(test1.test_results$t.test_BH < p_adj & abs(test1.test_results$LOG2FC) >= LOG2FC, 
                                             ifelse(test1.test_results$LOG2FC > LOG2FC ,'Up','Down'), "Not changed")
  
  cat(table(test1.test_results$change_t.test))
  ##Down Not changed          Up 
  ##777         651         503 
  cat('\n')
  
  # metabolites class and kegg numbers
  test1.test_results$cpd <- class$cpd[match(test1.test_results$Metabolite, class$metabolite)]
  test1.test_results$Class <- class$Class[match(test1.test_results$Metabolite, class$metabolite)]
  
  write.xlsx(test1.test_results,paste0(outDir, "/Metabolomics_",idents.1,"vs",idents.2,"_metabolite_FC.xlsx"))
  
  if (return_data) {
    return(test1.test_results)
  }
  
}


MarkerKEGG = function(df, path_metab = path_metab, path_intro = path_intro,
                      minGSSize = 2, pAdjustMethod = "fdr",
                      pvalueCutoff = 0.99, outDir = outDir, Clust = 'C1',
                      return_data = FALSE
) {
  
  cpd <- class$cpd[match(df, class$metabolite)]
  
  cpd[cpd==""] = 0
  
  if (all(cpd==0)) {
    print("++ No KEGG Number martch ++")
  }else {
    print("KEGG Enriching ...")
    x_enrich <- enricher(cpd, 
                         TERM2GENE = path_metab, TERM2NAME = path_intro,
                         minGSSize = minGSSize, pvalueCutoff = pvalueCutoff, pAdjustMethod = "fdr") 
    kegg_table <- as.data.frame(x_enrich)  
    kegg_table <- na.omit(kegg_table)
    if (dim(kegg_table)[1]!=0) {
      kegg_table$FoldEnrich <- apply(kegg_table, 1,
                                     function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                       as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                       as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                       as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
      kegg_table$Description <- factor(kegg_table$Description, levels = rev(kegg_table$Description))
      
      write.xlsx(kegg_table, file = paste0(outDir, 'MakerM_',Clust,'_KEGG_all.xlsx'))
      
      kegg_table = kegg_table[kegg_table$p.adjust<0.05, ]
      
      if (dim(kegg_table)[1]==0) {
        print("No Sinificant results!!!")
      }
      else {
        write.xlsx(kegg_table, file = paste0(outDir, 'MakerM_',Clust,'_KEGG_signif.xlsx'))
        
        colmax = max(-log10(kegg_table$p.adjust))
        colmin = min(-log10(kegg_table$p.adjust))
        sizemax = max(kegg_table$FoldEnrich)
        sizemin = min(kegg_table$FoldEnrich)
        
        ggplot(kegg_table, aes(x = -log10(p.adjust), y = Description)) + 
          geom_point(aes(fill = -log10(p.adjust),  size = FoldEnrich), shape = 21, color = "grey40") +
          scale_fill_viridis_c(direction = -1, 
                               end = 0.9, 
                               option = "C", 
                               limit = c(as.integer(colmin)-1, as.integer(colmax)+1), 
                               breaks = seq(as.integer(colmin)-1,as.integer(colmax)+1,as.integer((as.integer(colmax)+1)/2))) + 
          scale_size(range = c(1, 8), breaks = seq(as.integer(sizemin), as.integer(sizemax), as.integer((as.integer(sizemax))/4))) + 
          theme_bw() +
          theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 2) +
          labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (p.adjust)") +
          theme(panel.grid.major.x = element_line(size = 0), 
                panel.grid.minor.x = element_line(size = 0)) +
          theme(plot.title = element_text( size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
                axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
                axis.title.y = element_text(size = 6, vjust = 2, face = "plain"),
                axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
                axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
                legend.title = element_text(size = 6, face = "plain", colour = "black"),
                legend.text = element_text(size = 6), 
                legend.key.size = unit(0.4, "cm"), 
                legend.spacing = unit(0.4, "cm"))+
          xlim(c(colmin,colmax))
        
        ggsave(paste0(outDir, '/MakerM_',Clust,'_KEGG_signif.pdf'), width = 7,height = 7)
        
        kegg_table = arrange(kegg_table,pvalue)
        print("KEGG Enrichment Done! ")
      }
      if (return_data) {
        return(kegg_table)
      }
    }
  }
}

line_plot_of_group <- function(data, color, ylim = NULL, line_Tpos = 'D',
                               alpha = 0.05, size = 0.3) {
  # Find average expression over age of all metabolites in cluster
  avg <- data %>%
    group_by(group) %>%
    summarise(across(expression, c(mean = mean, sd = sd, se= ~sd(.)/sqrt(n())) ))
  
  # Define metabolite number annotation df
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = c(paste0("n=", length(unique(data$metab)))),
    hjustvar = c(-0.5) ,
    vjustvar = c(1))
  
  print(max(data$expression))
  print(min(data$expression))
  
  # Find average expression over group of all metabolites in cluster
  #data_avg <- aggregate(apply(data[,!(colnames(data) %in% c('group'))],2,function(x) as.numeric(x)),by=list(data$group), mean, na.rm=TRUE)
  data_avg <- data %>%
    group_by(group, metab) %>%
    summarise(across(expression, c(mean = mean, sd = sd, se= ~sd(.)/sqrt(n())) ))
  
  p1 = ggplot(data_avg, aes(x = group, y = expression_mean)) +
    theme_bw() +
    geom_line(#stat="smooth", method = "loess", span = 0.75, se = FALSE,
      aes(group = metab),
      alpha = alpha, color = color, linewidth = size) +
    scale_y_continuous(expand = c(0.5,0)) +
    geom_line(data = avg,
              aes(x = group, y = expression_mean,group = 1),
              #stat="smooth", method = "loess", span = 0.75, se = FALSE,
              color = darken(color, amount = 0.3), linewidth = size * 3) +
    theme(legend.position = "right",
          aspect.ratio = 1,
          text = element_text(size = 8),
          panel.spacing = unit(c(0, 0, 0, 0), "null")) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
    {if(!is.null(ylim))ylim(ylim)}
  
  if (!is.null(line_Tpos)) {
    # Metabolite  plot
    data_pos = data_avg[data_avg$group==line_Tpos,]
    data_pos$x_pos =  data_pos$group
    data_pos$y_pos =  data_pos$expression_mean
    
    # Metabolite plot
    plot <- p1+geom_text_repel(data =  data_pos,aes(x = x_pos,y = y_pos,label = metab), nudge_x = 1,nudge_y = 0.5,segment.color = 'gray30',
                      color = color,inherit.aes = T,max.overlaps = 20,min.segment.length = unit(0, 'lines'),
                      size = 3) 
  }else {
    plot = p1
  }
  
  return(plot)
}

run_loess <- function(metab, data) {
  # Extract metab specific data
  data <- data[, c(metab, non_metab_cols)]
  
  # Generate LOESS model and save to lo_data
  lo <- loess(get(metab) ~ Age, data, span = 0.75)
  
  # Generate predicted values
  lo_predict <- predict(lo, data.frame(Age = seq(age_min, age_max, 1))) %>%
    as.data.frame()
  colnames(lo_predict) <- metab
  
  return(lo_predict)
}

line_plot_of_loess <- function(data, color, ylim = NULL, 
                               alpha = 0.05, size = 0.3, plot_ind_genes = TRUE) {
  # Find average expression over age of all metabolites in cluster
  avg <- data %>%
    group_by(Age) %>%
    summarise(across(expression, c(mean = mean, sd = sd, se= ~sd(.)/sqrt(n())) ))
  
  # Generate LOESS model and save to lo_data
  lo <- loess(expression_mean ~ Age, avg, span = 0.75)
  lo_predict <- predict(lo, data.frame(Age = avg$Age))
  avg <- add_column(avg, lo = lo_predict)
  
  # Define metabolite number annotation df
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = c(paste0("n=", length(unique(data$metab)))),
    hjustvar = c(-0.5) ,
    vjustvar = c(1))
  
  print(max(data$expression))
  print(min(data$expression))
  
  # Metaboliterate plot
  plot <- ggplot(data, aes(x = Age, y = expression)) +
    theme_bw() +
    geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE,
              aes(group = metab),
              alpha = .05, color = color, linewidth = size) +
    scale_y_continuous(expand = c(0.5,0)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_line(data = avg,
              aes(x = Age, y = expression_mean),
              stat="smooth", method = "loess", span = 0.75, se = FALSE,
              color = darken(color, amount = 0.3), linewidth = size * 3) +
    theme(legend.position = "none",
          aspect.ratio = 1,
          text = element_text(size = 8),
          panel.spacing = unit(c(0, 0, 0, 0), "null")) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
    {if(!is.null(ylim))ylim(ylim)}
  
  return(plot)
}

cor.test.mod = function(x,y,method = "spearman",adjust = 'BH'){
  n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
  r <- cor(x,y,method = method)
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  p.adj = p
  p.adj[] = p.adjust(p,method = adjust)
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, p.adj,se)
  names(out) <- c("r", "n", "t", "p", "p.adj","se")
  return(out)
}
