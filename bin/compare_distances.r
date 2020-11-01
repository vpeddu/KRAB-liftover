library(ape)
library(ggplot2)
library(ggtree)

setwd("/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/testing")

trees<-list.files('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree',full.names = TRUE)

human_tree<-trees[which(grepl('human',trees))]

non_human_trees<-trees[-which(grepl('human',trees))]

human_tree<-read.tree(human_tree)

human_cophenetic<-cophenetic(human_tree)

#comparisons<-c()
for (i in non_human_trees){ 
  base = strsplit(basename(i),'[.]')[[1]][1]
  
  comparison_tree<-read.tree(i)
  comparison_cophenetic<-cophenetic(comparison_tree)
  in_comparison<-colnames(comparison_cophenetic)
  
  
  gene<-c()
  comparison_val<-c()
  human_val<-c()
  bad_list<-c()
  for(i in in_comparison){ 
    
    comparison_col_index<-which(colnames(comparison_cophenetic) == i)
    
    nearest_neighbor = sort(comparison_cophenetic[comparison_col_index,])[2]
    #print(i)
    #print(nearest_neighbor)
    human_col<-which(colnames(human_cophenetic) == i)
    #print(human_col)
    human_row<-which(rownames(human_cophenetic) == names(nearest_neighbor))
    
    human_index<-human_cophenetic[human_row,human_col]
    #print(names(nearest_neighbor))
    #print(human_index)
    #print('')
    
    if(identical(human_index, numeric(0))){ 
      bad_list<-append(bad_list, names(nearest_neighbor))
    }
    else{
      gene<-append(gene,names(nearest_neighbor))
      comparison_val<-append(comparison_val,nearest_neighbor)
      human_val<-append(human_val,human_index)
    }
  }
  
  human_comparison_df<-data.frame(gene,comparison_val,human_val)
  human_comparison_df$diff<-human_comparison_df$human_val - human_comparison_df$comparison_val
  
  write.csv(human_comparison_df,paste0(base,'_cophenetic_comparison.csv'))
  
  plot<-ggplot(human_comparison_df, aes( x = gene, y = diff)) + 
    geom_point() + 
    geom_text(aes(label=ifelse(diff>2,as.character(gene),'')),hjust=0,vjust=0) + 
    ylab(paste0(base,' cophenetic - Human cophenetic ')) + 
    xlab('gene') + 
    geom_hline(yintercept = .5) + 
    theme_classic() + 
    theme(axis.text.x = element_blank())
  plot
  
  
  }




