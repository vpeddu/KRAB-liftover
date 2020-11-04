library(ape)
library(ggplot2)
library(ggtree)

setwd("/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/testing")



human_tree<-read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree/human.tree.newick')

human_cophenetic<-cophenetic(human_tree)

chimp_tree<-read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree/chimp.tree.newick')

chimp_cophenetic<-cophenetic(chimp_tree)

in_chimp<-colnames(chimp_cophenetic)



gene<-c()
chimp_val<-c()
human_val<-c()
bad_list<-c()
for(i in in_chimp){ 
  
  chimp_col_index<-which(colnames(chimp_cophenetic) == i)
  
  nearest_neighbor = sort(chimp_cophenetic[chimp_col_index,])[2]
  #print(i)
  #print(nearest_neighbor)
  human_col<-which(colnames(human_cophenetic) == i)
  #print(human_col)
  human_row<-which(rownames(human_cophenetic) == names(nearest_neighbor))
  
  human_index<-human_cophenetic[human_row,human_col]
  print(names(nearest_neighbor))
  print(human_index)
  print('')
  
  if(identical(human_index, numeric(0))){ 
    bad_list<-append(bad_list, names(nearest_neighbor))
  }
  else{
    gene<-append(gene,names(nearest_neighbor))
    chimp_val<-append(chimp_val,nearest_neighbor)
    human_val<-append(human_val,human_index)
  }
}

human_chimp_df<-data.frame(gene,chimp_val,human_val)
human_chimp_df$diff<-human_chimp_df$human_val - human_chimp_df$chimp_val

interesting_genes<-human_chimp_df$gene[human_chimp_df$diff > 2]
for(i in interesting_genes){ 
  print(i)
  chimp_col_index<-which(colnames(chimp_cophenetic) == i)
  #print(chimp_col_index)
  nearest_neighbor = sort(chimp_cophenetic[chimp_col_index,])[2]
  #print(which(chimp_cophenetic == chimp_cophenetic[chimp_col_index]))
  print(nearest_neighbor)
  
}

plot<-ggplot(human_chimp_df, aes( x = gene, y = diff)) + 
  geom_point() + 
  geom_text(aes(label=ifelse(diff>2,as.character(gene),'')),hjust=1,vjust=0) + 
  ylab('Chimp cophenetic - Human cophenetic ') + 
  xlab('Gene') + 
  geom_hline(yintercept = 2, linetype = 'dashed') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())
plot
ggsave(plot = plot, 'patristic_diff_chimp-human.pdf', height = 10, width = 10)

especially_interesting <- human_chimp_df$gene[human_chimp_df$diff > 2]
especially_interesting<-append(c('ZNF473', 'ZNF28', 'ZNF496', 'ZNF597'), especially_interesting)
fill_col<-c(rep('Tycko',4), rep('Vikas',length((human_chimp_df$gene[human_chimp_df$diff > 2]))))

metadata_df<-data.frame(especially_interesting,fill_col)


#human tree
human_473_tree<-ggtree::read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree/human.tree.newick')
for (i in 1:length(human_473_tree$tip.label)){
  #print(human_473_tree$tip.label[i])
  if (human_473_tree$tip.label[i] %in% especially_interesting){ 
      print('ok')
    }
  else{ 
    human_473_tree$tip.label[i]<-NA
    }
  }
human_473_plot<-ggtree(human_473_tree) 

#plot(human_473_plot)

human_473_plot<-human_473_plot %<+% metadata_df + 
  geom_tiplab(aes(fill = factor(fill_col)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  theme(legend.position = c(0.5,0.2), 
        legend.title = element_blank()) # no title
        #legend.key = element_blank()) # no keys
ggsave(plot = human_473_plot, file = 'human_473_plot.pdf', height = 10, width = 10)


#chimp
chimp_473_tree<-ggtree::read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree/chimp.tree.newick')
for (i in 1:length(chimp_473_tree$tip.label)){
  #print(chimp_473_tree$tip.label[i])
  if (chimp_473_tree$tip.label[i] %in% especially_interesting){ 
    print('ok')
  }
  else{ 
    chimp_473_tree$tip.label[i]<-NA
  }
}
chimp_473_plot<-ggtree(chimp_473_tree) 

#plot(chimp_473_plot)

chimp_473_plot<-chimp_473_plot %<+% metadata_df + 
  geom_tiplab(aes(fill = factor(fill_col)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  theme(legend.position = c(0.5,0.2), 
        legend.title = element_blank()) # no title
#legend.key = element_blank()) # no keys

ggsave(plot = chimp_473_plot, file = 'chimp_473_plot.pdf', height = 10, width = 10)


#marmoset
marmoset_473_tree<-ggtree::read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/fastTree/marmoset.tree.newick')
for (i in 1:length(marmoset_473_tree$tip.label)){
  #print(marmoset_473_tree$tip.label[i])
  if (marmoset_473_tree$tip.label[i] %in% especially_interesting){ 
    print('ok')
  }
  else{ 
    marmoset_473_tree$tip.label[i]<-NA
  }
}
marmoset_473_plot<-ggtree(marmoset_473_tree) 

#plot(marmoset_473_plot)

marmoset_473_plot<-marmoset_473_plot %<+% metadata_df + 
  geom_tiplab(aes(fill = factor(fill_col)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  theme(legend.position = c(0.5,0.2), 
        legend.title = element_blank()) # no title
#legend.key = element_blank()) # no keys

ggsave(plot = marmoset_473_plot, file = 'marmoset_473_plot.pdf', height = 10, width = 10)
