library(ape)
library(ggplot2)
library(ggtree)
library(seqinr)
library(ips)
library(reshape2)
library(adephylo)

ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(seq_name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

# Usage example
df<-ReadFasta(args[1])
#df$header<-as.character(df$header)
#df$seq_name<-as.character(df$seq_name)


setwd("/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/testing/combined_trees/")



combined_tree<-read.tree('/Users/vikas/Documents/UCSC/rotations/Haussler/nf_workflow/KRAB-liftover/uniprot_output/combined_output/tree/combined.tree.newick')

distance_matrix<-as.matrix(distTips(combined_tree, tips = all, method = 'nNodes' ))

#get unique genes
unique_genes<-unique(sapply(combined_tree$tip.label, function(x) strsplit(x, "_")[[1]][2], USE.NAMES=FALSE))

not_found_list<-c()
bad_list<-c()

gene<-c()
human_dist<-c()
chimp_dist<-c()
marmoset_dist<-c()
for(i in unique_genes){ 
  marmoset_index<-which(colnames(distance_matrix) == paste0('marmoset_',i))
  
  if(identical(integer(0), marmoset_index)){ 
    not_found_list<-append(not_found_list,i)
    }
  else{ 
    human_index<-which(rownames(distance_matrix) == paste0('human_',i))
    chimp_index<-which(rownames(distance_matrix) == paste0('chimp_',i))
    
    if(identical(chimp_index, integer(0)) || identical(human_index, integer(0))){ 
        bad_list<-append(bad_list, i)
      }
    else{ 
      print(chimp_index)
      print(i)
      gene<-append(gene,i)
      human_dist<-append(human_dist, distance_matrix[human_index,marmoset_index])
      chimp_dist<-append(chimp_dist, distance_matrix[chimp_index,marmoset_index])
      }
    }
  }


distance_compared_to_marmoset<-data.frame(gene, human_dist, chimp_dist)
# distance_compared_to_marmoset$human_dist<-distance_compared_to_marmoset$human_dist - 1 
# distance_compared_to_marmoset$chimp_dist<- distance_compared_to_marmoset$chimp_dist - 1


melted_distance_df<-melt(distance_compared_to_marmoset)

melted_distance_df$variable<-as.factor(melted_distance_df$variable)

distance_plot<-ggplot(melted_distance_df,aes( x = gene , y = value, color = variable, shape = variable)) +
  geom_point(height =0, width = 1, size = 3) + 
  geom_text(aes(label=ifelse(value>10,as.character(gene),'')),hjust=1,vjust=0)+ 
  theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  ylab('Topological distance to equivalent Marmoset gene')
distance_plot

## MSA 

