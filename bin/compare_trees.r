library(ape)


newick_files<-list.files('.',full.names = TRUE, pattern = "*.newick")


conditions = c(c('human','chimp'), c('human','marmoset'))
conditions = factor('human', 'chimp', 'marmoset')

col1 = c()
col2 = c()
for (i in newick_files){ 
  for(j in newick_files){
      if (i != j){ 
        col1 = append(col1, i)
        col2 = append(col2,j)
      }
    }
}

combo_df<-data.frame(col1,col2)

comp1<-c()
comp2<-c()
oneNotInTwo<-c()
twoNotInOne<-c()
nodeDiff<-c()
tipDiff<-c()


for(i in 1:nrow(combo_df)){ 
  
  print(i)
  
  animal1 <- combo_df$col1[i]
  animal2 <- combo_df$col2[i]
  
  name1<-strsplit(combo_df$col1[i],'[.]')[[1]][1]
  assign(combo_df$col1[i], animal1)
  
  name2<-paste0(combo_df$col2[i]) 
  assign( combo_df$col2[i], animal2)

  x<-comparePhylo(read.tree(animal1), read.tree(animal2))
  
  
  comp1 = append(comp1,paste0(strsplit(basename(combo_df$col1[i]), '[.]')[[1]][1]))
  comp2 = append(comp2, paste0(strsplit(basename(combo_df$col2[i]), '[.]')[[1]][1]))
  
  oneNotInTwo = append(oneNotInTwo, x$messages[3])
  twoNotInOne = append(twoNotInOne, x$messages[4])
  nodeDiff = append(nodeDiff, x$messages[5])
  tipDiff = append(tipDiff, x$messages[2])
  
  
  }


info_df<-data.frame(comp1, comp2, oneNotInTwo, twoNotInOne, nodeDiff, tipDiff)

 
write.csv(info_df,"comparison_info.csv")


