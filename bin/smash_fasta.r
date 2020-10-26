# smash fasta
args = commandArgs(trailingOnly=TRUE)



# Function
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


df$letter<-sapply(strsplit(as.character(df$seq_name), "_"), `[`, 1)
df$header<-sapply(strsplit(as.character(df$seq_name), "_"), `[`, 2)
df$gene<-df$gene<-sapply(strsplit(as.character(df$header), "::"), `[`, 1)

smashed_exons<-c()
new_headers<-c()
seen<-c()

for( i in df$gene){ 
  matches<-which(df$gene == i)
  if(i %in% seen){ 
    next
  }
  if(df$letter[matches[1]] == "A"){
    smashed<-(paste0(df$sequence[matches][1], (df$sequence[matches][2])))
    smashed_exons<-append(smashed_exons,smashed )
    new_headers<-append(new_headers,i)
    seen<-append(seen,i)
  }
}
smashed_fasta_df<-data.frame(smashed_exons, new_headers)
colnames(smashed_fasta_df)<-c('seq','name')
#DNAString(smashed_fasta_df, names = smashed_fasta_df$new_headers, seq = smashed_fasta_df$smashed_exons)

# writeXStringSet(smashed_fasta_df, 'test.fasta', append=FALSE,
#                 compress=FALSE, compression_level=NA, format="fasta")

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# exampleData = dplyr::data_frame(name = smashed_fasta_df$new_headers,
#                                 seq = smashed_fasta_df$smashed_exons)
writeFasta(smashed_fasta_df, "smashed.fasta")
