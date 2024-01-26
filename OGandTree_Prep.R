library("tidyverse")
library("ape")


#Creating Directory for cleaned OGs
dir.create("./HyPhy_inputs/OG_alignments_cleaned")


############# Removing Stop Codons from Orthogroups #####################

#Getting list of orthogroups
og_files <- list.files("./HyPhy_inputs/OG_alignments/",pattern = "_aligned.fa")

#Making function
Removing_stops_from_OGs <- function(cdsFile) {
  system(paste("(echo 1; echo ./HyPhy_inputs/OG_alignments/",cdsFile,"; echo 1; echo ./HyPhy_inputs/OG_alignments_cleaned/",cdsFile,") | /usr/local/bin/hyphy ./hyphy/res/TemplateBatchFiles/CleanStopCodons.bf", sep = ""))
}

Possibly_Removing_stops_from_OGs <- possibly(Removing_stops_from_OGs, otherwise = "error", quiet = TRUE)

#Running function
map(og_files, Possibly_Removing_stops_from_OGs)


################ Preparing Trees for Analysis ###########################

#Getting list of trees
tree_files <- list.files("./HyPhy_inputs/Labeled_Trees/")

#Making function
Preparing_Trees <- function(tree_file){
  
  #Removing vertical bars(|) from branch tip names
  system(paste("sed -i'.original' -e 's|\\||_|g' ./HyPhy_inputs/Labeled_Trees/", tree_file, sep = ""))
  
  #Adding semicolon to end of trees (so they'll read in)
  treeText <- readr::read_file(paste("./HyPhy_inputs/Labeled_Trees/", tree_file, sep = ""))
  treeText <- gsub('$', ';', treeText)
  tree <- read.tree(text = treeText)
  
  #Removing periods from tip names in trees
  tree[["tip.label"]] <- gsub('\\.', '_', tree[["tip.label"]])
  write.tree(tree, file = paste("./HyPhy_inputs/Labeled_Trees/", tree_file, sep = ""))
  
}

Possibly_Preparing_Trees <- possibly(Preparing_Trees, otherwise = "error", quiet = TRUE)

#Running function
map(tree_files, Possibly_Preparing_Trees)
