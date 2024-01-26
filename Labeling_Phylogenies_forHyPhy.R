library("tidyverse")
library("googlesheets4")
library("ape")

#Reading in google sheet
Aggression_Phenotypes <- read_sheet("https://docs.google.com/spreadsheets/d/1L1h-0y_NwunlBJbj9SgTN2lqtvMBzEFoLeLBiYWZ0pM/edit?usp=sharing")


######### Getting Species Codes ####################

#Splitting genus and species
Aggression_Phenotypes$Genus <- str_split_i(Aggression_Phenotypes$Species, " ", 1)
Aggression_Phenotypes$Species <- str_split_i(Aggression_Phenotypes$Species, " ", 2)

#Isolating parts of each needed for species code
Aggression_Phenotypes$Genus <- str_sub(Aggression_Phenotypes$Genus, start = 1, end = 1)
Aggression_Phenotypes$Species <- str_sub(Aggression_Phenotypes$Species, start = 1, end = 3)

#Putting together to form species code and making lowercase
Aggression_Phenotypes$Species_code <- paste(Aggression_Phenotypes$Genus, Aggression_Phenotypes$Species, sep = "") %>% tolower()



######## Getting High Discrimination and Aggression Species ##################

High_Disc <- filter(Aggression_Phenotypes, `Discrimination phenotype`=="high")

High_Aggr <- filter(Aggression_Phenotypes, grepl("high", `Aggression phenotype`))



######### Making and applying labels ##########################

#Creating directory for labeled trees
dir.create("./HyPhy_inputs/Labeled_Trees")


#Making list of input tree files
input_trees <- list.files("./HyPhy_inputs/Tree_Inferences", full.names = TRUE)

#Making into function

Labeling_Phylogenies <- function(input_tree) {
  
  #Reading in tree
  tree <- read.tree(file = input_tree)
  
  #Making labels
  Tip_labels <- tree$tip.label
  Tl_Disc_ants <- Tip_labels[str_detect(Tip_labels, str_c(High_Disc$Species, collapse="|"))]
  Tl_Aggr_ants <- Tip_labels[str_detect(Tip_labels, str_c(High_Aggr$Species, collapse="|"))]
  
  ortho_number <- str_split_i(input_tree, "/", 4) %>% str_split_i("\\.", 1)
  
  #Saving to computer as text files w/ specific names
  write(Tl_Disc_ants, paste("./HyPhy_inputs/",
                            ortho_number,
                            "_Disc.txt",
                            sep = ""))
  write(Tl_Aggr_ants, paste("./HyPhy_inputs/",
                            ortho_number,
                            "_Aggr.txt",
                            sep = ""))
  
  system(paste("/usr/local/bin/hyphy ./hyphy-analyses/LabelTrees/label-tree.bf --tree HyPhy_inputs/Tree_Inferences/",ortho_number,".tree  --list HyPhy_inputs/",ortho_number,"_Disc.txt --output HyPhy_inputs/Labeled_Trees/", ortho_number, "_Disc.tree", sep = ""))
  system(paste("/usr/local/bin/hyphy ./hyphy-analyses/LabelTrees/label-tree.bf --tree HyPhy_inputs/Tree_Inferences/",ortho_number,".tree  --list HyPhy_inputs/",ortho_number,"_Aggr.txt --output HyPhy_inputs/Labeled_Trees/", ortho_number, "_Aggr.tree", sep = ""))
}


Possibly_Labeling_Phylogenies <- possibly(Labeling_Phylogenies, otherwise = "error", quiet = TRUE)

#Running Script
map(input_trees, Possibly_Labeling_Phylogenies)


