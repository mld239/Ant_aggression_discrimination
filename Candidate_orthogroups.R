library(Biostrings)
library('rBLAST')
library("purrr")
library("tidyverse")
library(stringr)



################### PART ONE: Finding Candidates in Genomes using BLAST ####################################

#Putting into blast database
makeblastdb(file = "./genomes_cleaned/all_genomes_cleaned.fasta", dbtype = "nucl")


#Making into function
Identify_candidates <- function(candidate_gene) {

  #Reading in candidate gene
  cand_gene <- readDNAStringSet(candidate_gene, format = "fasta")
  
  #Prep the search
  blastSearch <- blast(db = "./genomes_cleaned/all_genomes_cleaned.fasta", type = "blastn")
  
  #Run the search
  searchResults <- predict(blastSearch, cand_gene, BLAST_args = "-max_target_seqs 1") 
  
  #Tagging on gene name
  searchResults$gene_name <- candidate_gene
  
return(searchResults)}

#Making so function won't stop at errors
Possibly_Identify_candidates <- possibly(Identify_candidates, otherwise = "couldn't open file", quiet = TRUE)

#Making list of candidate genes
gene_list <- list.files("./gene_sequences", pattern = "*.fna", full.names = TRUE)

#Running function on list of genes
Blast_results <- map(gene_list, Possibly_Identify_candidates)

#Converting results to data frame
Blast_results_simplified <- as.data.frame(do.call(rbind, Blast_results)) 




################### PART 2: Finding Orthogroups of Candidates ##########################################

#Reading in N0.tsv and consolidating species columns into one column
N0 <- read_tsv("./OrthoFinder/amino_acid_sequences/OrthoFinder/Results_Sep23/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
N0 <- unite(N0, all_species, acep_genomeAssembly.fna.transdecoder:waur_genomeAssembly.fna.transdecoder, sep = ",")

#Replacing pipes (|) in species gene names with underscores (_)
N0$all_species <- gsub("[|]","_",as.character(N0$all_species))

#Making into function
Identify_candidate_orthogroups <- function(candidate_gene) {

  matching_orthogroups <- N0 %>%
    filter(str_detect(all_species, candidate_gene))
  
  matching_orthogroups$candidate_gene <- candidate_gene

  return(matching_orthogroups)
  }

#Making so function won't stop at errors
Possibly_Identify_candidate_orthogroups <- possibly(Identify_candidate_orthogroups, otherwise = "error", quiet = TRUE)

#Making list of candidate genes (and removing the pipes again)---don't know if I actually needed to do this
genes <- Blast_results_simplified$SubjectID
gene_list <- gsub("[|]","_",as.character(genes))
                       
#Running function on list
Candidates_by_orthogroup <- map(gene_list, Possibly_Identify_candidate_orthogroups)
Candidates_by_orthogroup_simplified <- as.data.frame(do.call(rbind, Candidates_by_orthogroup)) %>%
  distinct()

# Fix vertical bars in gene names to underscore:
Blast_results_simplified$SubjectID <- gsub(pattern = "\\|",
                                           replacement = "_",
                                           Blast_results_simplified$SubjectID)
candidate_gene_connectedto_orthogroup <- full_join(Blast_results_simplified, Candidates_by_orthogroup_simplified, by = c("SubjectID" = "candidate_gene"))
candidate_gene_connectedto_orthogroup <- select(candidate_gene_connectedto_orthogroup, c("OG", "gene_name"))
write_csv(candidate_gene_connectedto_orthogroup, file = "./candidates_by_orthogroup.csv")

########################### PART 3a: Getting Tree Inferences and Alignments of Orthogroups and putting them into HyPhy Inputs ###########################

#Making list of OG to search for
OG_list <- Candidates_by_orthogroup_simplified$OG

#Creating Hyphy directory and Tree and Alignment subdirectories
dir.create("HyPhy_inputs")
dir.create("./HyPhy_inputs/Tree_Inferences")
dir.create("./HyPhy_inputs/OG_alignments")

#TREES#

#Making function to move trees to HyPhy directory
Moving_HyPhy_Trees <- function(target_OG) {
  
  file.copy(from = paste("./Tree_Inferences/", target_OG,".tree", sep = ""), to = "./HyPhy_inputs/Tree_Inferences")
  
}

Possibly_Moving_HyPhy_Trees <- possibly(Moving_HyPhy_Trees, otherwise = "error", quiet = TRUE)

#Running Function
map(OG_list, Possibly_Moving_HyPhy_Trees)


#ALIGNMENTS#

#Reading in all genomes file
all_genomes <- read.fasta("genomes_cleaned/all_genomes_cleaned.fasta")
#Removing | and {species_code} things present in later gene names for some reason (DO EVERY TIME)
all_genomes$seq.name.cleaned <- gsub("[|]","_",as.character(all_genomes$seq.name))
all_genomes$seq.name.cleaned <- gsub("_\\{species_code\\}_","_",as.character(all_genomes$seq.name.cleaned))

#Other (possibly faster) way to filter
#filter(all_genomes, grepl())

#Trying to filter all genomes file with this list
#nt_seq <- lapply(gene_list, FUN = function(x) filter(all_genomes, grepl(x, all_genomes$seq.name.cleaned)))
#nt_seq_simp <- as.data.frame(do.call(rbind, nt_seq))

#Making into function
Converting_OG_to_nt <- function(orthogroup) {
  
  #Reading in orthogroup
  OG <- read.fasta(paste("./aligned_HOG/",orthogroup,".fa", sep = ""))
  #Removing ".p" and "|" from gene names and putting in new column
  OG$seq.name.cleaned <- str_split_i(OG$seq.name, "\\.p", 1)
  OG$seq.name.cleaned <- gsub("[|]","_",as.character(OG$seq.name.cleaned))
  
  #Making list of gene names
  gene <- OG$seq.name.cleaned
  
  #Finding matching gene names in all genomes file and making data frame with only the matches
  nt_seq <- lapply(gene, FUN = function(x) filter(all_genomes, grepl(x, all_genomes$seq.name.cleaned)))
  #Simplifying data frame (I might not need this part)
  nt_seq_simp <- as.data.frame(do.call(rbind, nt_seq))
  #Making gene names match names in orthogroup
  nt_seq_simp$seq.name <- OG$seq.name
  
  #Converting dataframe to fasta and saving to specific name
  output_filename <- paste("./HyPhy_inputs/OG_alignments/",orthogroup,".fa",sep = "")
  dat2fasta(dat=nt_seq_simp, outfile = output_filename)
  
}

Possibly_Converting_OG_to_nt <- possibly(Converting_OG_to_nt, otherwise = "error", quiet = TRUE)

#Make another list w/ OG already in HyPhy inputs OG subdirectory and get OGs that aren't in this list in the OG_list list
#setdiff() - will probably do it

#Running function
map(OG_list, Possibly_Converting_OG_to_nt)



