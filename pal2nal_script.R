library("tidyverse")

#Getting list of file names
input_file <-list.files("./HyPhy_inputs/OG_alignments", pattern = "*.fa", full.names = FALSE)
input_file <- gsub(".fa", "", as.character(input_file))

#Making function
Aligning_OGs <- function(input_file){
  
  pal2nal_input <- paste("./pal2nal.v14/pal2nal.pl  ./aligned_HOG/", input_file,".fa",
                         " ./HyPhy_inputs/OG_alignments/", input_file, ".fa",
                         " -output fasta > ./HyPhy_inputs/OG_alignments/", input_file, "_aligned", ".fa", 
                         sep="")
  system(pal2nal_input)
}

Possibly_Aligning_OGs <- possibly(Aligning_OGs, otherwise = "couldn't open file", quiet = TRUE)

#Running function
map(input_file, Possibly_Aligning_OGs)

