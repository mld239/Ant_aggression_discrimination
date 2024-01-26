library("tidyverse")
library(RJSONIO)
library("googlesheets4")

## Clear Objects ##
rm()

################# Processing BUSTED results #################

BUSTED_results <- list.files("./HyPhy_outputs", full.names = TRUE) %>% grep(pattern = 'relax',  invert = TRUE, value = TRUE)

Processing_BUSTED_Results <- function(input_file){

  result_file <- fromJSON(input_file)

  result_vector <- c(input_file, result_file[["input"]][["file name"]], result_file[["test results"]][["p-value"]], result_file[["test results background"]][["p-value"]], result_file[["test results shared distributions"]][["p-value"]])

  return(result_vector)
}

Possibly_Processing_BUSTED_Results <- possibly(Processing_BUSTED_Results, otherwise = "error", quiet = TRUE)

#Running function
BUSTED_processed_results <- map(BUSTED_results, Possibly_Processing_BUSTED_Results)

#Organizing it better
BUSTED_processed_results <- as.data.frame(do.call(rbind, BUSTED_processed_results))

#Naming column names
colnames(BUSTED_processed_results) <- c("file name w/ trait", "file name", "foreground p-value", "background p-value", "shared p-value")



####### Aligning results to candidates and their orthogroups ######

#Reading in candidates w/ orthogroups
candidate_orthogroups <- read.csv("./candidates_by_orthogroup.csv") %>% distinct()

BUSTED_processed_results$orthogroup <- str_split_i(BUSTED_processed_results$`file name`, pattern = "/", 8) %>% str_split_i(pattern = "_", 1)

#Concatenating information
BUSTED_results_wCandidates <- full_join(BUSTED_processed_results, candidate_orthogroups, by = c("orthogroup" = "OG"))
BUSTED_results_wCandidates$gene_name <- str_split_i(BUSTED_results_wCandidates$gene_name, pattern = "/", 3) %>% str_split_i(pattern = "\\.", 1)

#Reading in candidate info from google spreadsheet
candidate_info <- read_sheet("https://docs.google.com/spreadsheets/d/1OtOGhctFd3ReJdkB22yUGqr28V0COCLB50URQHiSjd8/edit?usp=sharing")
candidate_info$`NCBI Gene ID` <- as.character(candidate_info$`NCBI Gene ID`)

##### Just for my knowledge ####
#Total number of genes
number_of_genes <- unique(candidate_info$`NCBI Gene ID`)
#Genes associated with aggression
number_of_aggr_genes <- filter(candidate_info, `Specific trait` == "Aggression")
number_of_aggr_genes <- unique(number_of_aggr_genes$`NCBI Gene ID`)
#Genes associated with discrimination
number_of_disc_genes <- filter(candidate_info, `Specific trait` == "Nestmate Discrimination")
number_of_disc_genes <- unique(number_of_disc_genes$`NCBI Gene ID`)
#########################################################################

#Concatenating information
BUSTED_results_wCandidates <- full_join(BUSTED_results_wCandidates, candidate_info, by = c("gene_name" = "NCBI Gene ID"))

#Determining selection based on p-values
BUSTED_results_wCandidates <- BUSTED_results_wCandidates %>% mutate(selectionOn =
                                                                      case_when(as.numeric(as.character(`foreground p-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) > 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) <= 0.05 ~ "ForegroundOnly",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) <= 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) <= 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) > 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) > 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) <= 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) <= 0.05 ~ "BackgroundOnly",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) > 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) <= 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) > 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) > 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) <= 0.05 ~ "NoEvidenceOfSelection",
                                                                                
                                                                                as.numeric(as.character(`foreground p-value`)) > 0.05 & 
                                                                                  as.numeric(as.character(`background p-value`)) > 0.05 &
                                                                                  as.numeric(as.character(`shared p-value`)) > 0.05 ~ "NoEvidenceOfSelection"))


BUSTED_results_wCandidates$"trait analysed" <- str_split_i(BUSTED_results_wCandidates$`file name w/ trait`, pattern = "_", 3)

BUSTED_filtered_results <- filter(BUSTED_results_wCandidates, (`Specific trait` == "Aggression" & `trait analysed` == "Aggr.json") | (`Specific trait` == "Nestmate Discrimination" & `trait analysed` == "Disc"))



########## Processing Relax Results #################

Relaxed_results <- list.files("./HyPhy_outputs", full.names = TRUE, pattern = "*relax.json")

Processing_Relaxed_Results <- function(input_file){
  
  relaxed_result <- fromJSON(input_file)

  result_vector <- c(input_file, relaxed_result[["input"]][["file name"]], relaxed_result[["test results"]][["p-value"]], relaxed_result[["test results"]][["relaxation or intensification parameter"]])

  return(result_vector)
}

Possibly_Processing_Relaxed_Results <- possibly(Processing_Relaxed_Results, otherwise = "error", quiet = TRUE)

#Running function
Relaxed_processed_results <- map(Relaxed_results, Possibly_Processing_Relaxed_Results)

#Organizing it better
Relaxed_processed_results <- as.data.frame(do.call(rbind, Relaxed_processed_results))

#Naming column names
colnames(Relaxed_processed_results) <- c("result file name", "input file name", "P-value", "K-value")


####### Aligning results to candidates and their orthogroups ######

#Reading in candidates w/ orthogroups
candidate_orthogroups <- read.csv("./candidates_by_orthogroup.csv") %>% distinct()
#Removing spaces from gene names
candidate_orthogroups$gene_name <- str_remove_all(candidate_orthogroups$gene_name, " ")
#Getting name w/o file extension
candidate_orthogroups$file_name <- str_split_i(candidate_orthogroups$gene_name, pattern = "/", 3) %>% str_split_i(pattern = "\\.", 1)

Relaxed_processed_results$orthogroup <- str_split_i(Relaxed_processed_results$`input file name`, pattern = "/", 8) %>% str_split_i(pattern = "_", 1)

#Concatenating information
Relax_results_wCandidates <- full_join(Relaxed_processed_results, candidate_orthogroups, by = c("orthogroup" = "OG"))
Relax_results_wCandidates$gene_name <- str_split_i(Relax_results_wCandidates$gene_name, pattern = "/", 3) %>% str_split_i(pattern = "\\.", 1)

#Reading in candidate info from google spreadsheet (may not need to do again)
candidate_info <- read_sheet("https://docs.google.com/spreadsheets/d/1OtOGhctFd3ReJdkB22yUGqr28V0COCLB50URQHiSjd8/edit?usp=sharing")
candidate_info$`NCBI Gene ID` <- as.character(candidate_info$`NCBI Gene ID`)

mismatches <- base::setdiff(candidate_orthogroups$file_name, candidate_info$`NCBI Gene ID`)

#Concatenating information
Relax_results_wCandidates <- full_join(Relax_results_wCandidates, candidate_info, by = c("gene_name" = "NCBI Gene ID"))

#Determining selection based on p-values - make this look at foreground p value and k value
Relax_results_wCandidates <- Relax_results_wCandidates %>% mutate(selectionOn =
                                                                      case_when(as.numeric(as.character(`P-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`K-value`)) > 1 ~ "IntensifiedSelection",
                                                                                
                                                                                as.numeric(as.character(`P-value`)) <= 0.05 & 
                                                                                  as.numeric(as.character(`K-value`)) <= 1  ~ "RelaxedSelection",
                                                                                
                                                                                as.numeric(as.character(`P-value`)) > 0.05 ~ "Nonsignificant"))
                                                                                
                                                                               

Relax_results_wCandidates$"trait analysed" <- str_split_i(Relax_results_wCandidates$`result file name`, pattern = "_", 3)

Relax_filtered_results <- filter(Relax_results_wCandidates, (`Specific trait` == "Aggression" & `trait analysed` == "Aggr") | (`Specific trait` == "Nestmate Discrimination" & `trait analysed` == "Disc"))

######## Uploading results to google sheets ###############

BUSTED_filtered_results %>%
  write_sheet(
    ss = gs4_get(
      "https://docs.google.com/spreadsheets/d/1KHr6uqYVAMwYhMFqHlIcLJ2xp5ar7uDhpAEuGmdGBJ4/edit?usp=sharing"),
    sheet = "Sheet1"
  )


Relax_filtered_results %>%
  write_sheet(
    ss = gs4_get(
      "https://docs.google.com/spreadsheets/d/1yrFa-oRUB8nkKbxSgv_WO1C3xlZ3Toz4Kp1Qs6D_f1o/edit?usp=sharing"),
    sheet = "Sheet1"
  )


