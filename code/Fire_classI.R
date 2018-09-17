require(Peptides)
require(tidyverse)
require(stringr)


# Clean the table and prepare ingredients ---------
source("code/SILVI.R")

################################################
# CLASS I ----------------------------
################################################

################################################
# STEP A 
################################################

# Here add the path to prediction files for CLASS I epitopes:
PATH_TO_CSV_FOLDER = "datai_example"
cooked <- PATH_TO_CSV_FOLDER %>% cook_I() %T>% print


# Now you have a fasta file named 1_blastme_I.fasta of conserved epitopes
# Use blastponline against the proteome of your choice. Download and place the resulting text file in the working directory


################################################
# STEP B (after blast)
################################################

# Path and name of the blastp results:
PATH_TO_BLAST_TXT  = "alignment_I.txt"

digested <- cooked %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>%
  export_csv("3_blast_mismatches_I") %>% print


################################################
# STEP C 
################################################

final_Class1 <- merge_result_classI()  

write.table(file="res_classI.csv",final_Class1,sep=";",row.names=F)
