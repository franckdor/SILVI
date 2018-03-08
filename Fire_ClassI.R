# Clean the table and prepare ingredients ---------
source("code/Epitope.R")
silent <- list.files(pattern="csv|fasta")  %>% file.remove()

# CLASS I ----------------------------
PATH_TO_CSV_FOLDER = "data/classI"
cooked <- PATH_TO_CSV_FOLDER %>% cook_I() %T>% print

# then 
PATH_TO_BLAST_TXT  = "blast_ClassI.txt"
digested <- cooked %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>%
  export_csv("3_blast_mismatches_I") %>% print

source("FilterAndRankPeptides_ClassI.R")

# STEP3
