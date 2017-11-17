# Clean the table and prepare ingredients ---------
source("code/Epitope.R")
silent <- list.files(pattern="csv|fasta")  %>% file.remove()

# CLASS I ----------------------------
PATH_TO_CSV_FOLDER = "data"
cooked <- PATH_TO_CSV_FOLDER %>% cook_I() %T>% print
# save(cooked, file="cooked_temp.rda")
#load("cooked_temp.rda")
# beepr::beep()

# then 
PATH_TO_BLAST_TXT  = "data/blast_results_a7.txt"
digested <- cooked %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>%
  export_csv("3_blast_mismatches_I") %>% print
# save(digested, file="digested_temp.rda")
# load("digested_temp.rda")
# beepr::beep()

# fix

 # CLASS II ----------------------------
PATH_TO_CSV_FOLDER = "data/set2"
cookedtestII <- PATH_TO_CSV_FOLDER %>% cook_II() %T>% print
# then:
PATH_TO_BLAST_TXT  = "data/0Y4H86R4015-Alignment.txt"
digested <- cookedtestII %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>% export_csv("3_blast_mismatches_II") %>% print

PATH_TO_BLAST_TXT  = "seta_dataii-Alignment.txt"
digested %>% select(peptide, blast, middle, starts_with("mismatch"))

