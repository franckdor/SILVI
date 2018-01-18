# Clean the table and prepare ingredients ---------
source("code/Epitope.R")
silent <- list.files(pattern="csv|fasta")  %>% file.remove()

# CLASS I ----------------------------
PATH_TO_CSV_FOLDER = "data/classI"
cooked <- PATH_TO_CSV_FOLDER %>% cook_I() %T>% print
# save(cooked, file="cooked_temp.rda")
#load("cooked_temp.rda")
# beepr::beep()

# then 
PATH_TO_BLAST_TXT  = "blast_ClassI.txt"
digested <- cooked %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>%
  export_csv("3_blast_mismatches_I") %>% print
# save(digested, file="digested_temp.rda")
# load("digested_temp.rda")
# beepr::beep()

# fix

 # CLASS II ----------------------------
PATH_TO_CSV_FOLDER = "data/classII"
cookedtestII <- PATH_TO_CSV_FOLDER %>% cook_II() %T>% print
# then:
PATH_TO_BLAST_TXT  = "blast_ClassII.txt"
digested <- cookedtestII %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>% export_csv("3_blast_mismatches_II") %>% print


