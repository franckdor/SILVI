# Clean the table and prepare ingredients ---------
source("code/Epitope.R")
silent <- list.files(pattern="csv|fasta")  %>% file.remove()

 # CLASS II ----------------------------
PATH_TO_CSV_FOLDER = "data/classII"
cookedtestII <- PATH_TO_CSV_FOLDER %>% cook_II() %T>% print
# then:
PATH_TO_BLAST_TXT  = "blast_ClassII.txt"
digested <- cookedtestII %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>% export_csv("3_blast_mismatches_II") %>% print

source("FilterAndRankPeptides_ClassII.R")
