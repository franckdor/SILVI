# Clean the table and prepare ingredients ---------
source("code/Epitope.R")
silent <- list.files(pattern="csv|fasta")  %>% file.remove()

# CLASS I ----------------------------
PATH_TO_CSV_FOLDER = "data"
cooked <- PATH_TO_CSV_FOLDER %>% cook_I() %T>% print
# save(cooked, file="cooked_temp.rda")
load("cooked_temp.rda")
# beepr::beep()

# then
PATH_TO_BLAST_TXT  = "setb_datai-Alignment.txt"
digested <- cooked %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>%
  export_csv("3_blast_mismatches_I") %>% print
# save(digested, file="digested_temp.rda")
# load("digested_temp.rda")
# beepr::beep()

# fix
 path <- "data/setb_dataii/b1_iedbii.csv"
import_1_protein_predictor_class2(path)

 # CLASS II ----------------------------
PATH_TO_CSV_FOLDER = "data/setb_dataii_test"
cookedtestII <- PATH_TO_CSV_FOLDER %>% cook_II() %T>% print
# then:
PATH_TO_BLAST_TXT  = "seta_dataii-Alignment.txt"
digested <- cookedaII %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>% export_csv("3_blast_mismatches_I") %>% print

PATH_TO_BLAST_TXT  = "seta_dataii-Alignment.txt"
digested %>% select(peptide, blast, middle, starts_with("mismatch"))

