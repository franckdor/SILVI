# Clean the table and prepare ingredients ---------
source("code/SILVI.R")
require(Peptides)
require(tidyverse)
require(stringr)


################################################
# STEP A 
################################################

PATH_TO_CSV_FOLDER = "dataii_example"
cookedtestII <- PATH_TO_CSV_FOLDER %>% cook_II() %T>% print

# Now you have a fasta file named 1_blastme_II.fasta of conserved epitopes
# Use blastponline against the proteome of your choice. Download and place the resulting text file in the working directory

################################################
# STEP B
################################################

PATH_TO_BLAST_TXT  = "alignment_II.txt"
digested <- cookedtestII %>% digest_I_and_II(PATH_TO_BLAST_TXT) %T>% export_csv("3_blast_mismatches_II") %>% print


#################################################
# STEP C
#################################################

final_ClassII <- merge_result_classII()  

write.table(file="res_classII.csv",final_ClassII,sep=";",row.names=F)
