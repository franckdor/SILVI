# https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/
#

library(tidyverse)
library(magrittr) # nécessaire pour éviter une erreur avec  %<>% 
# Install package from CRAN only if not installed, and load the library

if (!require(testthat)) install.packages('testthat')
library(testthat)

# source("../../code/SILVI.R") # This is only needed if your project is not a package
# source("/Users/dorkeld/Documents/old_documents/travail_mac/projets/Code_LOIRE_SILVI/14_oct_test_fin/SILVI/code/SILVI.R")

source("SILVI.R") # contient la fct best_sliding utilisée par process_blast.R
source("process_blast.R")


# PREREQUIS POUR LES TESTS
# lecture d'un fichier fasta de test
# pour faire une dataframe de test
# puis lancement de add_best_From_blast pour avoir un dataframe correcte

FASTA_FILE="../tests_data/1_blastme_I_pf_csp.fasta"

peptide <- c()         #vector() # un vecteur de peptides

con = file(FASTA_FILE, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  #print(line)
  if (! startsWith(line, ">")  ) {
    # print (line)
    peptide <- c(peptide,line)
  }
}

close(con)

peptide_df <- data.frame(peptide, stringsAsFactors=FALSE)
# df$peptide <- factor(df$peptide)
peptide_df <- as_tibble(peptide_df)

cat(" len peptide : ",length(peptide_df$peptide),"\n")
expect_equal(length(peptide_df$peptide), 228)

# lecture  d'un fichier d'alignement après blast
PATH_TO_BLAST_TXT="../tests_data/pf_csp_alignment_I.txt"
add_best_df <- add_best_from_blast(peptide_df,PATH_TO_BLAST_TXT)
print(add_best_df)







#import_all_protein_predictor_class1




# TESTER cook_I
# TESTER cook_II




# test : digest_I_and_II
# on va tester les ss fction individuellemnt avnt de tester tte la fct
# les scores
test_that("Add one to 99", {
  # expect_equal(add_one(99), 100)
  
 
  res_df <- count_mismatches_fd(add_best_df)
  #  PRODUIT  PAR count_mismatches ou count_mismatches_fd appelé par à mettre dans test_SILVI.R
  
  # print(res_df$m_1[5])
  # print(res_df$m_2[5])
  # print(res_df$m_3[5])
  
})

# tester count_mismatch_fd (à remettre dans test_SILVI.R) # ou créer un fichier final_scoring.R

test_that("test des scores False true", {
   
   res_df <- count_mismatches_fd(add_best_df)
   
   # FALSE	TRUE	TRUE	FALSE	TRUE	TRUE	TRUE	TRUE	TRUE
   expect_equal(res_df$m_1[5],FALSE) 
   expect_equal(res_df$m_2[5],TRUE) 
   expect_equal(res_df$m_3[5],TRUE) 
   
   # FALSE	FALSE	TRUE	TRUE	TRUE	FALSE	FALSE	FALSE	FALSE
   expect_equal(res_df$m_1[23],FALSE) 
   expect_equal(res_df$m_2[23],FALSE)  
   expect_equal(res_df$m_3[23],TRUE) 
  
   # FALSE	FALSE	TRUE	FALSE	TRUE	TRUE	TRUE	TRUE	TRUE
   expect_equal(res_df$m_1[42],FALSE)   
   expect_equal(res_df$m_2[42],FALSE) 
   expect_equal(res_df$m_3[42],TRUE) 
 
   # TRUE	TRUE	TRUE	FALSE	FALSE	FALSE	FALSE	FALSE
   expect_equal(res_df$m_1[51],FALSE)     
   expect_equal(res_df$m_2[51],TRUE) 
   expect_equal(res_df$m_3[51],TRUE) 
   
 })

