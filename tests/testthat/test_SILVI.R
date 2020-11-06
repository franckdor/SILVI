# https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/
#

library(testthat)

# source("../../code/SILVI.R") # This is only needed if your project is not a package
# source("/Users/dorkeld/Documents/old_documents/travail_mac/projets/Code_LOIRE_SILVI/14_oct_test_fin/SILVI/code/SILVI.R")



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

#  PRODUIT  PAR count_mismatches ou count_mismatches_fd appelé par à mettre dans test_SILVI.R
#  print(add_best_df$m_1[5])
#  print(add_best_df$m_2[5])
#  print(add_best_df$m_3[5])


# TESTER cook_I
# TESTER cook_II


# digest_I_and_II


# test : digest_I_and_II
# on va tester 
test_that("Add one to 99", {
  # expect_equal(add_one(99), 100)
  
  # res_tible <- digest_I_and_II()

})
