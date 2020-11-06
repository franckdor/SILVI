#
# test_process_blast.R
# test des fonctionnalités du fichier process_blast_R
# remplace t_process_blast.R : à détruire
#
#
# Utilisation:
#
# se mettre au niveau du fichier code.
#
# source("../tests/testthat/test_process_blast.R")
#
# Ref Voir :
#        Doc sur mac
#        https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/
#

# il faut être au niveau du ss rep  code pour effectuer les tests
# setwd("/Users/dorkeld/Documents/old_documents/travail_mac/projets/Code_LOIRE_SILVI/Debuggage_alignements/code")

# source("../../code/SILVI.R") # This is only needed if your project is not a package
# source("/Users/dorkeld/Documents/old_documents/travail_mac/projets/Code_LOIR"E_SILVI/Debuggage_alignements/code/t_process_blast.R")


source("SILVI.R") # contient la fct best_sliding utilisée par process_blast.R
source("process_blast.R")


library(tidyverse)
library(magrittr) # nécessaire pour éviter une erreur avec  %<>% 
# Install package from CRAN only if not installed, and load the library
if (!require(testthat)) install.packages('testthat')
library(testthat)


  
  # lecture d'un fichier fasta de test
  # pour faire une dataframe de test
  
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
  # print(add_best_df$fd_middle2[1:10])
  # print(add_best_df$fd_middle2[23:42])
  # print(add_best_df$fd_middle2[51])



# ---------------------------------------------------
# "test des sorties de add_best_from_blast"  
  
  # utilser expect_setequal et comparer avec peptide_df$peptide
  test_that("test de recup des peptides mettre dans cooked", {
    
    expect_equal(add_best_df$peptide[5],"ANANSAVKN") 
    expect_equal(add_best_df$peptide[23],"DIEKKICKM")  
    expect_equal(add_best_df$peptide[42],"EPSDKHIKE")
    expect_equal(add_best_df$peptide[51],"GHNMPNDPN")
  })  
  
test_that("test de fd_middle2 completion à droite", {

  expect_equal(add_best_df$fd_middle2[5],"-NA-SA+KN") 
  expect_equal(add_best_df$fd_middle2[23],"--EKK--ICKM")  
  expect_equal(add_best_df$fd_middle2[42],"--S-KHIKE")
  expect_equal(add_best_df$fd_middle2[51],"-HN+-PNDPN")
  

})

