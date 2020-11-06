# https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/
#

library(testthat)

# source("../../code/SILVI.R") # This is only needed if your project is not a package
# source("/Users/dorkeld/Documents/old_documents/travail_mac/projets/Code_LOIRE_SILVI/14_oct_test_fin/SILVI/code/SILVI.R")

#import_all_protein_predictor_class1

#  PRODUIT  PAR count_mismatches ou count_mismatches_fd appelé par à mettre dans test_SILVI.R
#  print(add_best_df$m_1[5])
#  print(add_best_df$m_2[5])
#  print(add_best_df$m_3[5])


# TESTER cook_I
# TESTER cook_II


digest_I_and_II


# test : digest_I_and_II
# on va tester 
test_that("Add one to 99", {
  # expect_equal(add_one(99), 100)
  
  # res_tible <- digest_I_and_II()

})
