source("code/SILVI.R")
check_csvs <- function(PATH_TO_CHECK, classIorII){
  # domestic functions that prints summaries
  spit_diag <- function(tf){
    if (all(tf))
      cat("\tOK\n")
    else {
      cat("\tNOT OK because of:\n\t")
      cat(pretty_lf[which(!tf)], sep="\n\t")
    }
    sum(!tf)
  }

  # domestic function to build db_columns
  # check for allowed columns
  # bloody dirty
  build_db_df <- function(){
    dbl <- "code/db_headers" %>%
      readLines() %>%
      strsplit(";")

    db_columns <- dbl %>%
      sapply(`[`, -(1:2)) %>%
      data_frame()

    db_map <- dbl %>% sapply(`[`, 1:2) %>%
      t %>%
      as_data_frame()

    db_df <- cbind(db_map, db_columns)
    colnames(db_df) <- c("class", "db", "columns")
    db_df$class %<>% factor %>% as.numeric
    db_df
  }
  db_df <- build_db_df()
  db_df <- db_df[db_df$class==classIorII, ]

  # short names for csv files
  pretty_lf <- PATH_TO_CHECK %>% list.files(full.names=FALSE)

  # a list of all reads lines
  rls <- PATH_TO_CHECK %>%
    list.files(full.names=TRUE) %>%
    lapply(function(lf1) lf1 %>% readLines(n=2))

  # a list of data.frame attempts
  dfs <- rls %>% lapply(function(rls1) read.csv2(text=rls1, sep=";", dec="."))

  # then a stack of tests
  err <- 0

  cat("* Testing that the number of columns match the nb of data per row :")
  err <- err + ((rls %>% sapply(function(rl1) rl1 %>% sapply(strsplit, ";") %>%
                                 sapply(length) %>% unique %>% length() %>% `==`(1))) %>% spit_diag)

  cat("* Testing if all data.frame can be read : ")
  err <- err + (dfs %>% sapply(is.data.frame) %>% spit_diag())

  cat("* Testing if filenames have proper db encoding (eg ..._iedb.csv) : ")
  match_db_found <-  match(seq_along(pretty_lf), grep(paste(db_df$db, collapse="|"), pretty_lf))
  err <- err + (match_db_found %>% is.na %>% `!` %>% spit_diag())

  cat("* Testing if columns are nicely detected (;) : ")
  err <- err + (dfs %>% sapply(function(df1) ncol(df1) != 1)  %>% spit_diag())

  cat("* Testing if all required columns are present (as declared in code/db_headers) : ")
  # bloody dirty
  pretty_db <- pretty_lf %>% gsub("\\.csv", "", .) %>% strsplit("_") %>% sapply(`[`, 2)
  res <- vector("logical", length(pretty_lf))
  for (i in seq_along(pretty_lf)){
    exp <- db_df[db_df$db==pretty_db[i],]$columns %>% unlist
    obs <- dfs[[i]] %>% colnames()
    res[i] <- match(obs, exp) %>% is.na() %>% any() %>% `!`
  }
  err <- err + (res %>% spit_diag())

  # Final error counting
  cat(err, "errors found testing all .csv files in", PATH_TO_CHECK)
}

check_csvs(PATH_TO_CHECK = "data/setb_dataii",
           classIorII = 2)

check_csvs(PATH_TO_CHECK = "data/setb_datai_vb_cleanup/",
           classIorII = 1)

# Helpers (very helpful but very dangerous, work on a copy if you're not sure)
# if you need some changes to be done, globally on a folder, here are some recipes
# first change the path to fix below, then run one or more of the recipes below
PATH_TO_FIX="data/setb_datai_vb_cleanup"

rewrite_all_csvs_to_semicolons(PATH_TO_FIX)
rewrite_all_csvs_to_same_nb_data_as_nb_columns(PATH_TO_FIX)

################################################################################
# Recipe to change the "...","...","..."  and ...(tab)...(tab) format to ...;...;...
# in all files within a folder

rewrite_all_csvs_to_semicolons <- function(PATH_TO_FIX){
  lf <- list.files(PATH_TO_FIX, full=T)
  for (i in seq_along(lf))
    lf[i] %>% readLines() %>% gsub("\"", "", .) %>% gsub(",|\t", ";", .) %>% writeLines(con=lf[i])
}

################################################################################
# Recipe to remove any additional column on the right
# eg if the nb of data is different than the nb of columns

rewrite_all_csvs_to_same_nb_data_as_nb_columns <- function(PATH_TO_FIX){
  lf <- list.files(PATH_TO_FIX, full=T)

  remove_last_column_from_second_line_to_bottom <- function(x){
    rm_right <- function(x)
      x %>% strsplit(";") %>% sapply(function(l1) l1[1:(length(l1)-1)]) %>% paste(collapse=";")
    if (((x[1] %>% strsplit(";") %>% unlist %>% length())+1) == (x[2] %>% strsplit(";") %>% unlist %>% length()))
      c(x[1], lapply(x[-1], rm_right) %>% unlist)
    else
      x
  }

  for (i in seq_along(lf))
    lf[i] %>% readLines() %>% remove_last_column_from_second_line_to_bottom() %>% writeLines(con=lf[i])
}

