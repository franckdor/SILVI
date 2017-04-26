# Preliminaries -------------------------------------------

# rinse environmment
rm(list=ls())

# dependencies
library(plyr)
library(dplyr)
library(magrittr)

# Final pipes ---------------------------------------------
cook_I <- function(PATH_TO_CSV_FOLDER){
  PATH_TO_CSV_FOLDER %>%
    import_all_protein_predictor_class1() %>%
    map_supertypes_alleles("data/map_supertypes_alleles.csv", remove_these_alleles=c("", "B4101"))  %>%
    # only_peptides_repeated_seqnum_times()  %>%
    common_among_predictors_I() %>%
    reorder_column("supertype") %>%
    reorder_column("predictor") %T>%
    export_csv("1_common_I") %T>%
    export_FASTA() %>%
    return()
}

cook_II <- function(PATH_TO_CSV_FOLDER){
  PATH_TO_CSV_FOLDER %>%
    import_all_protein_predictor_class2() %>%
    common_among_predictors_II() %T>%
    export_csv("1_common_II") %T>%
    export_FASTA() %>%
    return()
}

digest_I_and_II <- function(cooked, PATH_TO_BLAST_TXT){
  cooked %>%
    add_best_from_blast(PATH_TO_BLAST_TXT) %T>%
    export_csv("2_common_blast_I") %>%
    count_mismatches() %>%
    after_two_spaces_all_mismatches() %>%
    final_polish() %>%
    return()
}

# Domestic functions --------------------------------------

# import class1 ====================
# given a .csv path, import them as data_frame
#   works whatever .csv use  "," or ";" as field separators
#   it would have been read.table(path, h=T, sep=",") if all files used ",".
import_1_protein_predictor_class1 <- function(path){
  # Rcpp solution
  # df <- read_file_cpp2(path) %>% gsub(";", ",", .) %>% textConnection()  %>% read.table(sep=",", header=T)
  # noRcpp solution
  # df <- readLines(path) %>% gsub(";", ",", .) %>% textConnection()  %>% read.table(sep=",", header=T)
  # now with homogeneous files
  df <- read.csv2(path)
  # retrieve and add protein name from path (precarious)
  protein <- strsplit(path, "/|\\.|_")  %>% sapply(function(.) .[length(.)-2])
  predictor <- strsplit(path, "/|\\.|_")  %>% sapply(function(.) .[length(.)-1])
  df %<>% mutate(protein=protein, predictor=predictor, file=path)
  # rename the score column into something consistent accross files
  colnames(df)[grep("ann_ic50|ic50|score|nm", colnames(df))] <- "score"
  # allele as character
  df %<>% mutate(allele=allele %>% as.character())
  # reorder columns
  df %>% select(protein, allele, predictor, seq_num, score, peptide, file) %>%
    # we filter seqnum
    only_all_seqnum() %>%
    # we refactor peptide to drop missing levels
    mutate(peptide=factor(peptide), score=score %>% as.character %>% as.numeric) %>%
    # and return a pretty data_frame
    # arrange(allele, seq_num, peptide) %>%
    as_data_frame()
}

# print summary
how_many_rows <- function(df, more)
  cat("*", nrow(df), "peptides", ifelse(missing(more), "\n", paste0(more, "\n")))

# imports a list of protein_predictor.csv
import_all_protein_predictor_class1 <- function(dir.path){
  # list all files paths
  lf_csv <- list.files(dir.path, pattern="csv", full.names=TRUE)
  # import all of them and return as a list
  # lapply(lf_csv, import_1_protein_predictor) would be an option
  # but we prefer a loop to spit messages and quickly identify problems with files
  res <- vector("list", length(lf_csv))
  cat("\n")
  cat("* Importing:", length(lf_csv), ".csv files from", dir.path, "\n")
  for (i in seq_along(lf_csv)){
    cat("\t-> ", lf_csv[i], "   ")
    res[[i]] <- import_1_protein_predictor_class1(lf_csv[i])
    cat("OK -", length(lf_csv)-i,"remaining\n")
  }
  # merge all dfs
  res %>% do.call(rbind, .) %T>% how_many_rows("imported")
  # return(res)
}

# import class2 ====================
# given a .csv path, import them as data_frame
#   works whatever .csv use  "," or ";" as field separators
#   it would have been read.table(path, h=T, sep=",") if all files used ",".
import_1_protein_predictor_class2 <- function(path){
  # df <- readLines(path)  %>%  gsub(";", ",", .) %>% textConnection()  %>% read.table(sep=",", header=T)
  # now with homogeneous files
  df <- read.csv2(path)
  # rename some columns, and create empty others to allow further do.call(rbind)
  netmhc_id <- df  %>% colnames  %>% grep("^C|core$", .)
  netmhc_flag <- netmhc_id %>% any()

  # if (netmhc_flag){
  #   df  %<>% rename(full_peptide=Peptide, peptide=Core, seq_num=Identity, allele=Allele) %>%
  #     mutate(start=NA, end=NA, percentile_rank=NA, score=NA)
  # } else {
  #   df  %<>% rename(full_peptide=peptide) %>% rename_(peptide=colnames(df)[netmhc_id])
  #     mutate(start=NA, end=NA, percentile_rank=NA, Pos=NA, Affinity.nM.=NA, X.Rank=NA, BindingLevel=NA, score=NA)
  # }
  df %<>% rename(peptide=core, score=ic50) %>% mutate
  # retrieve and add protein name from path (precarious)
  protein <- strsplit(path, "/|\\.|_")  %>% sapply(function(.) .[length(.)-2])
  predictor <- strsplit(path, "/|\\.|_")  %>% sapply(function(.) .[length(.)-1])
  df %<>% mutate(protein=protein, predictor=predictor, file=path)
  # rename the score column into something consistent accross files
  # colnames(df)[grep("ann_ic50|ic50|score|nm", colnames(df))] <- "score"
  if ("nn_align_ic50" %in% colnames(df)){
    df %<>% mutate(score=nn_align_ic50)
  }
  # allele as character
  df %<>% mutate(allele=allele %>% as.character())
  # reorder columns
  df %>%
    # select(protein, allele, predictor, seq_num, full_peptide, peptide,
                # start, end, percentile_rank, Pos, Affinity.nM., X.Rank, BindingLevel, score, file) %>%
    # we filter seqnum
    only_all_seqnum_II() %>%
    # score is sometimes a factor
    mutate(score=score %>% as.character %>% as.numeric) %>%
    # and return a pretty data_frame
    as_data_frame()
}

# imports a list of protein_predictor.csv
import_all_protein_predictor_class2 <- function(dir.path){
  # list all files paths
  lf_csv <- list.files(dir.path, pattern="csv", full.names=TRUE)
  # import all of them and return as a list
  # lapply(lf_csv, import_1_protein_predictor) would be an option
  # but we prefer a loop to spit messages and quickly identify problems with files
  res <- vector("list", length(lf_csv))
  cat("\n")
  cat("* Importing:", length(lf_csv), ".csv files from", dir.path, "\n")
  for (i in seq_along(lf_csv)){
    cat("\t↓ ", lf_csv[i], "   ")
    res[[i]] <- import_1_protein_predictor_class2(lf_csv[i])
    cat("\t✓\n")
  }
  # merge all dfs
  res %>% do.call(rbind, .) %T>%
    how_many_rows("imported")
}
# map with supertypes ===========

map_supertypes_alleles <- function(df,
                                   map.path="data/map_supertypes_alleles.csv",
                                   remove_these_alleles=NULL){
  # tidy %$£ allele names
  tidy_allele <- function(allele){
    allele %>%
      gsub("HLA-", "", .) %>% gsub("\\:", "", .) %>%
      gsub("_", "", .) %>% gsub("\\/", "-", .) %>%
      gsub("\\*", "", .) %>%  gsub("[[:space:]]", "", .) %>%
      gsub("\\(.*\\)", "", .)
  }
  map <- read.csv2(map.path, header=TRUE)
  map %<>% mutate(allele = tidy_allele(allele))
  df  %<>% mutate(allele = tidy_allele(allele))
  df %<>%
    # mutate(allele=tidy_allele(allele))  %>%
    left_join(map, by="allele")  %>%
    mutate(supertype=as.character(supertype)) %>%
    select(protein, supertype, allele, peptide, everything())
  # hook for the manual removing of some alleles
  remove_ids <- which(df$allele %in% remove_these_alleles)
  remove_n   <- sum(remove_ids)
  if (remove_n == 0){
    cat("* No allele removed\n")
  } else {
    remove_these_alleles[nchar(remove_these_alleles)==0] <- "''"
    cat("*", length(remove_ids), "alleles removed (", paste(remove_these_alleles, collapse=", "), ")\n")
    df %<>% slice(-remove_ids)
  }

  nas <- is.na(df$supertype) %>% sum
  if (nas == 0){
    cat("* All alleles mapped with their supertypes\n")
  } else {
    cat("*", nas, "alleles mismatched when mapping their supertypes:\n\t")
    df %>% filter(is.na(supertype)) %$% table(allele, predictor, file) %>% print
    cat("\n")
    df$supertype[is.na(df$supertype)] <- "NO_MAPPING"
  }
  df
}

# filter ====================

# seqnum filtering ---------
only_all_seqnum <- function(df){
  df  %>% group_by(peptide, allele)  %>% filter(n() == length(unique(df$seq_num)))  %>% ungroup()
}

only_all_seqnum_II <- function(df){
  retain <- df %$% table(peptide, seq_num) %>% apply(1, function(.) all(. > 0)) %>% which %>% names
  df2 <- filter(df, peptide %in% retain)
  df2$peptide %<>% droplevels()
  df2
}

# filter peptides present more than once among predictors
# and many other things
common_among_predictors_I <- function(df){
  df %<>%
    # only retains the first letter of the predictor
    mutate(predictor=predictor %>% as.character %>% substr(1, 1) %>% toupper()) %>%
    # retain only peptide among at least two tables
    group_by(peptide) %>% filter(n()>=2)  %>%
    # cosmectics
    arrange(peptide)
  # we erase supertype information in S (always returning all of them
  # yet we need this additional information)
  df$supertype <- paste0(df$predictor, "_", df$supertype)
  df[df$predictor=="S", "supertype"] <- NA

  df$score <- paste0(df$predictor, "_", df$score)
  df[df$predictor=="S", "score"] <- NA
  # just turn commas into points
  df$score %<>% gsub(",", ".", .)
  df %<>%
    # remember where in predictor
    mutate(predictor=paste(unique(predictor), collapse="+"),
           seq_num=paste(unique(seq_num), collapse="+"),
           score=paste(unique(score), collapse="+"),
           supertype=paste(unique(supertype), collapse="+"))

  # reclean supertype and score column
  df$supertype %<>% gsub("\\+NA", "", .) %>% gsub("NA", "", .)
  df$score %<>% gsub("\\+NA", "", .) %>% gsub("NA", "", .)

  df %>%
    # take only one peptide
    slice(1) %>% filter(nchar(predictor) != 1) %>% ungroup() %>%
    # drops predictor
    select(-allele) %>%
    # peptide finally as a character
    mutate(peptide=as.character(peptide)) %>%
    # return this beauty
    return()
}

# filter peptides present more than once among predictors
# and many other things
common_among_predictors_II <- function(df){
  df %<>%
    # only retains the first letter of the predictor
    mutate(predictor=predictor %>% as.character %>% substr(1, 1) %>% toupper()) %>%
    # retain only peptide among at least two tables
    group_by(peptide) %>% filter(n()>=2)  %>%
    # cosmectics
    arrange(peptide)
  # we erase supertype information in S (always returning all of them
  # yet we need this additional information)
  # df$supertype <- paste0(df$predictor, "_", df$supertype)
  # df[df$predictor=="S", "supertype"] <- NA

  df$score <- paste0(df$predictor, "_", df$score)
  df[df$predictor=="S", "score"] <- NA
  # df[df$predictor=="N", "score"] <- NA
  # just turn commas into points
  df$score %<>% gsub(",", ".", .)
  df %<>%
    # remember where in predictor
    mutate(predictor=paste(unique(predictor), collapse="+"),
           seq_num=paste(unique(seq_num), collapse="+"),
           allele=paste(unique(allele), collapse="+"))

  # reclean supertype and score column
  df$score %<>% gsub("\\+NA", "", .) %>% gsub("NA", "", .)

  df %>%
    # take only one peptide
    # slice(1) %>%
    filter(nchar(predictor) != 1) %>% ungroup() %>%
    # peptide finally as a character
    mutate(peptide=as.character(peptide)) %>%
    # we drop some useless columns
    # select(-Pos, -Affinity.nM., -X.Rank, -BindingLevel, -percentile_rank) %>%
    # we want a single pair of full_peptide x peptide
    group_by(full_peptide, peptide) %>%
    # we collapse scores here
    mutate(score=paste(unique(score), collapse="+")) %>%
    slice(1) %>% ungroup %>%
    # and we return this beauty
    return()
}

# simply handles predictor, or any other column, way of collapsing to make appearent that
# S+I is strictly equivalent to I+S
reorder_column <- function(df, column){
  df[, column] %<>% unlist() %>% strsplit("\\+") %>% lapply(sort) %>% lapply(paste0, collapse="+") %>% unlist()
  df
}

# export ===============

export_csv <- function(x, name){
  write.csv(x, file=paste0(name, ".csv"), row.names = FALSE)
  cat("* Intermediate file written: ", paste0(getwd(), "/", name, ".csv\n"))
}

export_FASTA <- function(x){
  paste0(">id", 1:length(x$peptide), "\n", x$peptide) %>%
    writeLines(con="1_blast_me.fasta")
  cat(">>> Blast the file: ", paste0(getwd(), "/1_blast_me.fasta\n"))
}

# after Blast =============
add_best_from_blast <- function(df, blasted.path){
  blasted <- readLines(blasted.path)
  # partition queries
  start_lines <- grep("Query=", blasted)
  end_lines   <- c((start_lines[-1]-1), length(blasted))
  if (length(start_lines)!= length(end_lines))
    stop("start_lines and end_lines are not of the same length")
  if (nrow(df) != length(start_lines))
    message("Note that fasta file and blasted files do not contain the same number of sequences (", nrow(df), "vs. ", length(start_lines), ")")

  # only fills available ids
  if (
    blasted[start_lines] %>% strsplit(" ") %>% sapply(function(x) x[length(x)]) %>%
    gsub("id", "", .) %>% as.numeric %>% diff %>% table %>% length() %>% equals(1)
  ){
    message("removing ", nrow(df) - length(start_lines), " PEPTIDES")
    # start_lines <- start_lines[1:nrow(df)]
    # end_lines   <- end_lines[1:nrow(df)]
    df <- df[1:length(start_lines), ]
  } else {
    stop("something is wrong with the input ids")
  }

  queries <- vector("list", length(start_lines))
  for (i in seq_along(start_lines)){
    queries[[i]] <- blasted[start_lines[i]:end_lines[i]]
  }
  # protein name is located between the line(s) where we have the first ">"
  # and the line before the second 'Length='
  # if more than one line, we collate them
  # and remove the ">"
  protein_name <- function(query){
    query[grep(">", query)[1] : (grep("^Length=", query)[2]-1)]  %>%
      paste0(collapse=" ") %>%
      gsub(">", "", .)
  }

  middle_alignment <- function(query) {
    best_score_line <- grep(" Score = ", query)[1]
    query[best_score_line+4]  %>% trimws()
  }
  middle <- queries %>% sapply(middle_alignment)

  # finds the top line of each alignment
  best_alignment <- function(query) {
    best_score_line <- grep(" Score = ", query)[1]
    # we grab query line and subject
    list(query_line   = query[best_score_line+3]  %>% trimws(),
         middle_line =  query[best_score_line+4])
  }

  regrow_query2 <- function(query, n_expected = 9){
    query_info <- query$query_line %>% strsplit(" +") %>% unlist %>% `[`(2:3)
    query_info[2] <- query$middle_line %>% gsub(" ", "", .)
    left <- paste0(rep("-", as.numeric(query_info[1])-1), collapse ="")
    middle <- query_info[2] %>% gsub("-", "", .)
    n_missing_right <- n_expected - nchar(left) - nchar(middle)
    if (n_missing_right>0)
      right <- paste0(rep("-", n_missing_right), collapse ="")
    else
      right <- ""
    final <- paste0(left, middle, right)
    return(final)
  }

  df %<>% mutate(middle=middle,
                 blast=queries %>%  lapply(best_alignment) %>% lapply(regrow_query2)  %>% unlist(),
                 blast_info=queries %>% sapply(protein_name))

  # turns into character and replace spaces with dashes
  fac_to_dashes <- function(x) x %>% as.character() %>% gsub(" ", "-", .)
  # use it to polish peptide columns
  df %<>% mutate(peptide = fac_to_dashes(peptide),
                 middle  = fac_to_dashes(middle),
                 blast   = fac_to_dashes(blast))
  # adds a new column a reslides any gaps
  df %<>% mutate(middle2=NA)
  for (i in 1:nrow(df)){
    df$middle2[i] <- best_sliding(df$blast[i], df$middle[i])
  }
  cat("* Best blasts successfully merged\n")
  return(df)
}

# given a fixed and a sliding strings, tries to add gaps
# on the left and or on the right and return the best alignment
best_sliding <- function(fixed, sliding){
  # if of same length simply return original string
  if (nchar(fixed)<=nchar(sliding))
    return(sliding)
  # explode strings
  fixed %<>% strsplit("") %>% unlist()
  sliding %<>% strsplit("") %>% unlist()
  # calculates length difference
  n <- abs(length(fixed) - length(sliding))
  # otherwise, calculate every way to distribute spaces
  df <- expand.grid(left=0:n, right=0:n)
  # retain only those summing to n
  df <- df[rowSums(df)==n, ]
  # prepare a numeric to store results
  res <- numeric(nrow(df))
  for (i in 1:nrow(df)){
    sliding_i <- c(rep("-", df[i, ]$left), sliding, rep("-", df[i, ]$right))
    res[i] <- sum(fixed == sliding_i[1:length(fixed)])
  }
  # get the best sliding combination
  res <- df[which.max(res),]
  # finally append spaces, collate and return it
  c(rep("-", res$left), sliding, rep("-", res$right)) %>%
    paste0(collapse="") %>% return()
}

# mismatches ==================
count_mismatches <- function(df){
  l1 <- df$peptide %>% strsplit("")
  l2 <- df$middle2 %>% strsplit("") %>% lapply("[", 1:9)
  m_df <- mapply(function(x, y) (x == y) | (y=="+"), l1, l2)  %>% t() %>% as.data.frame()
  colnames(m_df) <- paste0("m_", 1:9)
  return(cbind(df, m_df))
}

after_two_spaces_all_mismatches <- function(df){
  # domestic function for 1 line
  after_two_spaces_all_mismatches_1 <- function(x){
    # no more than 1 space
    if (x$middle %>% grep(" {2,}", .)  %>% length()  %>% equals(0)) {
      return(x)
    } else {
      # that's were the space begins
      left <- x$middle %>% strsplit("") %>% extract2(1) %>% equals(" ") %>% which() %>% min()
      x[grep("m_", names(x))[-(1:left)]] <- TRUE
      return(x)
    }
  }
  # then we apply it
  for (i in 1:nrow(df))
    df[i, ] <- after_two_spaces_all_mismatches_1(df[i, ])
  return(df)
}

# just reordering and renaming some columns
final_polish <- function(df){
  # quick and very dirty
  if( "supertype" %in% colnames(df)) {
    # reorder columns and rename blast on the fly
    df %>%
      select(peptide, middle, middle2, subject=blast, starts_with("m_"), protein, supertype, predictor, seq_num, score, blast_info, file) %>%
      as_data_frame() %>% return()
  } else {
    df %>%
      select(peptide, middle, middle2, subject=blast, starts_with("m_"), protein, predictor, seq_num, score, blast_info, file) %>%
      as_data_frame() %>% return()
  }
}
