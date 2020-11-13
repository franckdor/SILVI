# File : process_blast.R
require(stringr)
# source("best_sliding.R")

######## library(tidyverse) # fd à enlever

# process blastp output saved as text file
#

# in : df = cooked (tible)
# in : blasted.path (full path for blast)

# process blast output
# complete cooked tible with columns
# for each query peptide that match
# retrieve information about subject sequence
# retrieve

# 1 reads full blast files readlines output is (tibble/df)

add_best_from_blast <- function(df, blasted.path){

# first step :

  blasted <- readLines(blasted.path) # store all blast file in memory

  #  uses blasted struct to split queries (alignts output...) in blocs
  start_lines <- grep("Query #", blasted) # each alignment with new  query starts this way (new blast 2.10 ouptput)
  end_lines   <- c((start_lines[-1]-1), length(blasted))

  if (length(start_lines)!= length(end_lines))
    stop("start_lines and end_lines are not of the same length : AJOUT FD")

  if (nrow(df) != length(start_lines))
    message("Note that fasta file and blasted files do not contain the same number of sequences (", nrow(df), "vs. ", length(start_lines), ")")

  # FD  ccode ajouté

  print(
    blasted[start_lines] %>% strsplit(" ") %>% sapply(function(x) x[3])  %>%
    gsub("id", "", .) %>% as.numeric %>% diff %>% table %>% length() %>% equals(1)
      )

  # FD fin  code ajouté



  ##### PATCH TO ONLY FILL AVAILABLE IDS
  #  (FD  comment : les id sont continus, écart de 1 => à la fin test que tt equals 1)
  # on utilise == à la fin car testthat n'aime pas equals qui serait mieux
  if (
    #blasted[start_lines] %>% strsplit(" ") %>% sapply(function(x) x[length(x)]) %>%
    #gsub("id", "", .) %>% as.numeric %>% diff %>% table %>% length() %>% equals(1)

    blasted[start_lines] %>% strsplit(" ") %>% sapply(function(x) x[3])  %>%
    gsub("id", "", .) %>% as.numeric %>% diff %>% table %>% length() == 1

  ){
    message("PATCHING AND REMOVING ", nrow(df) - length(start_lines), " PEPTIDES")
    # start_lines <- start_lines[1:nrow(df)]
    # end_lines   <- end_lines[1:nrow(df)]
    df <- df[1:length(start_lines), ]
  } else {
    stop("SOMETHING IS WRONG")
  }
  #### PATCH ENDS

  vector("list", length(start_lines)) # ??????????? redondant

  queries <- vector("list", length(start_lines))
  for (i in seq_along(start_lines)){
    queries[[i]] <- blasted[start_lines[i]:end_lines[i]]

  }

# second part : analyse of query blocs and score alignments

  # inner functions to parse protein name (query line)
  # and middle alignt sequence
  # called later

  protein_name <- function(query){
    line_idx=grep(">", query)[1]
    query_line=query[line_idx]  %>% gsub(">", "", .)
    seq_id_line=query[line_idx + 1 ] %>% gsub("Sequence ID:\\s+","", .) %>% gsub("\\s+Length:\\s+\\d+","", .)
    paste0(seq_id_line,query_line)
  }

  # proccess one query record (to get the first alignt data)
  # find line starting first alignt : Score:22.7 bits(46), Expect:6.5,
  # return the 5 line ()
  middle_alignment <- function(query) {
    best_score_line <- grep("^Score:", query)[1] # get the first
    query[best_score_line+5]  %>% trimws()
  }
  # call the inner function and create a char vect/list with middle line
  # for all the query
  # get a vector of middle qlignment seq
  middle <- queries %>% sapply(middle_alignment)


  #  print(typeof(middle_fd))
  #  print(middle)

  # nlle fct : corrigée
  # in  : un bloc de query (vecteur de lignes) concerne  une seq quey
  # out : a list, query/middle
  #
  # grep : returns a vector of the matched elements indices (idx of line)
  # ici on cherche la première ligne de l'alignt
  # et la ligne du milieu
  # on rend une liste (query,middle)
  best_alignment_ORI <- function(query) {
    best_score_line <- grep("^Score:", query)[1]
    #print(best_score_line)
    # we grab query line and subject
    list(query_line   = query[best_score_line+4]  %>% trimws(),
         middle_line  = query[best_score_line+5])
  }

#  best_alignment_FD <- queries[1] %>% sapply(best_alignment)  # test fd
#  print(best_alignment_FD)

  # function to parse first alignment : query,middle and subject line
  # in : query : a bloc of lines concerning a query (Query #) and its alignemants
  # out : a list of token from query line (start and stop pos.)
  best_alignment <- function(query) {

    best_score_line <- grep("^Score:", query)[1]
    st_query =query[best_score_line+4]

    # regex to split query line in simple elements : ie  Query  2    SSSNTRVL  9
    res_expr <- str_match( st_query,"Query[:space:]*([:digit:]*)[:space:]*([:alpha:]*)[:space:]*([:digit:]*)"  )
    # cat(" ",res_expr[,2])
    # we grab query line and subject
    # add analyse of query_line
    # list returned
    list(query_line    = query[best_score_line+4]  %>% trimws(),
         query_start   = res_expr[,2],
         query_aligned = res_expr[,3],
         query_end     = res_expr[,4],
         middle_line   = query[best_score_line+5])
  }


  # inner function :
  # left completion with dashes : if query-start > 1
  # right completion with dashes if lg < n_expected
  # in : query a list with two elements : (query_line,query_middle)
  # in : n_expected
  # out a string :
  regrow_query2 <- function(query, n_expected = 9){
    query_info    <- query$query_line %>% strsplit(" +") %>% unlist %>% `[`(2:3)
    query_info[2] <- query$middle_line %>% gsub(" ", "", .)

    left <- paste0(rep("-", as.numeric(query_info[1])-1), collapse ="")
    middle <- query_info[2] %>% gsub("-", "", .)                          # on enlève les espaces de middle

    n_missing_right <- n_expected - nchar(left) - nchar(middle)
    if (n_missing_right>0)
      right <- paste0(rep("-", n_missing_right), collapse ="")
    else
      right <- ""

    final <- paste0(left, middle, right)

    return(final)
  }


  # on modifie cooke tibble with mutate : nlles colonnes

  # use regrow_query2 function
  df %<>% mutate(middle=middle,
                 blast=queries %>%  lapply(best_alignment) %>% lapply(regrow_query2)  %>% unlist(),  #
                 blast_info=queries %>% sapply(protein_name))


  # ----  process middle to add left dash :  ajout fd ****
  l_fd_alg_info <- queries %>%  lapply(best_alignment) # %>% print()

  # cat("\nType === ",typeof(l_fd_alg_info))
  # print(l_fd_alg_info[2])
  # cat("\nv_fd_alg_info\n",l_fd_alg_info$query_line)
  # cat("\nv_fd_alg_info[1]$query_start ",l_fd_alg_info[[1]]$query_start)
  # cat("\n->names",names(l_fd_alg_info))
  # cat("\n",lengths(l_fd_alg_info))
  # cat("\nv_fd_alg_info\n",l_fd_alg_info$query_start)
  # cat("\nv_fd_alg_info\n",l_fd_alg_info$query_aligned)
  # apply(v_fd_alg_info, 1, function(x){cat(x); cat("\n")})
  compute_fd_middle <- function (query_start,res_middle) {

                 # nb of dash for left part

                 left_dashes_count = as.integer(query_start) - 1           # attention provient de l'appel  de la fct (semble globale)
                 #cat("\nquery_start ",query_start," res_middle ",res_middle,"\n")
                 st_dash <- str_dup("-", left_dashes_count)                # create a string of -
                 st_middle_and_left_dashes <- paste0(st_dash,res_middle)   # middle string is completed by left dashes

    st_middle_and_left_dashes
  }
  # appel à cette fct.
  # apply(df,1,compute_fd_middle(x))

  # adds a new column a reslides any gaps
  # use apply to speed  up for/next

  df$fd_middle <- NA  # used to avoid warning
  # df %>% add_column(fd_middle)
  for (i in 1:nrow(df)){
      st_middle = l_fd_alg_info[[i]]$middle_line %>% trimws()
      #cat("\n==> :",st_middle,":")
      #cat(compute_fd_middle(l_fd_alg_info[[i]]$query_start,st_middle))
      df$fd_middle[i] <- compute_fd_middle(l_fd_alg_info[[i]]$query_start,st_middle)
  }

  # inner function
  # turns into character and replace spaces with dashes
  fac_to_dashes <- function(x) x %>% as.character() %>% gsub(" ", "-", .)

### Commenté par fd car la col. peptide n'existe pas ds la df vide passee
###  # use it to polish peptide columns
###    df %<>% mutate(peptide = fac_to_dashes(peptide),
###                   middle  = fac_to_dashes(middle),
###                   blast   = fac_to_dashes(blast))
###
###    # adds a new column a reslides any gaps
###    df %<>% mutate(middle2=NA)
###    for (i in 1:nrow(df)){
###      df$middle2[i] <- best_sliding(df$blast[i], df$middle[i])
###    }

  ###  VERSION - FD : we add fd_column  query column left completed with dash
  # use it to polish peptide columns
    df %<>% mutate(peptide = fac_to_dashes(peptide), # enlevé pou test, car on a pas de colonne peptide entrée (vient du fasta)
                   middle  = fac_to_dashes(middle), # on remplace les espaces internes par des -
                   blast   = fac_to_dashes(blast)  ,
                   fd_middle=fac_to_dashes(fd_middle) # Ajout de cette col.
                   )
                   # cat("\n === ============================== == \n") # FD : à enlever
                   # cat("\n === fd_middle == \n")                      # FD : à enlever
                   # print(df$fd_middle)                                # FD : à enlever


    # adds a new column a reslides any gaps
    df %<>% mutate(middle2=NA)
    for (i in 1:nrow(df)){
           df$middle2[i] <- best_sliding(df$blast[i], df$middle[i])


    }

    # FD : new column fd_middle2
    # process each row fd_middle of df matrix
    # use best_sliding to complete fd_middle with left dashes if necessary
    # and fill fd_middle2 column with result
    df %<>% mutate(fd_middle2=NA) # add new column and store score coming from best sliding
    for (i in 1:nrow(df)){
           df$fd_middle2[i] <- best_sliding(df$blast[i], df$fd_middle[i])
           #print(df$fd_middle2[i])

    }
    # fin ajout FD

    #cat("middle ",df$middle) # ligne fd
    #cat("middle2 ",df$middle2) # ligne fd

  cat("* Best blasts successfully merged\n")
  return(df)
} # end of add_best
