require(tidyverse)
require(stringr)

# load class II results


setapred = read.csv("3_blast_mismatches_II.csv")
# quick peek, expect 20 columns
# Add infos from Propas:
# First extract peptides sequences:

file.create("peptideII.fasta")
for(row in 1:nrow(setapred)){
	write(paste(">",row,sep=""),file = "peptideII.fasta",append=T)
	write(as.character(setapred[row,"full_peptide"]),file = "peptideII.fasta",append=T)
}

# perl code/ProPAS.commandline.1.03.pl -i  peptideII.fasta  -of csv
system("perl code/ProPAS.commandline.1.03.pl -i  peptideII.fasta -of csv")

propasdata = read.csv("peptideII.fasta.propas.csv",skip=2)

# Bind propas and prediction data
setapred = cbind(setapred,propasdata)

# How many blast matches ? 

mismatch = c()
for(row in 1:nrow(setapred)){
	mismatch = append(mismatch,sum(as.integer(!as.vector(setapred[row,6:14]))))
}
# Add matches number to data
setapred$mismatch = mismatch


get_duplicated_alleles = function(st){
	tmp = gsub('\\*', '', st) %>%  gsub("_","",.)  %>% gsub(':','',.) %>% gsub('HLA-','',.) %>% gsub('-','/',.)
	tmp = strsplit(tmp,"\\+")[[1]]
	tmp = tmp[duplicated(tmp)]
	return(tmp)
}


setapred =  setapred %>% rowwise %>%  mutate(promiscuity = length(get_duplicated_alleles(as.character(allele))))

# how many of predictor for the peptides  ? 

#setapred %>% rowwise %>% mutate(numpred = length(setdiff(str_split_fixed(predictor,"\\+",3),c("")))) %>% glimpse
setapred = setapred %>% rowwise %>% mutate(Num_predictor = length(setdiff(str_split_fixed(predictor,"\\+",3),c(""))))
setapred = setapred %>% rowwise %>% mutate(hasI = grepl("I",predictor)) 
setapred = setapred %>% rowwise %>% mutate(hasN = grepl("N",predictor)) 


# Now for the filtering step:
#setapred %>% arrange(desc(Num_predictor),desc(mismatch),desc(scoreN)) %>% glimpse

# convert supertype to vector of characters 
setapred$allele = as.vector(setapred$allele)

# create new dataframe to duplicate results by supertypes

setapred2 = data.frame()
HLA_restriction = c()
setapred$allele = as.vector(setapred$allele)
for(row in 1:nrow(setapred))
    {
    st = as.vector((get_duplicated_alleles(as.character(setapred[row,"allele"]))))
    if(length(st)==0){ 
	setapred2 = rbind(setapred2,setapred[row,])
    	HLA_restriction = append(HLA_restriction,NA)}
    for(s in st){
    	setapred2 = rbind(setapred2,setapred[row,])
    	HLA_restriction = append(HLA_restriction,s)
    }
}
setapred2$HLA_restriction = HLA_restriction 


# Score of N predictor (if any)
# home made function to extract scores from strings
get_allele_scores = function(st){
	sst = strsplit(st,"\\_")[[1]]
	score = sst[length(sst)]
	tmp = sst[-length(sst)]
	tmp = tmp[-1]
	allele  = paste(tmp,collapse ="_")
	return(c(allele,score))}

get_pred_allele_scores = function(st,pred){
	sst = strsplit(st,"\\+")[[1]]
	test = sst[grep(pred,sst)]
	test2 = sapply(test,get_allele_scores)
	return(test2) }


setapred2$score = as.vector(setapred2$score)
scorevals = c()
for(row in 1:nrow(setapred2)){
	score = as.character(setapred2[row,"score"])
	res = get_pred_allele_scores(score,"N")
	if(length(res)>0 & setapred2[row,"promiscuity"]>0){
	   	found = FALSE
		for(i in 1:dim(res)[2]){
	   	allele = as.character(res[1,i])
		scoreN = as.numeric(res[2,i])
		allele = gsub('\\*', '', allele) %>%  gsub("_","",.)  %>% gsub(':','',.) %>% gsub('HLA-','',.) %>% gsub('-','/',.)
			if(setapred2[row,"HLA_restriction"]==allele){
			found = TRUE
			scorevals=append(scorevals,scoreN)
			}
	   	}
	   	if(!found){scorevals =append(scorevals,NA)}
	}else{
		scorevals=append(scorevals,NA)}
}
setapred2$scoreN = scorevals

write.table(file="res_classII.csv",setapred2,sep=";",row.names=F)







