require(tidyverse)
require(stringr)

# load class I results

setapred = read.csv("3_blast_mismatches_I.csv")
# quick peek, expect 20 columns
# Add infos from Propas:
# Please run:
# perl code/ProPAS.commandline.1.03.pl -i  1_blastme_I.fasta -of csv
system("perl code/ProPAS.commandline.1.03.pl -i  1_blastme_I.fasta -of csv")

propasdata = read.csv("1_blastme_I.fasta.propas.csv",skip=2)

# Bind propas and prediction data
setapred = cbind(setapred,propasdata)

# How many blast matches ? 
mismatch = setapred %>% select(c(1,seq(5,13))) %>% gather(pos,val,2:10) %>% group_by(peptide) %>% summarize(mismatch = length(grep("FALSE",val)))


# Add matches number to data
setapred = setapred %>% left_join(mismatch,by="peptide")



# Define function to get number of duplicated supertype:
get_duplicated_supertype = function(st){
	sst = strsplit(st,"\\+")[[1]]
	supertypes  = substring(sst,3)
	res = supertypes[duplicated(supertypes)]
	return(res) }

						       
setapred = setapred %>% rowwise %>% mutate(promiscuity = length(get_duplicated_supertype(as.character(supertype))))


# how many of predictor for the peptides  ? 

#setapred %>% rowwise %>% mutate(numpred = length(setdiff(str_split_fixed(predictor,"\\+",3),c("")))) %>% glimpse
setapred = setapred %>% rowwise %>% mutate(Num_predictor = length(setdiff(str_split_fixed(predictor,"\\+",3),c(""))))
setapred = setapred %>% rowwise %>% mutate(hasI = grepl("I",predictor)) 
setapred = setapred %>% rowwise %>% mutate(hasN = grepl("N",predictor)) 
setapred = setapred %>% rowwise %>% mutate(hasS = grepl("S",predictor)) 


# Score of N predictor (if any)
# home made function to extract scores from strings
get_scores = function(st){
	sst = strsplit(st,"\\_")[[1]]
	return(sst[3]) }

get_pred_scores = function(st,pred){
	sst = strsplit(st,"\\+")[[1]]
	test = sst[grep(pred,sst)]
	test2 = sapply(test,get_scores)
	return(as.numeric(as.vector(test2))) }


setapred = setapred %>% rowwise %>% mutate(scoreN = max(get_pred_scores(as.character(score),"N")))
setapred$scoreN = ifelse(setapred$scoreN == -Inf,NA,setapred$scoreN)


# Now for the filtering step:
#setapred %>% arrange(desc(Num_predictor),desc(mismatch),desc(scoreN)) %>% glimpse

# convert supertype to vector of characters 
setapred$supertype = as.vector(setapred$supertype)

# create new dataframe to duplicate results by supertypes
setapred2 = data.frame()
HLA_restriction = c()
for(row in 1:nrow(setapred))
    {
#   print(setapred[row,])
    st = as.vector((get_duplicated_supertype(as.character(setapred[row,"supertype"]))))
    if(length(st)==0){ 
	setapred2 = rbind(setapred2,setapred[row,])
    	HLA_restriction = append(HLA_restriction,NA)}
    for(s in st){
    	setapred2 = rbind(setapred2,setapred[row,])
    	HLA_restriction = append(HLA_restriction,s)
    }
    }


setapred2$HLA_restriction = HLA_restriction 
allele2supertype = read.table("table/map_supertypes_alleles.csv",header =T,sep=';')

str(allele2supertype)
allele2supertype %>% select(allele)

allele2supertype = allele2supertype %>% rowwise %>% mutate(allele_simple =  gsub('\\*', '', allele))


# Score of N predictor (if any)
# home made function to extract scores from strings
get_allele_scores = function(st){
	sst = strsplit(st,"\\_")[[1]]
	return(sst[2:3]) }

get_pred_allele_scores = function(st,pred){
	sst = strsplit(st,"\\+")[[1]]
	test = sst[grep(pred,sst)]
	test2 = sapply(test,get_allele_scores)
	return(test2) }

setapred2$score = as.vector(setapred2$score)

allele2supertype$allele_simple = as.vector(allele2supertype$allele_simple)
allele2supertype$supertype = as.vector(allele2supertype$supertype)

scorevals = c()
for(row in 1:nrow(setapred2)){
	score = as.character(setapred2[row,"score"])
	res = get_pred_allele_scores(score,"N")
	if(length(res)>0 & setapred2[row,"promiscuity"]>0){
	   for(i in 1:dim(res)[2]){
	   	allele = as.character(res[1,i])
		scoreN = as.numeric(res[2,i])
		super = as.character(allele2supertype[allele2supertype$allele_simple==allele,"supertype"])
		if(setapred2[row,"HLA_restriction"]==super){scorevals=append(scorevals,scoreN)}
	   }
	}else{scorevals=append(scorevals,NA)}
}

setapred2$scoreN = scorevals

# Now look at supertype specific anchor positions in blast results ...
supertype = read.csv("table/anchorpositions.csv",sep=",")

anchormm = c()
for(i in 1:nrow(setapred2)){
	if(!is.na(setapred2[i,"HLA_restriction"])){
		pos = supertype %>% filter(supertype == setapred2[i,]$HLA_restriction) %>% mutate(pos = anchor + 4) %>% select(pos)
		anchormm = append(anchormm,sum(as.integer(!setapred2[i,as.vector(pos$pos)])))
		}
	else{
		anchormm = append(anchormm,NA)
		}
	}
setapred2$anchormm = anchormm


write.table(file="res_classI.csv",setapred2,sep=";",row.names=F)


