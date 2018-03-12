#STA 525

#####################
##  BiomaRt_intro  ##
#####################

## Load these from your directory!
load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Array.RData") #ArrayData
load("~/Google Drive/Stats for Bioinformatics/HW3/GeneMap.RData")    #GeneMap
load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Valid.RData") #ValidData



# Biomart - Database
# ========================================================
#   - A bioinformatics hub, where information from various databases on a certain gene, chromosome location, identifiers for different high throughput platforms etc. can be acquired
# 
# - **Query pipeline:**      
#   Dataset -> Filters -> Attributes -> Results
# 
# [Biomart Database](http://useast.ensembl.org/biomart/martview/7b7bc13e12cb5ea67f1438dde0fe05fa)


# biomaRt - R interface/package
# ========================================================
#   
#   * R package, providing a R interface that alows us to get data from Biomart data base without going to their web page.          
# 
# * It speeds up the process as the acquired data directly uploded to our R environment.
# 
# biomaRt - R interface/package
# ========================================================
#   

## 1. Choose your mart ##
#########################

# Install biomaRt if you haven't done so
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

library("biomaRt")
listEnsembl()
# or another way to go: 
listMarts()

#Select Ensembl Genes mart
ensembl=useMart("ensembl")

 
# ## 2. Choose your dataset ##
# ############################
# * Each mart has different datasets.     
# 
# * For example Ensembl Genes have a different dataset for every species in their catalog.        

listDatasets(ensembl)[c(1:2, 32),1:2]

#our samples are from humans!
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


## 3. How to build a biomaRt query  ##       
######################################  
#  Dataset -> Filters -> Attributes -> Results                

# The `getBM()` function is the main query function in biomaRt.               
# **It has four main arguments:**    
#   
#   1.  **attributes**: a vector of attributes to retrieve (the output of the query).  
# 2.  **filters**: a vector of filters that one will use as input to the query.  
# 3.  **values**: a vector of values for the filters (for multiple filters, a list)  
# 4.  **mart**: an object of class Mart, which is created by the useMart function.  

## Let's list the first 3 attributes and filters for our `ensembl` mart
listAttributes(ensembl)[1:3,]
listFilters(ensembl)[1:3,]


# biomaRt helper functions
# ========================================================
#   
#   * "`listAttributes()`"  gives a large list of attributes, can be very hard to navigate        
# 
# * We can list the filters but it is not always straight forward to figure out what values that specific filter contains.      
# 
# * Some filters that start with "with" require a boolean as a value rather than a character or numceric vector.   

### 1. Attribute Pages ###
###########################

attributePages(ensembl)[1:3]
listAttributes(ensembl, page="feature_page")[1,]
listAttributes(ensembl, page="sequences")[1,]

## 2. Determine filter type ##
##############################

filterType("phenotype_description",ensembl)
filterType("with_wikigene",ensembl)

## 3. Determine acceptable variables for the filter ###
########################################################

# this gives a very long string, separated with commas!
PhenOpt <- filterOptions("phenotype_description",ensembl) 
# so we split them
PhenOpt <- unlist(strsplit(PhenOpt, ","))
# We have a lot of variables to choose from! 
length(PhenOpt)
PhenOpt[c(64, 1989,1990)] # Random examples


## 4. Grepping the variables we want  ##   
########################################  
#  Since the list of accepted variables is very long, 
#grepping a word (or beginning of a word) and looking for the values we want is the way to go! 
  
# get phenotype variables related to 
# "inflammation" or "inflammatory""
infl <- grep("inflammat", PhenOpt)

#the list is shorter now, let's see couple of of them 
PhenOpt[infl[c(7,16, 18)]]


# biomaRt - An Example 
# ========================================================
  
#  Retrieve all EntrezGene IDs and HUGO gene symbols of genes which have a ”MAP kinase activity” GO term associated with it.
# 
# * GO identifier for MAP kinase activity is [GO:0004707](http://amigo.geneontology.org/amigo/term/GO:0004707).       
# * We will set `filter=go` and  `attributes=c(entrezgene, hgnc_symbol)`.      


MAP_genes <- getBM( attributes = c('entrezgene','hgnc_symbol'), 
                    filters='go_id', 
                    values='GO:0004707', mart=ensembl)


# **Our Results with the MAP kinase activity pathway**

MAP_genes[1:7,]


# An Example from Homework 3
# ========================================================  
#   
#   * We got the entrezgene IDs from the given Berry data            
# 
# * Now we want to see which ones are actually involved in inflammatory response,     since we would expect tuberculosis to activate some sort of inflammatory response pathways and therefore the genes involved in it.     
# 
# * We also want to see the actual gene symbols, as entrez IDs do not tell much to biologists.       
# 
# * Lastly, we want to know the GO IDs and GO names of these genes, to see what GO terms are associated with these genes.


## Remember the Query pipeline!!:     
  
# 1.  Dataset = ensembl          
# 2.  Filters = "entrez ID" and "phenotype description"               
# 3.  Values = "our genes' entrez IDs", and phenotype terms associated with inflammation as we found in the previous slide                   
# 4.  Attributes = "HGNC Gene symbol", "entrez IDs" , "GO ID" and "GO name"      
# 5. ------>  Results !!
  
### 1.  Let's load our Entrez IDs from the Berry data

# load the data from your directory! When you run the code in the GSEA part 
# provided by Dr. Gaile, you'll get these IDs

load("./myEntrezs.RData")
head(myEntrezs)


### 2.  Let's get the phenotype description variables we got in the previous slide

infl.Val <- PhenOpt[infl[c(7,16, 18)]]
infl.Val


### 3.  Finally results!          
# See how the values are assigned! We put the vectors containing our desired filter values 
#  in list in the same order with our filters.       


myresult<-getBM(attributes=c("entrezgene","hgnc_symbol", 
                             "go_id", "name_1006"), 
                filters= c("entrezgene", 
                           "phenotype_description"), 
                values= list(myEntrezs, infl.Val), mart=ensembl)

### Now let's look at our results!        
  
# 4th column, "name_1006" is our GO names

myresult[c(1,90,130),]

