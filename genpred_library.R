#' Anyela Camargo
#' 
#' 
#' 
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library('qtl2')
library('devtools')
library('data.table')
library('dplyr')
library('doParallel')


# Select RILS's marker profile
#' @param filename
#' @param markername
searchMarker <- function(sdat, markername){
  
  
  i = match(markername, colnames(sdat))
  m = sdat[, c(1,i), with=FALSE]
  return(m)
  
}


#' @param RILs profile for a given marker
convert_snp <- function(x) {
  #print(x)
  #convert to {-1,0,1,NA}
  alleles <- c('0', '2', '1'); # 0=AA, 2=BB
  y <- rep(NA,length(x))
  #print(alleles);
  y[which(x==alleles[1])] <- '1'
  y[which(x==alleles[2])] <- '3'
  y[which(x==alleles[3])] <- '2'
  #break();
  return(y)
}


#' process genodata
#' @param dataril, genetic map
#' @param datamap, map
get_rils <- function(sdat)
{
  
  clr = ncol(sdat)
  rr = nrow(sdat)# Cols dataframe
  
  # Convert snps to 1/0
  R <- apply(sdat[1:rr,13:clr],1, convert_snp)
  
  markername = as.vector(sdat[1:rr, 'Name'])
  genotypename = colnames(sdat)[13:clr]
  ril_data <- data.table(genotypename, R)
  colnames(ril_data) <- c('id', t(markername))
  
  fwrite(ril_data, file = 'MAGIC/magic_geno.csv', col.names=T,row.names=F,
         quote = F,sep=",")
  return(ril_data)
  
}

#' @param sdat raw data
get_parents <- function(sdat){
  
  colname = c("Name", "ALCHEMY",  "BROMPTON", "CLAIRE", "HEREWARD", "RIALTO", "ROBIGUS",
              "SOISSONS", "XI-19")  
  data_parent <- sdat[, colname, with=FALSE]
  
  clr = ncol(data_parent)
  rr = nrow(data_parent)# Cols dataframe
  markername = as.vector(data_parent[['Name']])
  
  R <- apply(data_parent[1:rr,2:clr],2, convert_snp)
  parent_data <- data.frame(markername, R)
  colnames(parent_data) <- colnames(data_parent)
  setDT(parent_data)
  fwrite(parent_data, file = 'MAGIC/magic_foundergeno.csv', col.names=T,row.names=F,
         quote = F,sep=",")
  
  
}

#' @param sdat raw data
get_map <- function(sdat){
  
  
  cname <- c('Name', 'Chr', 'Pos')
  map_data <- sdat[, cname, with= FALSE] 
  setDT(map_data)
  fwrite(map_data, file = 'MAGIC/magic_gmap.csv', col.names=T,row.names=F,
         quote = F,sep=",")
  
}

#' @param sdat raw data
get_pheno <- function(sdat){
  
  
  pheno_data = sdat[, lapply(.SD, mean, na.rm=TRUE), 
           by=geno, .SDcols=c(colnames(sdat)[3:20]) ]
  colnames(pheno_data)[1]  = 'ind'
  
  fwrite(pheno_data, file = 'MAGIC/magic_pheno.csv', col.names=T, row.names=F,
         quote = F,sep=",")
}


#' Get all the data ready for mapping
get_data <- function(){
  
  
  raw_geno <- fread('NIAB_MAGIC_ELITE_genotypes_mapped_markers.csv')
  raw_pheno <- fread('magic_gt.csv', header = TRUE)
  get_parents(raw_data)
  get_rils(raw_data)
  get_map(raw_data)
  get_pheno(raw_pheno)
  
  
}


marker_presence_abscence <- function(){
  ##process rils and produce regression plots
  
  fname <- 'MAGIC/magic_geno.csv'
  rils <- fread(fname, sep =',', header=TRUE)
  fname <- 'MAGIC/magic_pheno.csv' 
  pheno <- fread(fname, sep=',', header=TRUE)
    markername <- c('RAC875_rep_c105718_585', 'RHT2')
  m <- searchMarker(rils, markername)
  m <- merge(m, pheno, by.x= 'id', by.y='ind')
  
  cov_file <- read.table('MAGIC/magic_covar.csv', header=TRUE, sep=',')
  cov_file <- merge(cov_file, m[,1:3], by.x= 'ind', by.y='id')
  fwrite(cov_file, file='magic_covar1.csv', sep=',')
  
}


