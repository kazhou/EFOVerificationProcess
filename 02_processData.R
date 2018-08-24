# Add EFO IDs to Therapeutic Targets Database & GWAS catalog downloads

library(rJava)
library(rols)
library(jsonlite)
library(rjson)
library(data.table)
library(dplyr)
library(tidyr)
library(biomaRt)

load("~/src/ttd_df.rds")  # orginal TTD download
load("~/src/gwascat.rds")   #original GWAS catalog
gwascat<-gwascat[gwascat$P.VALUE<=5E-8,]
disease.list<-as.character(unique(sort(gwascat_ext$MAPPED_TRAIT)))

# EFO
load("~/src/efo_terms.rds") #created 7/9/18; update by reassigning and saving efo_terms<-rols::terms(x="efo")
load("~/src/efo_df.rds") # efo_df<-as(efo_terms, "data.frame")


TTDtoEFO <- function(df=ttd_df) {
  # Converts TTD disease traits to EFO ID through name, ICD, and/or keyword query; call if you need to regenerate ttd_df_ext
  # 
  # Args: 
  #   df: table with target-indication data from TTD; post web-scraping
  # Returns:
  #   TTD data table with EFO IDs column
  temp <- unique(dplyr::select(ttd_df,Disease, ICD9, ICD10))
  dis_efo <- as.data.frame(cbind(temp, matrix(nrow=length(temp$Disease),ncol=2)))
  rownames(dis_efo) <- NULL
  colnames(dis_efo) <- c("Disease", "ICD9","ICD10","EFO", "EFO_Src")
  for (i in 1:length(dis_efo$Disease)) {
    names <- unlist(strsplit(as.character(dis_efo$Disease[i]), "; "))  # for multiple diseases in same line
    efos <- c()
    efo.src <- c()
    for(x in 1:length(names)) {
      name <- names[x]
      icd9 <- unlist(strsplit(as.character(dis_efo$ICD9[i]),", "))
      if(isEmpty(icd9)){icd9<-"XXXXXX"}  # no matches
      icd10 <- unlist(strsplit(as.character(dis_efo$ICD10[i]),", "))
      if(isEmpty(icd10)){icd10<-"XXXXXX"}  # no matches
      res <- c(paste('ICD9:',icd9,sep=""), (paste('ICD10:',icd10,sep="")))
      
      if (!isEmpty(efo_df[grep(paste("^",name,"$",sep=""), efo_df$label, ignore.case = TRUE),])) {  
        #Check exact name match
        efos <- c(efos, efo_df$id[grep(paste("^",name,"$",sep=""), efo_df$label, ignore.case = TRUE)])
        efo.src <-c (efo.src, "Name match - exact")
      } else if (!isEmpty(efo_df[grep(paste("^",name,"$",sep=""), efo_df$first_synonym, ignore.case = TRUE),])){
        #Check exact synonym match
        efos <- c(efos, efo_df$id[grep(paste("^",name,"$",sep=""), efo_df$first_synonym, ignore.case = TRUE)])
        efo.src <- c(efo.src,"Synonym match - exact")
      }  else if (!isEmpty(res)) {
        #Query ICD
        j <- 1
        while (j <= length(res)) {
          qry <- OlsSearch(q=res[j], ontology="EFO")
          if(qry@numFound > 0)  {
            qry <- olsSearch(qry)
            qtrms <- as(qry, "data.frame")
            efos <- c(efos,unlist(lapply(keywords, function(k) {qtrms$obo_id[grep(k, qtrms$label, ignore.case = TRUE)]})))
            rm(qtrms)
            efo.src <- c(efo.src,rep("Query match - ICD", length(efos)))
            #case, british v american spelling, typos
          }
          rm(qry)
          if (length(efo.src) > 0) {
            j <- length(res)+1
            # End loop
          }
        }
      } else {
        # Query keywords
        keywords <- unlist(strsplit(name," "))
        qry <- OlsSearch(q=name, ontology="EFO")
        if (qry@numFound > 0)  {
          qry <- olsSearch(qry)
          qtrms <- as(qry, "data.frame")
          efos <- c(efos,unlist(lapply(keywords, function(k) {qtrms$obo_id[grep(k, qtrms$label, ignore.case = TRUE)]})))
          rm(qtrms)
          efo.src <- c(efo.src,rep("Query match - keywords", length(efos)))
        }
        rm(qry)
      }
      
      if (!isEmpty(efos)){
        dis_efo$EFO[i]<-paste(efos, collapse=", "); dis_efo$EFO_Src[i]<-paste(efo.src,collapse=", ")
      } else {
        next
      } 
    }
  }
  
  ttd_df_ext<-left_join(ttd_df, dis_efo, by=c("Disease"="Disease", "ICD9"="ICD9", "ICD10"="ICD10"))
  # save(ttd_df_ext, file="~/../e0366873/work/ANDIE/src/ttd_df_ext.rds")
  return(ttd_df_ext)
}

GWAStoEFO <- function(gwascat, disease.list){
  # Adds EFO ID column to GWAS catalog; call if you need to regenerate gwascat.ext
  # 
  # Args: 
  #   gwascat: original GWAS catalog
  #   disease.list: unique diseases in GWAS catalog
  # Returns:
  #   GWAS catalog with column for EFO IDs
  dis.extra.list <- as.data.frame(cbind(disease.list, matrix(nrow=length(disease.list), ncol=2)), stringsAsFactors=FALSE)
  colnames(dis.extra.list) <- c("Disease", "EFO", "EFO_Src")
  for (i in 1:length(disease.list)) {
    name <- disease.list[i]
    if (!isEmpty(grep(',',name))) {  # split disease listed in same line
      names <- unlist(strsplit(as.character(name), ", "))
    } else {
      names <- name
    }
    
    efoIDs <- c()
    efo.src <- c()
    for (x in 1:length(names)){
      name <- names[x]
      efoID <- efo_df$id[efo_df$label==name]
      if (length(efoID) > 1) {
        efoID <- paste(efoID, collapse=", ")
      }
      efoIDs <- c(efoIDs, efoID)
      if (!isEmpty(efoIDs)) {
        efo.src <- c(efo.src, "Name match - exact")
      } else {
        qry <- OlsSearch(q=name, ontology="EFO")
        if (qry@numFound > 0)  {
          qry <- olsSearch(qry)
          qtrms <- as(qry, "data.frame")
          efoIDs <- c(efoIDs,unlist(lapply(keywords, function(k) {qtrms$obo_id[grep(k, qtrms$label, ignore.case = TRUE)]})))
          rm(qtrms)
          efo.src <- c(efo.src,"Query match - keywords")
          #case, british v american spelling, typos
        }
        rm(qry)
      }
    }
    if (isEmpty(efoIDs)) {
      dis.extra.list$EFO[i] <- NA
      dis.extra.list$EFO_Src[i] <- NA
      next
    }
    dis.extra.list$EFO[i] <- paste(efoIDs, collapse=", ")
    dis.extra.list$EFO_Src[i] <- paste(efo.src, collapse=", ")
  }
  #  save(dis.extra.list,file="~/src/dis_extra_list.rds")
  
  # gwascat.ext<-left_join(gwascat,dis.extra.list,by=c("MAPPED_TRAIT"="Disease"))
  # save(gwascat.ext, file="~/src/gwascat_ext.rds")
  return(gwascat.ext)
}

load("~/src/ttd_df_ext.rds")
load("~/src/gwascat_ext.rds")

