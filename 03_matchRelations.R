# Script to compare EFO IDs of predicted and expected indications

library(dplyr)
library(data.table)
library(rJava)
library(rols)
library(jsonlite)
library(openxlsx)
library(S4Vectors)
library(tidyr)
library(doParallel)
library(foreach)
cl<-makeCluster(detectCores()-1) 
registerDoParallel(cl)

#bash shell arguments
args<-commandArgs(TRUE)

# test if there is at least one argument: if not, return an error
if (!grepl("all|approved|^[:upper:]*([,]?[:upper:]*)*", args[1])) { 
  stop("Targets format incorrect", call.=FALSE)
} 
if (!grepl("[0-9]*", args[2])) {
  stop("Score threshold format incorrect", call.=FALSE)
}
if (!grepl("^[0-9]*E?-?[0-9]*$", args[3])) {
  stop("P-value format incorrect", call.=FALSE)
}
if (!grepl("normal|FDR*", args[4])) {
  stop("P-value option format incorrect", call.=FALSE)
}
if (!grepl("all|[0-9]*", args[5])) {
  stop("Rows format incorrect", call.=FALSE)
}

input <- as.character(args[1])
threshold <- as.integer(args[2])
p.val <- as.numeric(args[3])
pstat <- as.character(args[4])
nrows <- as.character(args[5])
date.1 <- as.character(args[6])

load("~/src/ttd_df_ext.rds")
load("~/efo_terms.rds")
load("~/efo_df.rds")

date <- gsub("-","", Sys.Date())
files <- list.files(path=paste("~/files/",date.1,"_FET/txt/",sep=""), pattern = paste(threshold,"_[0-9]*.txt", sep=""))
completed <- unlist(lapply(as.list(files), function(f) {substr(f, start=1, stop=as.integer(unlist(gregexpr(pattern="_",f))[1])-1)}))

term.IDs <- function(terms) {
  # Get EFO IDs of a Terms object
  #
  # Args:
  #   terms: a Terms object (rols package)
  #
  # Returns:
  #   EFO IDs of terms
  p1<-c()
  for(l in 1:length(terms)) {
    if(is.null(terms[[l]])) {next}
    p1<-c(p1, names(terms[[l]]@x))
  }
  return(p1)
}

GetSibIDs <- function(terms) {
  # Extract EFO IDs of siblings of a Terms object
  #
  # Args:
  #   terms: a Terms object (rols package)
  #
  # Returns:
  #   EFO IDs of sibling terms
  s1 <- c()
  for (l in 1:length(terms)) {
    if (is.null(terms[[l]])) {
      next
    }
    s1 <- c(s1, names(terms[[l]]@x))
  }
  return(s1)
}

GetSibTerms <- function(term) {
  # Generates list of Terms representing siblings of a EFO ID
  #
  # Args:
  #   term: the EFO ID to find siblings of 
  #
  # Returns:
  #   List of Terms containing siblings
  return(lapply(names((parents(efo_terms@x[[term]]))@x), function(t) {if(!is.null(efo_terms@x[[t]])) {children(efo_terms@x[[t]])}})) #if(inherits(t, "Term")) {
}

MatchRelated <- function(file, p.val, nrows, status) {
  # Generates table of matches between predictd and expected indications based on p-value cutoff
  #
  # Args:
  #   file: filename of predicted indications list of target
  #   p.val: p-value cutoff
  #   nrows: "all" or number of rows to consider
  #
  # Returns:
  #   Table of matches
  
  target <- substr(file, start=1, stop=as.integer(unlist(gregexpr(pattern="_",file))[1])-1)
  dat <- read.delim(file=paste("~/files/",date.1,"_FET/txt/",file,sep=""), header=T, sep="\t")
  
  #sort by p-value
  dat2 <- dat[dat$P_Value<=p.val,] #accept user input
  dat2 <- dplyr::arrange(dat2, P_Value)
  
  if (nrows=="all") {
    top_dis<-dat2
  } else {
    top_dis<-dplyr::slice(dat2, 1:min(as.integer(nrows), length(dat2$Disease_Trait)))
  }
  
  n.pred<-length(top_dis$Disease_Trait)
  
  #"expected"
  ttd <- arrange(ttd_df_ext[ttd_df_ext$HGNC_Symbol==target & !is.na(ttd_df_ext$HGNC_Symbol),], Status)
  if(status=="approved") {
    ttd.appr <- filter(ttd, Status=="Approved") #FDA approved targets
  } else if(status=="phase3") {
    ttd.appr <- filter(ttd, Status %in% c("Phase 3","Phase 4"))
  } else {
    ttd.appr <- filter(ttd, Status %in% c("Phase 3","Phase 4", "Approved"))
  }
  
  # Predicted EFOs
  efo.res <- unique(unlist(strsplit(as.character(top_dis$EFO),", ")))
  if (is.null(efo.res)|| is.na(efo.res)) {
    conf <- as.data.frame(cbind(target, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected", "Status", "EFO_Src")
    return(conf)
  }
  
  # Expected EFOs and their source
  efo.exp.df <- as.data.frame(unlist(lapply(ttd.appr$EFO, 
        function(e) {paste(unlist(strsplit(as.character(e), ",")), as.character(ttd.appr$EFO_Src[grep(e, ttd.appr$EFO)]), sep = "/")})))
  colnames(efo.exp.df) <- c("IDSRC")
  efo.exp.df <- separate(efo.exp.df, IDSRC, sep = "/", into = c("EFO", "EFO_Src"))
  
  # Unique, non-NA
  efo.exp <- unique(as.character(efo.exp.df$EFO))
  efo.res <- efo.res[!is.na(efo.res)]
  efo.exp <- efo.exp[!is.na(efo.exp)]
  
  relations.df <- as.data.frame(cbind(efo.res,matrix(nrow=length(efo.res),ncol=3)), stringsAsFactors=FALSE)
  colnames(relations.df) <- c("Predicted_Indication", "Parents", "Siblings", "Children")
  
  for (i in 1:length(efo.res)) {
    e <- efo.res[i]
    if (!is.null(parents(efo_terms@x[[e]]))) {
      relations.df$Parents[i] <- paste(unique(names(parents(efo_terms@x[[e]])@x)), collapse=", ") 
      
      sibs <- unique(GetSibIDs(GetSibTerms(e)))
      sibs <- sibs[sibs!=e]
      if (is.null(sibs)) { 
        relations.df$Siblings[i] <- NA
      } else {
        relations.df$Siblings[i] <- paste(sibs, collapse=", ")
      }
    }
    
    if (!is.null(children(efo_terms@x[[e]]))) { 
      relations.df$Children[i] <- paste(unique(names(children(efo_terms@x[[e]])@x)),collapse=", ")
    }
  }
  
  all.results <- union(union(relations.df$Predicted_Indication, unlist(strsplit(as.character(relations.df$Siblings),", "))),
                     union(unlist(strsplit(as.character(relations.df$Parents),", ")), unlist(strsplit(as.character(relations.df$Children),", "))))
  
  matches <- unique(intersect(all.results,efo.exp)) 
  total <- length(ttd.appr$Disease)
  match.rate <- as.character(length(matches)/total)
  match.ratio <- paste(length(matches),":",total,sep="")
  
  if (length(matches) == 0) {
    conf <- as.data.frame(cbind(target, match.rate, match.ratio, n.pred, NA, NA, NA, NA, NA, NA, NA))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected","Status", "EFO_Src")
    return(conf)
  }
  
  match.dis <- unlist(lapply(matches, function(m){efo_df$label[efo_df$id==m&!is.na(efo_df$id)]}))
  
  match.df <- as.data.frame(matrix(nrow=0,ncol=11))
  colnames(match.df) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected","Status", "EFO_Src")
  for (m in 1:length(matches)) {
    match <- matches[m]
    
    ttd.dis <- as.character(ttd.appr$Disease[grep(match,ttd.appr$EFO)])
    if (length(ttd.dis[grep(match.dis[m], unlist(strsplit(as.character(ttd.dis),'; ')), ignore.case=TRUE)])>0) {
      ttd.dis<-unique(grep(match.dis[m], unlist(strsplit(as.character(ttd.dis),'; ')), ignore.case=TRUE, value=TRUE))
      ttd.dis<-ttd.dis[!is.na(ttd.dis)]
    }
    
    relation.predicted <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
    colnames(relation.predicted) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
    if (TRUE %in% (relations.df$Predicted_Indication %in% match)) {
      res.dis <- c(as.character(top_dis$Disease_Trait[grep(match,top_dis$EFO)])) 
      rrow <- matrix(c(as.character("Exact"), paste(res.dis, sep="/ ")), nrow=1,ncol=2)
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Parents),", ")))
    if (TRUE %in% (tmp %in% match)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Parents)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Parent",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Siblings),", ")))
    if (TRUE %in% grepl(match, relations.df$Siblings)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Siblings)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Sibling",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Children),", ")))
    if (TRUE %in% grepl(match, relations.df$Children)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Children)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Child",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    }
    
    stat <- unique(ttd.appr$Status[ttd.appr$Disease==ttd.dis])
    src <- unique(unlist(strsplit(as.character(efo.exp.df$EFO_Src[efo.exp.df$EFO==match]),", ")))  
    
    conf <- cbind(target, match.rate, match.ratio, n.pred, match, match.dis[m], relation.predicted,paste(ttd.dis,collapse="/"), paste(stat, collapse="/"),
                  paste(src,collapse="/"))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected", "Status","EFO_Src")

    match.df<-unique(as.data.frame(rbind(match.df, conf)))
    
  }
  
  return(match.df)
}

MatchRelatedADJ <- function(file, p.val, nrows, status) {
  # Generates table of matches between predictd and expected indications based on FDR-adjusted p-value cutoff
  #
  # Args:
  #   file: filename of predicted indications list of target
  #   p.val: FDR p-value cutoff
  #
  # Returns:
  #   Table of matches
  target <- substr(file, start=1, stop=as.integer(unlist(gregexpr(pattern="_",file))[1])-1)
  dat <- read.delim(file=paste("~/files/",date.1,"_FET/txt/",file,sep=""), header=T, sep="\t")
  
  #sort by p-value
  dat2 <- dat[dat$P_Adjusted_FDR<=p.val,] #accept user input
  dat2 <- dplyr::arrange(dat2, P_Value)
  
  if (nrows=="all") {
    top_dis<-dat2
  } else {
    top_dis<-dplyr::slice(dat2, 1:min(as.integer(nrows), length(dat2$Disease_Trait)))
  }
  
  n.pred<-length(top_dis$Disease_Trait)
  
  #"expected"
  ttd <- arrange(ttd_df_ext[ttd_df_ext$HGNC_Symbol==target & !is.na(ttd_df_ext$HGNC_Symbol),], Status)
  if(status=="approved") {
    ttd.appr <- filter(ttd, Status=="Approved")  #FDA approved targets
  } else if(status=="phase3") {
    ttd.appr <- filter(ttd, Status %in% c("Phase 3","Phase 4"))
  } else {
    ttd.appr <- filter(ttd, Status %in% c("Phase 3","Phase 4", "Approved"))
  }
  
  # Predicted EFOs
  efo.res <- unique(unlist(strsplit(as.character(top_dis$EFO),", ")))
  if (is.null(efo.res)|| is.na(efo.res)) {
    conf <- as.data.frame(cbind(target, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected", "Status", "EFO_Src")
    return(conf)
  }
  
  # Expected EFOs and their source
  efo.exp.df <- as.data.frame(unlist(lapply(ttd.appr$EFO, 
           function(e) {paste(unlist(strsplit(as.character(e), ",")), as.character(ttd.appr$EFO_Src[grep(e, ttd.appr$EFO)]), sep = "/")})))
  colnames(efo.exp.df) <- c("IDSRC")
  efo.exp.df <- separate(efo.exp.df, IDSRC, sep = "/", into = c("EFO", "EFO_Src"))
  
  # Unique, non-NA
  efo.exp <- unique(as.character(efo.exp.df$EFO))
  efo.res <- efo.res[!is.na(efo.res)]
  efo.exp <- efo.exp[!is.na(efo.exp)]
  
  relations.df <- as.data.frame(cbind(efo.res,matrix(nrow=length(efo.res),ncol=3)), stringsAsFactors=FALSE)
  colnames(relations.df) <- c("Predicted_Indication", "Parents", "Siblings", "Children")
  
  for (i in 1:length(efo.res)) {
    e <- efo.res[i]
    if (!is.null(parents(efo_terms@x[[e]]))) {
      relations.df$Parents[i] <- paste(unique(names(parents(efo_terms@x[[e]])@x)), collapse=", ") 
      
      sibs <- unique(GetSibIDs(GetSibTerms(e)))
      sibs <- sibs[sibs!=e]
      if (is.null(sibs)) { 
        relations.df$Siblings[i] <- NA
      } else {
        relations.df$Siblings[i] <- paste(sibs, collapse=", ")
      }
    }
    
    if (!is.null(children(efo_terms@x[[e]]))) { #is there a hasCHildren funct
      relations.df$Children[i] <- paste(unique(names(children(efo_terms@x[[e]])@x)),collapse=", ")
    }
  }
  
  all.results <- union(union(relations.df$Predicted_Indication, unlist(strsplit(as.character(relations.df$Siblings),", "))),
                       union(unlist(strsplit(as.character(relations.df$Parents),", ")), unlist(strsplit(as.character(relations.df$Children),", "))))
  
  matches <- unique(intersect(all.results,efo.exp)) 
  total <- length(ttd.appr$Disease)
  match.rate <- as.character(length(matches)/total)
  match.ratio <- paste(length(matches),":",total,sep="")
  
  if (length(matches) == 0) {
    conf <- as.data.frame(cbind(target, match.rate, match.ratio, n.pred, NA, NA, NA, NA, NA, NA, NA))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected","Status", "EFO_Src")
    return(conf)
  }
  
  match.dis <- unlist(lapply(matches, function(m){efo_df$label[efo_df$id==m&!is.na(efo_df$id)]}))
  
  match.df <- as.data.frame(matrix(nrow=0,ncol=11))
  colnames(match.df) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                          "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected","Status", "EFO_Src")
  for (m in 1:length(matches)) {
    match <- matches[m]
    
    ttd.dis <- as.character(ttd.appr$Disease[grep(match,ttd.appr$EFO)])
    if (length(ttd.dis[grep(match.dis[m], unlist(strsplit(as.character(ttd.dis),'; ')), ignore.case=TRUE)])>0) {
      ttd.dis<-unique(grep(match.dis[m], unlist(strsplit(as.character(ttd.dis),'; ')), ignore.case=TRUE, value=TRUE))
      ttd.dis<-ttd.dis[!is.na(ttd.dis)]
    }
    
    relation.predicted <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
    colnames(relation.predicted) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
    if (TRUE %in% (relations.df$Predicted_Indication %in% match)) {
      res.dis <- c(as.character(top_dis$Disease_Trait[grep(match,top_dis$EFO)])) 
      rrow <- matrix(c(as.character("Exact"), paste(res.dis, sep="/ ")), nrow=1,ncol=2)
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Parents),", ")))
    if (TRUE %in% (tmp %in% match)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Parents)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Parent",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Siblings),", ")))
    if (TRUE %in% grepl(match, relations.df$Siblings)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Siblings)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Sibling",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    } 
    tmp <- as.vector(unlist(strsplit(as.character(relations.df$Children),", ")))
    if (TRUE %in% grepl(match, relations.df$Children)) {
      res.dis <- unique(as.character(top_dis$Disease_Trait[grep(relations.df$Predicted_Indication[grep(match, relations.df$Children)], top_dis$EFO)]))
      rrow <- matrix(c(rep("Child",length(res.dis)), res.dis), ncol=2,nrow=length(res.dis))
      colnames(rrow) <- c("Relation_To_Predicted", "GWAS_Disease_Predicted")
      relation.predicted <- as.data.frame(rbind(relation.predicted, rrow))
    }
    
    stat <- unique(ttd.appr$Status[ttd.appr$Disease==ttd.dis])
    src <- unique(unlist(strsplit(as.character(efo.exp.df$EFO_Src[efo.exp.df$EFO==match]),", ")))  
    
    conf <- cbind(target, match.rate, match.ratio, n.pred, match, match.dis[m], relation.predicted,paste(ttd.dis,collapse="/"), paste(stat, collapse="/"),
                  paste(src,collapse="/"))
    colnames(conf) <- c("Target", "Match_Rate", "Match_Ratio", "N_Predictions", "Matched_EFO", "EFO_Label",
                        "Relation_To_Predicted", "GWAS_Disease_Predicted", "TTD_Disease_Expected", "Status","EFO_Src")
    
    match.df<-unique(as.data.frame(rbind(match.df, conf)))
    
  }
  
  return(match.df)
}

# Function calls
date<-gsub("-","", Sys.Date())

if(!dir.exists(paste("~/files/",date,"_results", sep=""))){
  dir.create(path = paste("~/files/",date,"_results", sep="")) 
}
approved<-filter(ttd_df_ext, Status=="Approved")
approved_hgnc<-unique(approved$HGNC_Symbol[!is.na(approved$HGNC_Symbol) & !grepl('[[:lower:]]|[[:punct:]]', approved$HGNC_Symbol)]) #432->404
clinical <- filter(ttd_df_ext, Status %in% c("Phase 3","Phase 4"))#,"Approved")) 
clinical_hgnc <- unique(clinical$HGNC_Symbol[!is.na(clinical$HGNC_Symbol) & !grepl('[[:lower:]]|[[:punct:]]', clinical$HGNC_Symbol)])
clin_hgnc_only <- setdiff(clinical_hgnc,approved_hgnc)

if(input == 'all') {
  targets<-unique(ttd_df_ext$HGNC_Symbol[!is.na(ttd_df_ext$HGNC_Symbol) & !grepl('[[:lower:]]|[[:punct:]]', ttd_df_ext$HGNC_Symbol)])
  files1<-unique(grep(paste("^",targets,"_",sep="",collapse="|"), 
                      files, value=TRUE))
  if(pstat=="normal") {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelated(files1[i], p.val, nrows, "all"))
    })
  } else {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelatedADJ(files1[i], p.val, nrows, "all"))
    })
  }
} else if(input == 'phase3') {
  targets<-clinical_hgnc
  files1<-unique(grep(paste("^",targets,"_",sep="",collapse="|"), 
                      files, value=TRUE))
  if(pstat=="normal") {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelated(files1[i], p.val, nrows, "phase3"))
    })
  } else {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelatedADJ(files1[i], p.val, nrows, "phase3"))
    })
  }
} else if(input=='approved') {
  targets<-approved_hgnc
  files1<-unique(grep(paste("^",targets,"_",sep="",collapse="|"), 
                      files, value=TRUE)) #removes targets not in files
  if(pstat=="normal") {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelated(files1[i], p.val, nrows, "approved"))
    })
  } else {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelatedADJ(files1[i], p.val, nrows, "approved"))
    })
  }
  
} else {
  targets<-unlist(strsplit(input, ","))
  files1<-unique(grep(paste("^",targets,"_",sep="",collapse="|"), 
                      files, value=TRUE)) #removes targets not in files
  if(pstat=="normal") {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelated(files1[i], p.val, nrows, "approved"))
    })
  } else {
    match.rel <- as.data.frame(foreach(i= 1:length(files1), .combine = rbind) %dopar% {
      library(dplyr)
      library(data.table)
      library(rJava)
      library(rols)
      library(jsonlite)
      library(openxlsx)
      library(S4Vectors)
      library(tidyr)
      as.data.frame(MatchRelatedADJ(files1[i], p.val, nrows, "approved"))
    })
  }
  
}

stopCluster(cl)

# Write files
name.1<-paste("~/files/",date,"_results/",date,"_",input,"_",nrows, "_",threshold,"_",as.character(p.val),sep="")
write.table(match.rel, file = paste(name.1,".txt",sep=""), quote=F, sep="\t", row.names = F)
openxlsx::write.xlsx(match.rel, file = paste(name.1,".xlsx",sep=""), quote=F, col.names = T, row.names = F)
