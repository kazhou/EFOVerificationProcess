# Code to scrape TTD website for target-indication information, as it is not available in the downloads

library(dplyr)
library(biomaRt)
library(xlsx)
#for html scraping
library(stringr)
library(rvest)
library(data.table)
library(doParallel)
library(foreach)
library(rols)
cl<-makeCluster(detectCores()-1) 
registerDoParallel(cl)

#get raw target->disease data
target_dis <- read.delim(file="~/src/P1-05-Target_disease.txt", header=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE) #6547
# save(target_dis, file="~/src/target_dis_full.rds")
target_dis0 <- target_dis
target_dis <- dplyr::select(target_dis, TTDTargetID, Target_Name, Indication)

#for each row/TTD_ID, find status of target-disease pair
targIDs <- as.character(pull(target_dis, TTDTargetID)) #6547
uniqIDs <- unique(targIDs) #1900

drug.info <- data.frame(matrix(vector(),ncol=4, nrow=0))
colnames(drug.info) <- c("Target", "Drug", "Status", "Disease")  
drug.info <- foreach (i=1:length(uniqIDs), .combine=rbind) %dopar% {
  library(dplyr)
  library(UniProt.ws)
  library(AnnotationDbi)
  library(biomaRt)
  library(xlsx)
  #for html scraping
  library(stringr)
  library(rvest)

  target <- uniqIDs[[i]] #uniqIDs[[i]]
  url <- paste("https://db.idrblab.org/ttd/target/", target, sep="") #view-source:, ttd_row[1]
  pg <- read_html(url)
  stat <- pg%>% html_nodes(xpath='//*[@width="11%"]')%>% html_text() #status
  dd <- pg%>% html_nodes(xpath='//*[@width="30%"]')%>% html_text() #odd = drug, even = disease

  if (length(stat)!=0) {   
    drug <- unlist(rep(NA, length(stat)))
    dis <- unlist(rep(NA, length(stat)))
    for (i in 1:length(stat)) {
      drug[i] <- dd[2*i-1]
      dis[i] <- dd[2*i]
    }
    tg <- as.vector(unlist(rep(target, length(stat))))
    df <- as.data.frame(cbind(tg, drug, stat, dis))
  }
  colnames(df) <- c("Target", "Drug", "Status", "Disease") 
  return(df)
}

stopCluster(cl)

#1808 unique targets
#10838 obs in drug.info

drug.info <- unique(drug.info)
save(drug.info, file="~/src/drug_info.rds")
diff <- setdiff(uniqIDs, drug.info$Target) #92

# drug.info<-filter(drug.info, !is.na(Drug))

#merge drug.info and final cut, merge NA disease w Indication
target.indications <- full_join(drug.info, finalCut, by=c("Target"="TTD_Target_ID", "Disease"="Indication")) 
success <- filter(target.indications, Target_Type=="Successful target")
clinical <- filter(target.indications, Target_Type=="Clinical Trial target") 
research <- filter(target.indications, Target_Type=="Research target")

# post-edit/modification
redo < -c(success$Target,clinical$Target) 
redo <- unique(redo) 
target.indications <- dplyr::select(target.indications, Target, HGNC_Symbol, Target_Name, Disease, Status)
df <- target.indications[order(target.indications[,'Target']),]
df <- df[!duplicated(df[c("Target_Name","Disease")]),] 
na_name <- filter(df, is.na(Target_Name))
# length(unique(na_name$Target)) #167
# df <- setDT(df)[df[is.na(Target_Name)], Target_Name := i.Target_Name, on = c("Target", "Target_Name")]
# df$Target_Name[is.na(df$Target_Name)] <- as.character(df$Target[is.na(df$Target_Name)])
for (i in 1:length(unique(df$Target))) {
  target <- unique(df$Target)[[i]]
  names <- df$Target_Name[df$Target==target]
  hgnc <- df$HGNC_Symbol[df$Target==target]
  len <- length(names)  
  len1 <- length(hgnc)  
  names <- unique(names)
  name <- names[!is.na(names)]
  new_names <- rep(name, len)
  hgnc <- unique(hgnc)
  hgnc <- hgnc[!is.na(hgnc)]
  new_hgnc <- rep(hgnc, len1)
  df$Target_Name[df$Target==target] <- new_names
  if (length(new_hgnc)!=0) {
    df$HGNC_Symbol[df$Target==target] <- new_hgnc
  }
}
df$Status <- as.character(df$Status)
df$Status[is.na(df$Status)] <- "Research"

approved <- filter(df, Status == "Approved") #1567 -> 1570
approved_hgnc <- filter(approved, !is.na(HGNC_Symbol)) #890 -> 1219
df_hgnc<- filter(df, !is.na(HGNC_Symbol)) #4623 ->5490
approved_na <- filter(approved, is.na(HGNC_Symbol)) #890 ->351
df_na <- filter(df, is.na(HGNC_Symbol)) #2207 -> 1359

# reappend ICD columns to df
icd <- dplyr::select(target_dis0, Indication, ICD9, ICD10)
# load("180702_df_final")
df <- unique(left_join(df, icd, by=c("Disease"="Indication"))) #6849
# save(df, file="180703_df_icd")

# webscrape for df_na
redoNA <- unique(df_na$Target) 
for (i in 1:length(redoNA)) { 
  target <- redoNA[[i]] 
  url <- paste("https://db.idrblab.org/ttd/target/", target, sep="") #view-source:, ttd_row[1]
  pg <- read_html(url)
  gene <- pg%>% html_nodes(xpath='//*[@class="target__field-gene-name"]')%>% html_text() #gene name, if exists
  if (!isEmpty(gene)) {
    df$HGNC_Symbol[df$Target==target] <- gene
  }
}

df_inhuman <- df[grep('[a-z]', df$HGNC_Symbol),] #200


df <- unique(df)
# rownames(df)<-NULL
save(df, file=".~/src/df_icdFinal")

stats <- sort(levels(as.factor(df$Status)))           

save(target.indications, file="~/src/target_indications.rds")
save(df, file="~/src/df_final")
write.table(df, file="~/src/df_final.txt", sep="\t", quote=F, row.names=F)
write.xlsx2(df, file="~/src/df_icd.xlsx",col.names = TRUE, row.names = FALSE)
