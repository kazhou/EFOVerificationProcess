# Build sh script for 03_matchRelations.R to submit jobs to queue

library(dplyr)

args<-commandArgs(TRUE)

targets <- as.character(args[1])
threshold <- as.character(args[2])
p.vals <- as.character(args[3]) 
pstat <- as.character(args[4])
nrows <- as.character(args[5])
date.1 <- as.character(args[6]) #YYYYMMDD

date<-gsub("-","", Sys.Date())
filename<-paste("~/sh/",date,"_",threshold,"_",p.vals,"_ALL", ".sh", sep="")
l<-paste("\n#$ -N ","Par_", p.vals,"_",threshold,"_",nrows,"_", date,sep="")
endL<-paste("\nRscript ~/scripts/03_matchRelations.R",targets,threshold,p.vals, pstat, nrows, date.1)
cat("#!/bin/sh
# Grid Engine options",l,
      "\n#$ -cwd
#$ -S /bin/bash
#$ -V
module load sge
module load R/3.5.0\n", 
    as.character(endL), file=filename, append=FALSE)
system(paste("qsub", filename))
