#Dowloaded CC founder contribution tables from:
#http://csbio.unc.edu/CCstatus/index.py?run=FounderProbs
#Build38
#The two errors are as follows: 
#On chromosome 5, there are two markers (SAbGeoEUCOMM001 SAbGeoEUCOMM002) that should not be given genome position. 
#On chromosome 13, there is a problem with the last set of approximately 80 markers that result in an inconsistent 
#pattern of founder haplotype reconstruction. (Beginning near UNC23486670 CH13:118447298 to UNC23498758 CH13:119480991)
#Mask two positions.

#Submit: sbatch --mem=50GB --time=12:00:00 -o submit_FounderStrTable.out --wrap "R --file=GenerateFounderStrTable.R" 

setwd("/pine/scr/j/x/jingxue/CCSNP_founder110218")

library(data.table)

generateFdStr<-function(fdpath){
        fd=fread(fdpath)
        #Convert founder probabilities to 0-1 table with a cutoff of 0
        convert<-function(CC){
                CC_bino=apply(CC[,-c(1:3)],1:2,function(x) ifelse(x>0,1,0))
                cbind(CC[,1:3],data.frame(CC_bino))
        }
        fd_bino=convert(fd)
        
        #Convert the founder probability table to letter table from the 0-1table
        letter<-function(CC_bino){
                name<-colnames(CC_bino[,-c(1:3)])
                len=nrow(CC_bino)
                letter=vector(mode="character",length=len)
                for (i in 1:len){
                        x=paste(rep(name,CC_bino[i,-c(1:3)]),collapse="")
                        letter[i]=paste(sort(unique(unlist(strsplit(x,"")))),collapse="")
                        rm(x)
                }
                data.frame(CC_bino[,1:3],letter)
        }
        fd_letter=letter(fd_bino) 
        fd_letter
}

founder.list=list.files(path=paste(getwd(),"/b38FounderTable",sep = ""))
founder.list
filenames=unlist(lapply(as.list(founder.list),function(x) substring(x,1,5)))
fullpaths=paste(getwd(),"b38FounderTable",founder.list,sep="/")
fdStr.list=lapply(fullpaths,generateFdStr)
names(fdStr.list)=filenames
save(fdStr.list,file="fdStr.list.RData")