setwd("/pine/scr/j/x/jingxue/snp142_CC001_CC011")

#formatted snp142 file under /proj/folamilb/projects/Jan16MethylSeq/dbSNP142/snp142_formatted/dbSNP142_strain_noHeader_simple.bed
#colnames: "CHROM","POS","POS","ID","REF","ALT",LETTERS[1:8]
#file does not have colnames, will add in the function
#The goal is to annotate a formatted snp142 file with CC001 and CC011 genotypes 
#so that SNPs between reference mm10  and CC001/CC011newly can be removed from recently called c1c2 mm10 methylseq CpGs

###The following part saved as Rscript file "Annot_snp142.R" and submitted to longleaf.
#Submit: sbatch --mem=160GB -n 4 --time=7-00:00:00 -o submit_Annot_snp142.out --wrap "R --file=Annot_snp142.R"

library(data.table)
library(plyr)
#Make a function to annotate founder strains to SNP file
annotatFd_bino=function(fd,SNP_fullpath){
  #Load SNP file at VD genes, mm10; replace B6 genotypes from REF genotype letter to "0/0" for later steps
  SNP=fread(SNP_fullpath)
  colnames(SNP)=c("CHROM","POS","POS","ID","REF","ALT",LETTERS[1:8])
  SNP[, B := "0/0"]
  
  
  len.SNP=nrow(SNP)
  
  #annotate SNP file with founders at the CC strains
  locationList=list()
  for (i in 1:len.SNP) {
    chr.match=paste("chr",fd$chromosome,sep="")==SNP$CHROM[i]
    coor.match.bfr=fd$position.B38.<SNP$POS[i]
    coor.match.aft=fd$position.B38.>SNP$POS[i]
    rows.bfr=which(chr.match&coor.match.bfr)
    rows.aft=which(chr.match&coor.match.aft)
    
    if(length(rows.bfr)==0)
    {row.aft=min(rows.aft)
    row.bfr=row.aft}
    else if(length(rows.aft)==0)
    {row.bfr=max(rows.bfr)
    row.aft=row.bfr}
    else{   
      row.bfr=max(rows.bfr)
      row.aft=min(rows.aft)
    }
    #fd: marker chromosome position.B38. letter
    
    toAttach=cbind(fd[row.bfr,3:4], fd[row.aft,3:4])
    locationList[[i]]=toAttach
    
    rm(chr.match,coor.match.bfr,coor.match.aft,row.bfr,row.aft,toAttach)
  }
  
  location=ldply(locationList,data.frame) #convert list to a data.frame
  colnames(location)=c("fd.B38.bfr","fd.bfr","fd.B38.aft","fd.aft")
  SNP_fd=data.table(SNP,location) #Now SNP file is annotated with CC founder before and after SNP.
  
  #determine if genotye at each SNP is the same with ref 
  RefCompare<-function(SNP_fd){
    x1=as.character(unlist(SNP_fd[16]))
    x2=as.character(unlist(SNP_fd[18]))      
    fd.gt=unlist(SNP_fd[7:14])
    fd.gt.bino=unlist(lapply(fd.gt,function(x) ifelse((substr(x,1,1)==0)&(substr(x,3,3)==0),1,0)))

    #gt is the output 0/1 genotype, 0 means different than ref, 1 means having reference genotype
    if(prod(fd.gt.bino)==1){gt=1}else if(nchar(x1)*nchar(x2)==0){gt=0} else{
      # str has the collapsed founder strains in letters at markers before and after
      x=paste(unlist(SNP_fd[16]),unlist(SNP_fd[18]),collapse="") 
      str=sort(unique(unlist(strsplit(x,""))))
      gt=prod(fd.gt.bino[which(LETTERS %in% str)])
      }
    gt
  }     
  
  # Not necessary here: Coerce data.table to data.frame
  # setDF(SNP_fd, rownames=NULL)
  SNP_gt=apply(SNP_fd,1,RefCompare)
  data.table(SNP_fd,genotype=SNP_gt)
}

#load formatted founder contribution table and select CC001 and CC011
load("/proj/folamilb/projects/Jan16MethylSeq/dbSNP142/founderContribution.B38/fdStr.list.RData")

#Note genotype column: 1 means having reference genotype, 0 means otherwise
snp142_CC001=annotatFd_bino(fd=fdStr.list$CC001,SNP_fullpath="/proj/folamilb/projects/Jan16MethylSeq/dbSNP142/snp142_formatted/dbSNP142_strain_noHeader_simple.bed")
snp142_CC011=annotatFd_bino(fd=fdStr.list$CC011,SNP_fullpath="/proj/folamilb/projects/Jan16MethylSeq/dbSNP142/snp142_formatted/dbSNP142_strain_noHeader_simple.bed")

save(list(snp142_CC001,snp142_CC011),"snp142_CC001_CC011.list.RData")

CC001vsCC011=apply(data.table(snp142_CC001$genotype,snp142_CC011$genotype),1,prod)
table(CC001vsCC011)
snp142_c1c2=data.table(snp142_CC001[,1:6],CC001vsCC011)
snp142_c1c2_toexclude=snp142_c1c2[CC001vsCC011==0,]
colnames(snp142_c1c2_toexclude)=c("chr","start","end","SNP.ID","REF","ALT","CC001vsCC011")
fwrite(snp142_c1c2_toexclude,file=("snp142_c1c2_toexclude"))

