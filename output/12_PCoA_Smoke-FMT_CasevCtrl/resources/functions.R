
##### This function parses out concatenated taxa names. #####
parseQIIME2Taxa<-function(table){
  
  bugNameSplit=gsub(pattern =  "d__", x = table$bugName,replacement = "")
  bugNameSplit=gsub(pattern =  "[.].__", x = bugNameSplit,replacement = "_/_")
  other <- paste(replicate(count,"Other"), collapse = "_/_")
  bugNameSplit=gsub(pattern ="^Other$",x = bugNameSplit,replacement = other)
  string <- strsplit(as.character(bugNameSplit),split = "_/_")
  
  temp_string=do.call(rbind,string)
  table <- cbind(temp_string,table)
  
  colnames(table)[colnames(table)=="1"] <- "Domain"
  colnames(table)[colnames(table)=="2"] <- "Phylum"
  colnames(table)[colnames(table)=="3"] <- "Class"
  colnames(table)[colnames(table)=="4"] <- "Order"
  colnames(table)[colnames(table)=="5"] <- "Family"
  colnames(table)[colnames(table)=="6"] <- "Genus"
  
  return(table)
  
  }


##### This function rounds p values based on value #####
roundP <- function(pvalue) {
  
  if (length(pval) == 0) {
    p = "p = NA"
  } else if (pvalue==0){
    p = "p < 2.2e-16"
  } else if (pvalue < 0.001) {
    p <- "p < 0.001"
  } else {
    p <- paste0("p = ", round(pvalue, digits = 3))
  }
  
}


##### This function assigns significance *s based on p value #####
sigStars <- function(pvalue) {
  
  if (pvalue < 0.05) {
    pStar <- "*"
  } else if (pvalue < 0.01){
    pStar <- "**"
  } else if (pvalue < 0.001){
    pStar <- "***"
  } else {
    pStar <- ""
  }
  
}


##### This function corrects the sign of p-values based on coefficient sign #####
getlog10p<-function(pval,coefficient){
  
  log10p<-sapply(1:length(pval), function(x){
    
    if (coefficient[x]>0) return(-log10(pval[x]))
    else return(log10(pval[x]))
  })
}


##### This function gives the log10 p-values from studies or time points that we want to compare #####
compareStudies<-function(level, file1, study1Name, study1Tumor, study1Term, file2, study2Name, study2Tumor, study2Term){
  
  myT1<-read.table(file1, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  myT1 <- myT1[myT1$Study == study1Name,]
  myT1 <- myT1[myT1$Tumor == study1Tumor,]
  myT2<-read.table(file2, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  myT2 <- myT2[myT2$Study == study2Name,]
  myT2 <- myT2[myT2$Tumor == study2Tumor,]
  
  common<-intersect(myT1$BugName,myT2$BugName) 
  myT1c<-myT1[myT1$BugName %in% common,]
  myT2c<-myT2[myT2$BugName %in% common,]
  myT2c<-myT2c[match(myT1c$BugName,myT2c$BugName),] # reorder taxa so in the same order as myT1c
  
  p1Col <- paste0("p_", study1Term)
  p2Col <- paste0("p_", study2Term)
  naColumns<-c(which(is.na(myT1c[,p1Col])),which(is.na(myT2c[,p2Col])))
  
  if(length(naColumns)>0){
    myT1c<-myT1c[-naColumns,]
    myT2c<-myT2c[-naColumns,]
  }
  
  p_1<-paste0("p_", study1Term)
  p_2<-paste0("p_", study2Term)
  s_1<-paste0("s_", study1Term)
  s_2<-paste0("s_", study2Term)
  
  pval1<-getlog10p(myT1c[,p_1],myT1c[,s_1])
  pval2<-getlog10p(myT2c[,p_2],myT2c[,s_2])
  
  df<-data.frame(pval1,pval2,bugName=myT2c$BugName)
  
  return(df) 
}


##### This function performs Spearman test between studies #####
correlationBetweenStudies<-function(df){
  
  myList<-list()
  test<-cor.test(df[,"pval1"],df[,"pval2"],method = "spearman")
  myList[[1]]<-test$estimate
  myList[[2]]<-test$p.value
  return(myList)
}



##### This function gives the log10 p-values from studies or time points that we want to compare  for FMT studies#####
compareStudiesFMT<-function(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term){
  
  myT1<-read.table(file1, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  myT1 <- myT1[myT1$Study == study1Name,]
  myT1 <- myT1[myT1[,TermType] == study1Term,]
  
  myT2<-read.table(file2, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  myT2 <- myT2[myT2$Study == study2Name,]
  myT2 <- myT2[myT2[,TermType] == study2Term,]
  
  common<-intersect(myT1$BugName,myT2$BugName) 
  myT1c<-myT1[myT1$BugName %in% common,]
  myT2c<-myT2[myT2$BugName %in% common,]
  myT2c<-myT2c[match(myT1c$BugName,myT2c$BugName),] # reorder taxa so in the same order as myT1c
  
  p1Col <- paste0("p_", comparison)
  p2Col <- paste0("p_", comparison)
  naColumns<-c(which(is.na(myT1c[,p1Col])),which(is.na(myT2c[,p2Col])))
  
  if(length(naColumns)>0){
    myT1c<-myT1c[-naColumns,]
    myT2c<-myT2c[-naColumns,]
  }
  
  p_1<-paste0("p_", comparison)
  p_2<-paste0("p_", comparison)
  s_1<-paste0("s_", comparison)
  s_2<-paste0("s_", comparison)
  
  pval1<-getlog10p(myT1c[,p_1],myT1c[,s_1])
  pval2<-getlog10p(myT2c[,p_2],myT2c[,s_2])
  
  df<-data.frame(pval1,pval2,bugName=myT2c$BugName)
  
  return(df) 
}


##### This function generate p-value versus p-value plot #####
plotPairwiseSegments<-function(df, xlab, ylab, coeficient, p, arrowRight, arrowLeft){
  
  df1<-df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]<log10(0.05)) | (df[,"pval1"]>-log10(0.05) & df[,"pval2"]>-log10(0.05)),]
  
  if (p==0){
    p.1 = "< 2.2e-16"
  } else {
    p.1 = paste("=",format(p,digits = 3))
  }
  
  theme_set(theme_classic())
  
  lim = ceiling(max(abs( c( df$pval1, df$pval2) )))
  inc <- lim / 2
  dec <- -(lim / 2)
  
  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+
    geom_point(size=1, shape = 1)+
    geom_segment(aes(x = log10(0.05), y = -lim, xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = -lim, y = log10(0.05), xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = -log10(0.05), y = lim, xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = lim, y = -log10(0.05), xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
    
    geom_segment(aes(x = log10(0.05), y = lim, xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = -lim, y = -log10(0.05), xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = -log10(0.05), y = -lim, xend = -log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = -log10(0.05), y = log10(0.05), xend = lim, yend = log10(0.05)), linetype = "dashed", color = "black")+
    labs(x=xlab,y=ylab,
         title = paste0("Spearman coefficient = ",format(coeficient,digits =3),"\nAdjusted p ",p.1))+
    geom_text_repel(data=df1,aes(x=,y=pval2,label=bugName),segment.colour="red",size=2.5,min.segment.length = 0,
                    segment.color="grey",segment.size=0.2)+
    geom_point(data=df1, aes(x=pval1,y=pval2), color='blue', size=2)+
    annotate("segment", x=-log10(0.05), y=-(lim+1), xend=lim, yend=-(lim+1),
             col="black", arrow=arrow(length=unit(0.3, "cm"))) +
    annotate("text", x=inc, y=-(lim+1.5), label = arrowRight) +
    annotate("segment", x=log10(0.05), y=-(lim+1), xend=-lim, yend=-(lim+1),
             col="black", arrow=arrow(length=unit(0.3, "cm")))+
    annotate("text", x=dec, y=-(lim+1.5), label = arrowLeft)
  
  
}


