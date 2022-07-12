get_first_number<-function(x,n=1){
  #used to get the value of a variable in input file, where there is an example value (or any other value) in the same line which we need to ignore
  xx<-as.double(unlist(str_split(x," |\t|,")))
  thisNum<-xx
  if(length(is.na(xx))==1){
    if(is.na(xx)){
      yy<-gsub("([^\\d*\\.?\\d*$])","#",x,perl=TRUE)
      yyy<-unlist(str_split(yy,"#"))
      xPos<-grep("[^\\d]",yyy)[n]
      thisNum<-as.numeric(yyy[xPos])
    }
  }else{
    if(n=="all"){
      thisNum<-as.double(unlist(str_split(xx[!is.na(xx)]," ")))
    }else{
      thisNum<-as.double(unlist(str_split(xx[!is.na(xx)]," "))[n])
    }
  }
  return(thisNum)
}


# 
# 
# 
# 
# sci1<-paste(allnumbers,collapse="e-")
# test1<-grep(sci1,xsplit)
# if(length(test1)>0){
#   if(test1==1){
#     thisNum<-as.double(sci1)
#   }
# }
# sci2<-paste(allnumbers,collapse="e+")
# test2<-grep(sci2,xsplit)
# if(length(test2)>0){
#   if(test2==1){
#     thisNum<-as.double(sci2)
#   }
# }
# if(is.na(thisNum)){
# # yy<-gsub("([^\\d*\\.?\\d*$])","#",x,perl=TRUE)
#   yyy<-unlist(str_split(yy,"#"))
#   xPos<-grep("[^\\d]",yyy)[n]
#   thisNum<-as.numeric(yyy[xPos])
# }