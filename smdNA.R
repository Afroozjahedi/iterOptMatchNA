smd = function(DATA,Var,GROUP){

  #pars <- as.list(match.call())
  if(is.factor(DATA[,Var])){
    #if (is.factor(DATA[ASD,as.character(Var)])){
    mASD = mean(as.numeric(DATA[DATA[,GROUP]==1,Var],na.rm = TRUE))
    mTD = mean(as.numeric(DATA[DATA[,GROUP]!=1,Var],na.rm = TRUE))
    varASD = var(as.numeric(DATA[DATA[,GROUP]==1,Var],na.rm = TRUE))
    varTD = var(as.numeric(DATA[DATA[,GROUP]!=1,Var],na.rm = TRUE))
    stdMeanDiff = (mASD - mTD)/(sqrt((varASD + varTD)/2))
    #cat("SMD",pars$Data,pars$Var,"=", stdMeanDiff,"\n")
    return(stdMeanDiff)
  }else {
  mASD = mean(DATA[DATA[,GROUP]==1,Var],na.rm = TRUE)
  mTD = mean(DATA[DATA[,GROUP]!=1,Var],na.rm = TRUE)
  varASD = var(DATA[DATA[,GROUP]==1,Var],na.rm = TRUE)
  varTD = var(DATA[DATA[,GROUP]!=1,Var],na.rm = TRUE)
  stdMeanDiff = (mASD - mTD)/(sqrt((varASD + varTD)/2))
  #cat("SMD",pars$Data,pars$Var,"=", stdMeanDiff,"\n")}
  return(stdMeanDiff)}
}
#debug(smd)
 #smd(nonOutlierLowADOS,RMSD.PRE.censoring)
#smdGender = smd(nonOutlierLowADOS,"Gender")
 
 #b <-  function(data,name) {
   
   ## match.call return a call containing the specified arguments 
   ## and the function name also 
   ## I convert it to a list , from which I remove the first element(-1)
   ## which is the function name
   
  # pars <- as.list(match.call()[-1])
  # data[,as.character(pars$name)]
   #cat(pars$name,"smd is",mean( data[,as.character(pars$name)]))
   #sd(data[,as.character(pars$name)])
   
# }

 