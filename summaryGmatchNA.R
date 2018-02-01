summaryGmatch <- function(DATA,GROUP,Var){
    if(!is.factor(DATA[,Var])){
        pvalue <- with(DATA, t.test(DATA[,Var] ~ DATA[,GROUP]))$p.value
        summary <- aggregate(DATA[,Var]~DATA[,GROUP],DATA,function(x) c(mean = mean(x) ,var = var(x)))
        #by(DATA[,Var],DATA[,GROUP], function(x) mean(x,na.rm = T))
        smd <- smd(DATA, Var, GROUP)
    }else{
        summary <- table(DATA[,GROUP],DATA[,Var])
        pvalue<- chisq.test(DATA[,GROUP], DATA[,Var])$p.value
        smd <- smd(DATA, Var, GROUP)
    }
    

# #=== Tabling pvalues and SMD for output ====
# 
# pvals <-
#         matrix(c(
#                 pvalMotion,
#                 pvalAge,
#                 pvalNVIQ,
#                 pvalGender,
#                 pvalHandedness
#         ),
#         ncol = 5)
# colnames(pvals) <-
#         c("RMSD.PRE.censoring",
#           "Age",
#           "WASI.NVIQ",
#           "Gender",
#           "Handedness"
#         )
# rownames(pvals) <- "Pvals"
# pvals <- as.table(pvals)
# pvals

# Calculate standardized mean difference
# smdMotion <- smd(DATA, RMSD.PRE.censoring)
# smdAge <- smd(DATA, "Age")
# smdNVIQ <- smd(DATA, "WASI.NVIQ")
# smdGender <- smd(DATA, "Gender")
# smdHandedness <- smd(DATA, "Handedness")
# 
# #===Tabling smd for output ===
# smd <-
#         matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
#                ncol = 5)
# colnames(smd) <-
#         c("RMSD.PRE.censoring",
#           "Age",
#           "WASI.NVIQ",
#           "Gender",
#           "Handedness"
#         )
# rownames(smd) <- "SMD"
# smd <- as.table(smd)
# #return(smdMotion )#= smdMotion)#,smdAge = smdAge, smdOpt = smdOpt))
return(list(smd = print(smd),print(summary),pval=print(pvalue)))
#cat("pvals", "smdOpt","table(DATA$group)","\n")
}
