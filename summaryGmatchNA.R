summaryGmatch <- function(DATA){
        GROUP = DATA$group
        summary <- aggregate(.~GROUP,DATA,function(x) c(mean = mean(x), var = var(x),range = range(x)))
        table(GROUP,DATA$Handedness)
        table(GROUP,DATA$Gender)
        
#=== Tabling pvalues and SMD for output ====
pvalMotion <-
        (t.test(RMSD.PRE.censoring ~ group, data = DATA))$p.value
pvalAge <- (t.test(Age ~ group, data = DATA))$p.value
pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = DATA))$p.value
pvalGender <- chisq.test(DATA$group, DATA$Gender)$p.value
pvalHandedness <-
        chisq.test(DATA$group, DATA$Handedness)$p.value

pvals <-
        matrix(c(
                pvalMotion,
                pvalAge,
                pvalNVIQ,
                pvalGender,
                pvalHandedness
        ),
        ncol = 5)
colnames(pvals) <-
        c("RMSD.PRE.censoring",
          "Age",
          "WASI.NVIQ",
          "Gender",
          "Handedness"
        )
rownames(pvals) <- "Pvals"
pvals <- as.table(pvals)
pvals

# Calculate standardized mean difference
smdMotion <- smd(DATA, RMSD.PRE.censoring)
smdAge <- smd(DATA, "Age")
smdNVIQ <- smd(DATA, "WASI.NVIQ")
smdGender <- smd(DATA, "Gender")
smdHandedness <- smd(DATA, "Handedness")

#===Tabling smd for output ===
smd <-
        matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
               ncol = 5)
colnames(smd) <-
        c("RMSD.PRE.censoring",
          "Age",
          "WASI.NVIQ",
          "Gender",
          "Handedness"
        )
rownames(smd) <- "SMD"
smd <- as.table(smd)
#return(smdMotion )#= smdMotion)#,smdAge = smdAge, smdOpt = smdOpt))
return(list(smdMotion = smdMotion,smdAge = smdAge, smdAge = smdAge,
            smdNVIQ =smdNVIQ, smdGender= smdGender, smdHandedness = smdHandedness,
            print(summary),print(pvals),print(smd),print(table(DATA$group))))
#cat("pvals", "smdOpt","table(DATA$group)","\n")
}