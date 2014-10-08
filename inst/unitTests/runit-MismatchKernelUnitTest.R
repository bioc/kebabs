##
## RUnit test for mismatch kernel
##

library(RUnit)

testGetExplicitRepAll <- function(xc, kernel, expectedResult, silent, sel,
                                  features=NULL, colSubset=NULL, 
                                  annCharset=NULL, ann1=NULL, ann2=NULL)
{
    dnaColNames <- paste(rep(c("A","C","G","T"),each=4),
                         rep(c("A","C","G","T"),4), sep="")
    rnaColNames <- paste(rep(c("A","C","G","U"),each=4),
                         rep(c("A","C","G","U"),4), sep="")
    aaColNames <- paste(rep(c("A","C","D","E","F","G","H","I","K","L","M",
                              "N","P","Q","R","S","T","U", "V","W","Y"),
                            each=21),
                        rep(c("A","C","D","E","F","G","H","I","K","L","M",
                            "N","P","Q","R","S","T","U", "V","W","Y"),21),
                        sep="")
    
    testGetExplicitRepU(xc=xc, kernel=kernel, silent=silent, 
                        expectedResult=expectedResult, sel=sel,
                        features=features, colSubset=colSubset, 
                        annCharset=annCharset, ann1=ann1, ann2=ann2,
                        dnaColNames=dnaColNames, rnaColNames=rnaColNames,
                        aaColNames=aaColNames)
}

testPosIndepKernelKernelMatrix <- function(silent)
{
    testName(silent, "Mismatch Kernel Position Independent")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("GGAACGCCCGCTT", "CCGATACCTA", "GACCGAGTTATACTCT")
    names(xc) <- namesX
    yc <- c("ACCCGCTT", "CCGGTTAATA")
    names(yc) <- namesY
    zc <- "ACACGGAGAAGACCATT"
    names(zc) <- namesZ
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - y")
    
    symres1u <- defMat(matrix(c(482,341,547,341,271,426,547,426,709), nrow=3),
                       namesX, namesX)
    symres2u <- defMat(matrix(c(187,184,184,261), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(287,209,328,319,250,417), nrow=3), 
                       namesX, namesY)

    symres1au <- defMat(matrix(c(1978,1327,2043,1327,1087,1633,2043,1633,2715),
                               nrow=3), namesX, namesX)
    
    symres2au <- defMat(matrix(c(867,643,643,1009), nrow=2), namesY, namesY)
    
    asymresau <- defMat(matrix(c(1205,855,1263,1135,930,1590), nrow=3), 
                        namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu),
                                    as.KernelMatrix(symres1au), 
                                    as.KernelMatrix(symres2au), 
                                    as.KernelMatrix(asymresau))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    testGetKernelMatrixAll(xc, yc, mm, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized x - y",
                               "with subsetting"))
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), mm,
                           expectedResUnnormalized, silent, 
                           subsetX=c(3:5), subsetY=c(3:4))
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - y")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)),
                                  normalizeMatrix(symres1au), 
                                  normalizeMatrix(symres2au), 
                                  normalizeMatrix(asymresau, diag(symres1au),
                                                             diag(symres2au)))
    
    mm <- mismatchKernel(k=2, m=1, normalized=TRUE)
    testGetKernelMatrixAll(xc, yc, mm, expectedResNormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - z")

    symres1u <- defMat(matrix(c(482,341,547,341,271,426,547,426,709), 
                              nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(852), nrow=1), namesZ, namesZ)
    
    asymresu <- defMat(matrix(c(595,446,734), nrow=3), namesX, namesZ)
    
    symres1au <- defMat(matrix(c(1978,1327,2043,1327,1087,1633,2043,1633,2715),
                               nrow=3), namesX, namesX)
    
    symres2au <- defMat(matrix(c(3470), nrow=1), namesZ, namesZ)
    
    asymresau <- defMat(matrix(c(2278,1704,2740), nrow=3), namesX, namesZ)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu),
                                    as.KernelMatrix(symres1au), 
                                    as.KernelMatrix(symres2au), 
                                    as.KernelMatrix(asymresau))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
    testGetKernelMatrixAll(xc, zc, mm, expectedResUnnormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - z")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)),
                                  normalizeMatrix(symres1au), 
                                  normalizeMatrix(symres2au), 
                                  normalizeMatrix(asymresau, diag(symres1au),
                                                             diag(symres2au)))
    
    mm <- mismatchKernel(k=2, m=1, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, zc, mm, expectedResNormalized, silent)
    
## TODO $$$ Remove 4 comments for kernel matrix with presence
#    testStepName(silent, "Pos Indep - kmer presence unnormalized")
    
    symres1u <- defMat(matrix(c(16,63,105,63,16,105,105,105,16), nrow=3),
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(16,63,63,16), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(49,49,49,63,63,63), nrow=3), namesX, namesY)

    symres1au <- defMat(matrix(c(152,369,615,369,152,615,615,615,152), nrow=3),
                        namesX, namesX)
    
    symres2au <- defMat(matrix(c(135,318,318,152), nrow=2), namesY, namesY)
    
    asymresau <- defMat(matrix(c(287,287,287,369,369,369), nrow=3), 
                        namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu),
                                    as.KernelMatrix(symres1au), 
                                    as.KernelMatrix(symres2au), 
                                    as.KernelMatrix(asymresau))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE, presence=TRUE)
    
#    testGetKernelMatrixAll(xc, yc, mm, expectedResUnnormalized, silent)
    

#    testStepName(silent, "Pos Indep - kmer presence normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)),
                                  normalizeMatrix(symres1au), 
                                  normalizeMatrix(symres2au), 
                                  normalizeMatrix(asymresau, diag(symres1au),
                                                             diag(symres2au)))
    
    mm <- mismatchKernel(k=2, m=1, normalized=TRUE, presence=TRUE)
    
#    testGetKernelMatrixAll(xc, yc, mm, expectedResNormalized, silent)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "tanslate lower to upper"))
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("GGaacGCcCGctT", "cCGAtacCTA", "GACcGAgTTATactct")
    names(xc) <- namesX
    yc <- c("aCCcgcTT", "ccggtTAaTA")
    names(yc) <- namesY
    zc <- "ACACgGAGaagACCATT"
    names(zc) <- namesZ
    
    symres1u <- defMat(matrix(c(482,341,547,341,271,426,547,426,709), nrow=3),
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(187,184,184,261), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(287,209,328,319,250,417), nrow=3), 
                       namesX, namesY)
    
    symres1au <- defMat(matrix(c(1978,1327,2043,1327,1087,1633,2043,1633,2715),
                        nrow=3), namesX, namesX)
    
    symres2au <- defMat(matrix(c(867,643,643,1009), nrow=2), namesY, namesY)
    
    asymresau <- defMat(matrix(c(1205,855,1263,1135,930,1590), nrow=3), 
                        namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu),
                                    as.KernelMatrix(symres1au), 
                                    as.KernelMatrix(symres2au), 
                                    as.KernelMatrix(asymresau))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE, ignoreLower=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mm, expectedResUnnormalized, silent)

    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "ignore lowercase"))
    
    symres1u <- defMat(matrix(c(41,37,46,37,60,75,46,75,132), nrow=3), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(18,12,12,28), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(16,24,35,12,30,46), nrow=3), namesX, namesY)
    
    symres1au <- defMat(matrix(c(211,139,131,139,264,279,131,279,574), nrow=3),
                        namesX, namesX)
    
    symres2au <- defMat(matrix(c(86,46,46,164), nrow=2), namesY, namesY)
    
    asymresau <- defMat(matrix(c(50,92,120,12,132,216), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu),
                                    as.KernelMatrix(symres1au), 
                                    as.KernelMatrix(symres2au), 
                                    as.KernelMatrix(asymresau))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE, ignoreLower=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mm, expectedResUnnormalized, silent,
                           biovectorsOnly=TRUE )
}

testExplicitRepresentation <- function(silent)
{
    ## UniProt IDs
    namesX <- c("P53554", "Q50EK3", "P15150", "Q64408", "P15538", "P97720", 
                "Q29527", "Q29552", "P15393", "P51663")
    
    ac <- c("MTIASSTASS", "MDVNILTMFV", "MALWAKARVR", "MAFRLKSDVR", 
            "MALRAKAEVC", "MALRAKADVW", "MALRAKAEVC", "MAIWAKAEAW", 
            "MALRVTADVW", "MALWAKARVW")
    
    dc <- c("ATGACAATAGCAAGCAGCACAGCAAGCAGC",
            "ATGGACGTAAACATACTAACAATGTTTGTA",
            "ATGGCACTATGGGCAAAAGCAAGAGTAAGA",
            "ATGGCATTTAGACTAAAAAGCGACGTAAGA",
            "ATGGCACTAAGAGCAAAAGCAGAAGTATGC",
            "ATGGCACTAAGAGCAAAAGCAGACGTATGG",
            "ATGGCACTAAGAGCAAAAGCAGAAGTATGC",
            "ATGGCAATATGGGCAAAAGCAGAAGCATGG",
            "ATGGCACTAAGAGTAACAGCAGACGTATGG",
            "ATGGCACTATGGGCAAAAGCAAGAGTATGG")
    
    names(ac) <- namesX
    names(dc) <- namesX
    
    dna <- DNAStringSet(dc)
    aa <- AAStringSet(ac)
    taa <- translate(dna)
    
    if (!all(taa == aa))
        stop("Error in DNA -> AA translation in Biostrings")
    
    testName(silent, "Explicit Representation for Mismatch Kernel")
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - unnormalized"))

    namesY <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG",
                "GT","TA","TC","TG","TT")
    
    erd1u <- defMat(matrix(
        c(22, 19, 14, 13, 12, 15, 14,  9, 18,  9, 14,  9, 13, 10,  8,  4,
          18, 11, 16, 17, 13,  8,  8, 12, 15,  9,  9, 11, 16, 13, 11, 16, 
          19, 15, 17, 14, 13,  8, 13,  7, 19, 10, 15, 12, 14,  8, 11,  8, 
          19, 14, 16, 16, 15,  8,  9,  9, 16,  9, 13, 12, 15, 10, 12, 10, 
          20, 17, 16, 15, 13,  9, 12,  7, 18,  9, 15, 11, 14,  9, 10,  8, 
          19, 15, 17, 14, 13, 10, 13,  8, 17, 10, 15, 11, 13,  9, 11,  8, 
          20, 17, 16, 15, 13,  9, 12,  7, 18,  9, 15, 11, 14,  9, 10,  8, 
          18, 16, 19, 12, 11,  8, 14,  8, 19,  9, 15, 13, 14,  8, 11,  8, 
          19, 13, 16, 14, 12, 10, 13,  9, 16, 11, 15, 11, 12, 10, 12, 10, 
          17, 14, 18, 13, 11,  8, 14,  8, 18, 10, 15, 13, 13,  9, 12, 10),
          nrow=10, byrow=TRUE), namesX, namesY)
    
    namesY2 <- 
        c("AA","AC","AD","AE","AF","AG","AH","AI","AK","AL","AM","AN","AP",
          "AQ","AR","AS","AT","AU","AV","AW","AY","CA","CC","CD","CE","CF",
          "CI","CK","CL","CM","CN","CR","CS","CT","CV","CW","DA","DC","DD",
          "DE","DF","DG","DH","DI","DK","DL","DM","DN","DP","DQ","DR","DS",
          "DT","DU","DV","DW","DY","EA","EC","ED","EE","EF","EG","EH","EI",
          "EK","EL","EM","EN","EP","EQ","ER","ES","ET","EU","EV","EW","EY",
          "FA","FC","FD","FE","FF","FG","FH","FI","FK","FL","FM","FN","FP",
          "FQ","FR","FS","FT","FU","FV","FW","FY","GA","GC","GD","GE","GF",
          "GI","GK","GL","GM","GN","GR","GS","GT","GV","GW","HA","HC","HD",
          "HE","HF","HI","HK","HL","HM","HN","HR","HS","HT","HV","HW","IA",
          "IC","ID","IE","IF","IG","IH","II","IK","IL","IM","IN","IP","IQ",
          "IR","IS","IT","IU","IV","IW","IY","KA","KC","KD","KE","KF","KG",
          "KH","KI","KK","KL","KM","KN","KP","KQ","KR","KS","KT","KU","KV",
          "KW","KY","LA","LC","LD","LE","LF","LG","LH","LI","LK","LL","LM",
          "LN","LP","LQ","LR","LS","LT","LU","LV","LW","LY","MA","MC","MD",
          "ME","MF","MG","MH","MI","MK","ML","MM","MN","MP","MQ","MR","MS",
          "MT","MU","MV","MW","MY","NA","NC","ND","NE","NF","NG","NH","NI",
          "NK","NL","NM","NN","NP","NQ","NR","NS","NT","NU","NV","NW","NY",
          "PA","PC","PD","PE","PF","PI","PK","PL","PM","PN","PR","PS","PT",
          "PV","PW","QA","QC","QD","QE","QF","QI","QK","QL","QM","QN","QR",
          "QS","QT","QV","QW","RA","RC","RD","RE","RF","RG","RH","RI","RK",
          "RL","RM","RN","RP","RQ","RR","RS","RT","RU","RV","RW","RY","SA",
          "SC","SD","SE","SF","SG","SH","SI","SK","SL","SM","SN","SP","SQ",
          "SR","SS","ST","SU","SV","SW","SY","TA","TC","TD","TE","TF","TG",
          "TH","TI","TK","TL","TM","TN","TP","TQ","TR","TS","TT","TU","TV",
          "TW","TY","UA","UC","UD","UE","UF","UI","UK","UL","UM","UN","UR",
          "US","UT","UV","UW","VA","VC","VD","VE","VF","VG","VH","VI","VK",
          "VL","VM","VN","VP","VQ","VR","VS","VT","VU","VV","VW","VY","WA",
          "WC","WD","WE","WF","WG","WH","WI","WK","WL","WM","WN","WP","WQ",
          "WR","WS","WT","WU","WV","WW","WY","YA","YC","YD","YE","YF","YI",
          "YK","YL","YM","YN","YR","YS","YT","YV","YW")
    
    erd2u <- defMat(matrix(
        c(4,2,2,2,2,2,2,3,2,2,2,2,2,2,2,4,4,2,2,2,2,2,0,0,0,0,1,0,0,0,0,0,4,
          2,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,2,0,0,0,0,0,0,1,0,
          0,0,0,0,0,0,4,2,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,
          2,0,0,0,0,1,0,0,0,0,0,4,2,0,0,2,0,0,0,0,1,0,0,0,0,0,4,2,0,0,2,1,1,
          1,1,1,1,2,1,1,1,1,1,1,1,5,3,1,1,1,1,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
          4,2,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,3,1,1,1,1,1,
          1,2,1,1,1,1,1,1,1,5,2,1,1,1,1,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,
          0,0,0,2,0,0,0,0,1,0,0,0,0,0,4,2,0,0,2,0,0,0,0,1,0,0,0,0,0,4,2,0,0,
          2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,5,3,3,3,3,3,3,4,3,3,3,3,
          3,3,3,5,4,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,6,4,2,2,2,2,2,0,0,
          0,0,1,0,0,0,0,0,4,2,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,
          2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,4,2,0,0,0,0,2,0,0,0,0,1,0,0,0,0,0,4,
          2,0,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,2,0,0,0,0,1,0,1,1,0,1,1,
          1,0,0,1,2,0,1,1,2,1,2,1,1,2,1,2,2,2,1,1,1,1,2,1,2,1,1,0,0,1,0,1,0,
          0,1,0,1,1,1,0,0,0,0,1,0,2,0,0,1,1,2,1,2,1,1,2,1,2,2,2,1,1,1,1,2,1,
          2,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,2,0,0,0,1,0,1,1,0,1,1,1,0,0,1,2,0,
          1,1,2,1,2,1,1,2,1,1,2,2,1,1,1,1,2,1,3,1,1,0,0,1,0,1,0,0,1,0,1,1,1,
          0,0,0,0,1,0,2,0,0,1,1,2,1,2,1,1,2,1,2,2,2,1,1,1,1,1,1,3,1,1,2,2,2,
          2,2,2,2,3,2,3,3,3,2,2,2,2,3,2,4,2,2,1,1,2,1,2,1,1,1,1,2,2,2,1,1,1,
          1,2,1,3,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,2,0,0,0,1,0,1,1,0,1,1,1,0,0,
          1,2,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,2,0,0,0,0,1,0,1,0,0,1,0,
          1,1,1,0,0,0,0,1,0,2,0,0,1,1,2,1,2,1,1,2,1,2,1,2,1,1,1,1,2,1,3,1,1,
          0,0,1,0,1,1,0,1,1,1,0,0,1,2,0,1,1,2,1,2,1,1,2,1,2,2,1,1,1,1,1,2,1,
          3,1,1,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,2,0,0,0,0,1,0,1,1,0,1,1,
          1,0,0,1,2,0,6,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,4,4,3,3,0,0,0,0,0,
          1,1,0,0,2,0,0,1,1,3,0,0,0,0,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1,0,3,0,0,
          0,0,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1,0,3,0,0,0,0,0,0,0,1,1,0,0,0,0,2,
          0,0,0,1,1,0,3,0,0,0,0,0,1,1,0,0,2,0,0,1,1,3,0,0,0,0,0,1,1,0,0,2,0,
          0,1,1,3,0,0,0,0,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1,0,3,1,1,1,1,1,1,1,2,
          2,1,1,1,1,3,1,1,1,2,2,1,4,1,1,1,1,1,1,1,2,2,1,1,1,1,3,1,1,1,2,1,1,
          3,1,1,1,1,1,1,1,2,2,1,1,1,1,3,1,1,1,2,2,1,3,0,0,0,0,0,0,0,1,1,0,0,
          0,0,2,0,0,0,1,1,0,3,0,0,0,0,0,1,1,0,0,2,0,0,1,1,3,0,0,0,0,0,1,1,0,
          0,2,0,0,1,1,4,1,1,1,1,1,1,1,2,2,1,1,1,1,3,1,1,1,1,2,1,3,0,0,0,0,0,
          0,0,1,1,0,0,0,0,2,0,0,0,1,1,0,3,0,0,0,0,0,0,0,1,1,0,0,0,0,2,0,0,0,
          1,1,0,3,0,0,0,0,0,1,1,0,0,2,0,0,1,1,4,1,1,1,1,1,1,1,2,2,1,1,1,1,2,
          1,1,1,2,2,1,3,1,1,1,1,1,1,1,2,2,1,1,1,1,3,1,1,1,2,2,1,3,0,0,0,0,0,
          1,1,0,0,2,0,0,1,1,2,1,2,1,1,1,1,1,2,2,1,1,1,1,3,2,1,1,2,1,1,1,0,1,
          0,1,0,1,1,0,0,2,1,0,1,0,2,1,2,1,2,1,1,1,2,2,1,1,1,1,3,2,1,1,1,1,1,
          1,0,1,0,1,0,0,0,1,1,0,0,0,0,2,1,0,0,1,0,0,2,1,2,1,2,1,1,1,2,2,1,1,
          1,1,2,2,1,1,2,1,1,1,0,1,0,1,0,1,1,0,0,2,1,0,1,0,1,0,1,0,1,0,1,1,0,
          0,2,1,0,1,0,1,0,1,0,1,0,0,0,1,1,0,0,0,0,2,1,0,0,1,0,0,2,1,2,1,2,1,
          1,1,2,2,1,1,1,1,3,1,1,1,2,1,1,2,1,2,1,2,1,1,1,1,2,1,1,1,1,3,2,1,1,
          2,1,1,1,1,2,1,2,1,1,1,2,2,1,1,1,1,3,2,1,1,2,1,1,1,0,1,0,1,0,0,0,1,
          1,0,0,0,0,2,1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,0,2,1,0,1,0,1,0,1,0,1,0,
          1,1,0,0,2,1,0,1,0,2,1,2,1,2,1,1,1,2,1,1,1,1,1,3,2,1,1,2,1,1,2,1,1,
          1,2,1,1,1,2,2,1,1,1,1,3,2,1,1,2,1,1,1,0,1,0,1,0,0,0,1,1,0,0,0,0,2,
          1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,0,2,1,0,1,0,2,1,2,1,2,1,1,1,2,2,1,1,
          1,1,2,2,1,1,2,1,1,1,0,1,0,1,0,0,0,1,1,0,0,0,0,2,1,0,0,1,0,0,1,0,1,
          0,1,0,1,1,0,0,2,1,0,1,0,6,4,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,4,3,3,
          3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
          1,0,0,4,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,1,1,1,3,1,0,1,0,0,0,0,1,
          1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,1,0,1,0,0,
          1,1,0,0,1,0,0,1,0,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,2,1,
          2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,4,2,1,2,1,1,1,1,2,2,1,1,1,1,1,
          1,1,1,2,1,1,3,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,3,1,0,1,0,0,
          0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,1,0,
          1,0,0,1,1,0,0,1,0,0,1,0,3,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,
          3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,0,0,1,1,0,0,
          0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,4,1,1,2,1,1,1,1,2,
          2,1,1,1,1,2,1,1,1,2,1,1,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,
          3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,6,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,
          4,4,3,3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,4,1,2,1,1,1,1,1,2,2,1,1,1,1,2,
          1,1,1,1,2,1,3,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,3,0,1,0,0,0,
          0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,3,0,1,
          0,0,0,1,1,0,0,1,0,0,1,1,3,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,
          3,1,2,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,2,1,4,1,2,1,1,1,1,1,2,2,1,1,
          1,1,1,1,1,1,2,2,1,3,1,2,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,2,1,3,0,1,
          0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
          3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,3,1,2,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,
          2,2,1,3,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,3,0,1,0,0,0,0,0,1,
          1,0,0,0,0,1,0,0,0,1,1,0,3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,4,1,2,1,1,1,
          1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,3,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
          1,1,0,3,0,1,0,0,0,1,1,0,0,1,0,0,1,1,6,4,3,3,3,3,3,3,3,3,3,3,3,3,4,
          3,3,3,4,3,3,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,1,0,1,0,0,0,0,1,1,0,0,
          0,0,1,0,0,0,1,0,0,4,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,1,1,1,3,1,0,
          1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,
          3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
          1,0,0,3,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,4,2,1,2,1,1,1,1,2,
          2,1,1,1,1,1,1,1,1,2,1,1,3,2,1,2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,
          3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,
          0,1,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,3,2,1,2,1,1,1,1,2,2,1,1,1,1,2,
          1,1,1,2,1,1,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,
          0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,4,1,1,
          2,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,1,1,3,1,0,1,0,0,0,0,1,1,0,0,0,0,1,
          0,0,0,1,0,0,3,1,0,1,0,0,1,1,0,0,1,0,0,1,0,8,4,4,4,4,4,4,4,4,4,4,4,
          4,4,4,4,4,4,4,5,4,4,0,0,1,0,1,1,0,0,0,0,0,0,0,2,4,0,0,1,0,0,0,1,1,
          0,0,0,0,0,0,0,0,0,0,2,0,4,1,1,2,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,3,1,
          4,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,0,0,1,0,1,1,0,0,0,0,0,
          0,0,2,4,0,0,1,0,1,1,0,0,0,0,0,0,0,2,5,1,1,2,1,1,1,2,2,1,1,1,1,1,1,
          1,1,1,1,2,1,4,1,1,2,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,3,1,4,0,0,1,0,0,
          0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,1,1,2,1,1,1,2,2,1,1,1,1,1,1,1,1,1,
          1,3,1,4,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,0,0,1,0,1,1,0,0,
          0,0,0,0,0,2,4,0,0,1,0,1,1,0,0,0,0,0,0,0,2,4,0,0,1,0,0,0,1,1,0,0,0,
          0,0,0,0,0,0,0,2,0,4,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,0,0,
          1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,0,0,1,0,1,1,0,0,0,0,0,0,0,2,
          4,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,2,0,4,1,1,2,1,1,1,2,2,1,1,1,
          1,1,1,1,1,1,1,3,1,4,0,0,1,0,1,1,0,0,0,0,0,0,0,2,4,2,2,2,2,2,2,2,2,
          2,2,2,2,2,3,2,3,2,4,3,2,2,0,1,0,0,0,0,1,0,0,1,0,1,2,1,3,1,2,1,1,1,
          1,1,1,2,1,1,1,1,2,1,2,1,2,2,1,2,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,
          2,1,0,2,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,2,1,0,2,0,1,0,0,0,0,1,0,
          0,1,0,1,2,1,2,0,1,0,0,0,0,1,0,0,1,0,1,2,1,2,0,1,0,0,0,0,0,0,1,0,0,
          0,0,1,0,1,0,2,1,0,2,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,2,1,0,3,1,2,
          1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,3,2,1,2,1,2,1,1,1,1,1,1,2,1,1,1,1,2,
          1,2,1,3,2,1,2,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,2,1,0,2,0,1,0,0,0,
          0,1,0,0,1,0,1,2,1,2,0,1,0,0,0,0,1,0,0,1,0,1,2,1,3,1,2,1,1,1,1,1,1,
          2,1,1,1,1,2,1,2,1,2,2,1,2,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,2,1,0,
          2,1,2,1,1,1,1,1,1,2,1,1,1,1,2,1,2,1,3,2,1,2,0,1,0,0,0,0,1,0,0,1,0,
          1,2,1,4,2,3,2,2,2,2,2,2,3,2,2,2,2,3,2,2,2,4,2,2,2,0,1,0,0,0,0,0,0,
          1,0,0,0,0,1,0,1,0,2,1,0,2,0,1,0,0,0,0,1,0,0,1,0,1,2,1,6,3,3,3,3,3,
          3,3,3,3,3,3,3,3,3,3,3,3,4,5,3,3,0,0,0,0,0,1,1,0,0,1,0,0,1,2,3,0,0,
          0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,2,0,3,0,0,0,0,0,0,0,1,1,0,0,0,0,1,
          0,0,0,1,2,0,3,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,2,0,3,0,0,0,0,0,
          1,1,0,0,1,0,0,1,2,3,0,0,0,0,0,1,1,0,0,1,0,0,1,2,3,0,0,0,0,0,0,0,1,
          1,0,0,0,0,1,0,0,0,1,2,0,3,1,1,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,3,1,
          4,1,1,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,2,1,3,1,1,1,1,1,1,1,2,2,1,1,
          1,1,2,1,1,1,2,3,1,3,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,2,0,3,0,0,
          0,0,0,1,1,0,0,1,0,0,1,2,3,0,0,0,0,0,1,1,0,0,1,0,0,1,2,4,1,1,1,1,1,
          1,1,2,2,1,1,1,1,2,1,1,1,1,3,1,3,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
          1,2,0,3,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,2,0,3,0,0,0,0,0,1,1,0,
          0,1,0,0,1,2,4,1,1,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,2,2,1,3,1,1,1,1,1,
          1,1,2,2,1,1,1,1,2,1,1,1,2,3,1,3,0,0,0,0,0,1,1,0,0,1,0,0,1,2), 
         nrow=10, byrow=TRUE), namesX, namesY2)
    
    expectedResUnnormalized <- list(erd1u, erd2u) 
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)

    testGetExplicitRepAll(dc, mm, expectedResUnnormalized, silent)

    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), mm, expectedResUnnormalized, 
                         silent, sel=c(3:12))
    

    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    mm <- mismatchKernel(k=2, m=1, normalized=TRUE)
    
    testGetExplicitRepAll(dc, mm, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - unnormalized"))
    
    erd1up <- erd1u
    erd1up[which(erd1up > 0)] <- 1
    erd2up <- erd2u
    erd2up[which(erd2up > 0)] <- 1
    
    expectedResPresenceUnnormalized <- list(erd1up, erd2up) 

    mm <- mismatchKernel(k=2, m=1, normalized=FALSE, presence=TRUE)
    
    testGetExplicitRepAll(dc, mm, expectedResPresenceUnnormalized, silent)

    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - normalized"))
        
    mm <- mismatchKernel(k=2, m=1, normalized=TRUE, presence=TRUE)
    
    testGetExplicitRepAll(dc, mm, expectedResPresenceUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Representation for feature",
                               "subset - unnormalized"))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
#    features <- list(c("AC","CA","GT"),c("AC","CA","GU"),c("AC","CA","GT"))
    features <- list(c("AC","CA","GT"),c("AC","CA","GT"))
    
#    expectedResFeaturesUnnormalized <- list(erd1u[,features[[1]]], 
#                                            erd2u[,features[[3]]]) 

    testGetExplicitRepAll(dc, mm, expectedResUnnormalized, silent,
                          features=features, colSubset=list(c(2,5,12),c(2,5,12),
                          c(2,5,12),c(2,22,112),c(2,22,122),c(2,22,122)))
    
    
    testStepName(silent, paste("Explict Representation for feature",
                               "subset - normalized"))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
#    features <- list(c("AC","CA","GT"),c("AC","CA","GU"),c("AC","CC","GT"))
    features <- list(c("AC","CA","GT"),c("AC","CA","GT"))

#    expectedResFeaturesUnnormalized <- list(erd1u[,features[[1]]], 
#                                            erd2u[,features[[3]]]) 
    
    testGetExplicitRepAll(dc, mm, expectedResUnnormalized, silent,
                          features=features, colSubset=list(c(2,5,12),c(2,5,12),
                          c(2,5,12),c(2,22,112),c(2,22,122),c(2,22,122)))
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))

    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)

    testERAgainstKernelMatrixAll(dc, mm, silent)

    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))

    mm <- mismatchKernel(k=2, m=1, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, mm, silent)
    
    testStepName(silent, paste("Explict Rep Sparse Subsetting against",
                               "Explicit Rep Dense"))
    
    mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
    
    testSubsettingOfSparseERAgainstDenseERAll(dc, mm, silent)
}

testKernelMatrix <- function(silent)
{
    testPosIndepKernelKernelMatrix(silent)
}

test_MismatchKernel <- function(seed=789, silent=TRUE)
{
    set.seed(seed)
    testKernelMatrix(silent)
    testExplicitRepresentation(silent)
}

## test is started by runit command in kebabs: runUnitTests.R
