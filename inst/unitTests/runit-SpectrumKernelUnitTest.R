##
## Unit Test via RUnit for spectrum kernel
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
    
    if (!is.null(annCharset))
    {
        dnaColNames <- paste(rep(dnaColNames, each=(nchar(annCharset))^2), 
                             paste(rep(unlist(strsplit(annCharset, split="")),
                                       each=nchar(annCharset)), 
                                   rep(unlist(strsplit(annCharset, split="")),
                                       nchar(annCharset)), sep=""), sep="")
        rnaColNames <- paste(rep(rnaColNames, each=(nchar(annCharset))^2), 
                             paste(rep(unlist(strsplit(annCharset, split="")),
                                       each=nchar(annCharset)), 
                                   rep(unlist(strsplit(annCharset, split="")),
                                       nchar(annCharset)), sep=""), sep="")
        aaColNames <- paste(rep(aaColNames, each=(nchar(annCharset))^2), 
                             paste(rep(unlist(strsplit(annCharset, split="")),
                                       each=nchar(annCharset)), 
                                   rep(unlist(strsplit(annCharset, split="")),
                                       nchar(annCharset)), sep=""), sep="")
    }
    
    testGetExplicitRepU(xc=xc, kernel=kernel, silent=silent, 
                        expectedResult=expectedResult, sel=sel,
                        features=features, colSubset=colSubset, 
                        annCharset=annCharset, ann1=ann1, ann2=ann2,
                        dnaColNames=dnaColNames, rnaColNames=rnaColNames,
                        aaColNames=aaColNames)
}

testPosIndepKernelKernelMatrix <- function(silent)
{
    testName(silent, "Spectrum Kernel Position Independent")
    
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
    
    symres1u <- defMat(matrix(c(18,9,11,9,13,14,11,14,23), nrow=3), 
                              namesX, namesX)
    symres2u <- defMat(matrix(c(9,4,4,11), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(11,7,8,7,8,9), nrow=3), namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized x - y",
                               "with subsetting"))
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), speck,
                           expectedResUnnormalized, silent, 
                           subsetX=c(3:5), subsetY=c(3:4))
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - y")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, normalized=TRUE)
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - z")

    symres1u <- defMat(matrix(c(18,9,11,9,13,14,11,14,23), nrow=3), 
                              namesX, namesX)
    symres2u <- defMat(matrix(c(32), nrow=1), namesZ, namesZ)
    asymresu <- defMat(matrix(c(13,10,18), nrow=3), namesX, namesZ)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    
    testGetKernelMatrixAll(xc, zc, speck, expectedResUnnormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - z")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, zc, speck, expectedResNormalized, silent)
    

    testStepName(silent, "Pos Indep - kmer presence unnormalized")
    
    symres1u <- defMat(matrix(c(9,5,6,5,7,7,6,7,11), nrow=3), namesX, namesX)
    symres2u <- defMat(matrix(c(6,3,3,8), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(6,4,5,5,4,6), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, normalized=FALSE, presence=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)
    

    testStepName(silent, "Pos Indep - kmer presence normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, normalized=TRUE, presence=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    
    
    testStepName(silent, 
                 "Pos Indep - kmer counts unnormalized tanslate lower to upper")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("GGaacGCcCGctT", "cCGAtacCTA", "GACcGAgTTATactct")
    names(xc) <- namesX
    yc <- c("aCCcgcTT", "ccggtTAaTA")
    names(yc) <- namesY
    zc <- "ACACgGAGaagACCATT"
    names(zc) <- namesZ
    
    symres1u <- defMat(matrix(c(18,9,11,9,13,14,11,14,23), nrow=3), 
                              namesX, namesX)
    symres2u <- defMat(matrix(c(9,4,4,11), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(11,7,8,7,8,9), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, normalized=FALSE, ignoreLower=FALSE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)

    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "ignore lowercase"))
    
    symres1u <- defMat(matrix(c(3,1,0,1,4,3,0,3,8), nrow=3), namesX, namesX)
    symres2u <- defMat(matrix(c(2,0,0,4), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(0,0,1,0,2,2), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, normalized=FALSE, ignoreLower=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent,
                           biovectorsOnly=TRUE )
}

testPosSpecKernelMatrix <- function(silent)
{
    testName(silent, "Spectrum Kernel Position Specific")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AACCCA", "AAACCA", "ACCCAA")
    names(xc) <- namesX
    yc <- c("CCCCCA", "ACCCAA")
    names(yc) <- namesY
    
    
    testStepName(silent, "Pos Spec - unnormalized")
    
    symres1u <- defMat(matrix(c(5,3,1,3,5,0,1,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(5,2,2,5), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(3,2,2,1,0,5), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, distWeight=1, normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)

    
    testStepName(silent, "Pos Spec - unnormalized with subsetting")
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    speck <- spectrumKernel(k=2, distWeight=1, normalized=FALSE)

    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), speck,
                           expectedResUnnormalized, silent, subsetX=c(3:5),
                           subsetY=c(3:4))
    
    
    testStepName(silent, "Pos Spec - normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    speck <- spectrumKernel(k=2, distWeight=1, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    
    testStepName(silent, "Pos Spec - attr pos on x")

    symres1u <- defMat(matrix(c(5,3,4,3,5,2,4,2,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(5,2,2,5), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(2,1,2,4,2,5), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, distWeight=1, normalized=FALSE)
    
    posX <- c(1,2,0)
    posY <- NULL
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, 
                           silent=silent, posX=posX, posY=posY)
    
    testStepName(silent, "Pos Spec - attr pos on y")

    symres1u <- defMat(matrix(c(5,3,1,3,5,0,1,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(5,2,2,5), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,0,3,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, distWeight=1, normalized=FALSE)
    
    posX <- NULL
    posY <- c(1,2)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, 
                           silent=silent, posX=posX, posY=posY)
    
    testStepName(silent, "Pos Spec - attr pos on x and y")
    
    symres1u <- defMat(matrix(c(5,3,4,3,5,2,4,2,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(5,2,2,5), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(3,1,3,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, distWeight=1, normalized=FALSE)
    
    posX <- c(1,2,0)
    posY <- c(1,2)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, 
                           silent=silent, posX=posX, posY=posY)
}

testAnnSpecKernelMatrix <- function(silent)
{
    testName(silent, "Spectrum Kernel Annotation Specific")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    
    xc <- c("GGAACGCCCGCTT", "CCGATACCTA", "GACCGAGTTATACTCT")
    names(xc) <- namesX
    yc <- c("ACCCGCTT", "CCGGTTAATA")
    names(yc) <- namesY
    
    
    testStepName(silent, "Ann Spec - unnormalized - annotation 1")
    
    symres1u <- defMat(matrix(c(14,5,8,5,11,10,8,10,19), nrow=3), 
                              namesX, namesX)
    symres2u <- defMat(matrix(c(7,2,2,11), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(6,2,4,5,6,4), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    annX1 <- c("eeeeeeeeeiiii", "eeeeiiiiii", "eeeeeeeeeeiiiiii")
    annY1 <- c("eeeiiiii", "eeeeiiiiii")
    annCharset1 <- "ei"
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent,
                           annCharset=annCharset1, annX=annX1, annY=annY1)
    
    
    testStepName(silent, "Ann Spec - normalized - annotation 1")
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent,
                           annCharset=annCharset1, annX=annX1, annY=annY1)
    
    testStepName(silent, "Ann Spec - normalized - annotation 1 with subsetting")
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    annX2 <- c("eeeiii", "iiieee")
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), speck,
                        expectedResUnnormalized, silent, 
                        annCharset=annCharset1, annX=c(annX2,annX1,annX2), 
                        annY=c(annX2,annY1,annX2), subsetX=c(3:5), 
                        subsetY=c(3:4))

    
    testStepName(silent, "Ann Spec - unnormalized - annotation 2")
    
    symres1u <- defMat(matrix(c(18,4,3,4,11,10,3,10,19), nrow=3), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(7,2,2,11), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(6,2,4,3,6,4), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    annX2 <- c("eeiiiiiiiiiii", "eeeeiiiiii", "eeeeeeeeeeiiiiii")
    annY2 <- c("eeeiiiii", "eeeeiiiiii")
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent,
                           annCharset=annCharset1, annX=annX2, annY=annY2)
    
    
    testStepName(silent, "Ann Spec - normalized - annotation 2")
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent,
                           annCharset=annCharset1, annX=annX2, annY=annY2)
    
}

testDistWeightKernelMatrix <- function(silent)
{
    testName(silent, "Spectrum Kernel with Distance Weights")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AACCCA", "AAACCA", "ACCCAA")
    names(xc) <- namesX
    yc <- c("CCCCCA", "ACCCAA")
    names(yc) <- namesY


    testStepName(silent, "Dist Weight - Explicit weight vector")
    
    speck <- spectrumKernel(k=2, distWeight=c(1,1,1,1,1), normalized=FALSE)
    speck1 <- spectrumKernel(k=2, normalized=FALSE)
    
    compareEqualKernelsAll(xc, yc, speck, speck1, silent)
    
    
    testStepName(silent, "Dist Weight - Predefined weight functions")
    
    symres1u <- defMat(matrix(c(6.6, 5.4, 5.0,
                                5.4, 6.6, 3.4,
                                5.0, 3.4, 6.6), nrow=3, byrow=TRUE), 
                              namesX, namesX)
    
    symres2u <- defMat(matrix(c(13.0,7.2,7.2,6.6), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(7.0,3.8,7.2,5.0,3.4,6.6), nrow=3), 
                              namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    speck <- spectrumKernel(k=2, distWeight=linWeight(sigma=5), 
                            normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    speck <- spectrumKernel(k=2, distWeight=linWeight(sigma=5), 
                            normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    
    symres1u <- defMat(matrix(c(6.637462, 5.456192, 5.394572,
                                5.456192, 6.637462, 3.976242,
                                5.394572, 3.976242, 6.637462), nrow=3,
                              byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(13.691288,7.434294,7.434294,6.637462), nrow=2),
                              namesY, namesY)
    
    asymresu <- defMat(matrix(c(7.345644, 5.394572,
                                4.037862, 3.976242,
                                7.434294, 6.637462), nrow=3, byrow=TRUE),
                              namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    speck <- spectrumKernel(k=2, distWeight=expWeight(sigma=5), 
                            normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, distWeight=expWeight(sigma=5), 
                            normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    symres1u <- defMat(matrix(c(6.921579, 5.882368, 6.222594,
                                5.882368, 6.921579, 4.850835,
                                6.222594, 4.850835, 6.921579), nrow=3,
                              byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(15.568664,8.508235,8.508235,6.921579), nrow=2),
                              namesY, namesY)
    
    asymresu <- defMat(matrix(c(8.284332, 6.222594,
                                4.510610, 4.850835,
                                8.508235, 6.921579), nrow=3, byrow=TRUE),
                              namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    speck <- spectrumKernel(k=2, distWeight=gaussWeight(sigma=5),
                            normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResUnnormalized, silent)
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    speck <- spectrumKernel(k=2, distWeight=gaussWeight(sigma=5),
                            normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, speck, expectedResNormalized, silent)
    
    
    testStepName(silent, "Dist Weight - User defined weight functions")
    
    speck <- spectrumKernel(k=2, 
                        distWeight=function(d, sigma=5){exp(-abs(d)/sigma)}, 
                        normalized=FALSE)
    speck1 <- spectrumKernel(k=2, distWeight=expWeight(sigma=5),
                             normalized=FALSE)
    
    compareEqualKernelsAll(xc, yc, speck, speck1, silent)
}

testExplicitRepresentation <- function(silent)
{
    ## UniProt IDs
    namesX <- c("P53554", "Q50EK3", "P15150", "Q64408", "P15538", "P97720", 
                "Q29527", "Q29552", "P15393", "P51663")
    
    ## first 10 AAs of Cytochrome P450
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
    
    testName(silent, "Explicit Representation for Spectrum Kernel")
    
    testStepName(silent, 
        "Explict Representation with/without zero features - unnormalized")

    namesY <- c("AA","AC","AG","AT","CA","CG","CT","GA","GC","GG","GT",
                "TA","TG","TT")
    
    erd1u <- defMat(matrix(c(3, 2, 6, 2, 7, 0, 0, 1, 6, 0, 0, 1, 1, 0,
                             4, 4, 0, 3, 2, 1, 1, 1, 0, 1, 3, 4, 3, 2,
                             5, 1, 4, 2, 3, 0, 1, 2, 3, 3, 1, 2, 2, 0,
                             5, 2, 3, 2, 1, 2, 1, 3, 2, 1, 1, 3, 1, 2,
                             5, 1, 5, 2, 3, 0, 1, 2, 4, 1, 1, 2, 2, 0,
                             4, 2, 4, 2, 3, 1, 1, 2, 3, 2, 1, 2, 2, 0,
                             5, 1, 5, 2, 3, 0, 1, 2, 4, 1, 1, 2, 2, 0,
                             5, 0, 3, 4, 4, 0, 0, 1, 4, 4, 0, 1, 3, 0,
                             2, 3, 4, 2, 3, 1, 1, 2, 2, 2, 2, 3, 2, 0,
                             4, 1, 3, 3, 3, 0, 1, 1, 3, 4, 1, 2, 3, 0), 
                            nrow=10, byrow=TRUE), namesX, namesY)
    
    namesY2 <- c("AD","AE","AF","AI","AK","AL","AR","AS","AW","DV","EA",
                 "EV","FR","FV","IA","IL","IW","KA","KS","LK","LR","LT",
                 "LW","MA","MD","MF","MT","NI","RA","RL","RV","SD","SS",
                 "ST","TA","TI","TM","VC","VN","VR","VT","VW","WA")
    
    erd2u <- defMat(matrix(
        c(0,0,0,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2,
          1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,0,
          0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,
          0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,
          0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,
          0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,
          0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,
          0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,
          0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,
          0,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,
          0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,
          0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,
          1), nrow=10, byrow=TRUE), namesX, namesY2)
    
    expectedResUnnormalized <- list(erd1u, erd2u) 
    
    
    speck <- spectrumKernel(k=2, normalized=FALSE)

    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent)

    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), speck, expectedResUnnormalized, 
                          silent, sel=c(3:12))
    

    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - unnormalized"))
    
    erd1up <- erd1u
    erd1up[which(erd1up > 0)] <- 1
    erd2up <- erd2u
    erd2up[which(erd2up > 0)] <- 1
    expectedResPresenceUnnormalized <- list(erd1up, erd2up)
    
    speck <- spectrumKernel(k=2, normalized=FALSE, presence=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResPresenceUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - normalized"))
    
    speck <- spectrumKernel(k=2, normalized=TRUE, presence=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResPresenceUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    speck <- spectrumKernel(k=2, normalized=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Representation for feature subset -",
                               "unnormalized"))
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent, 
                          features=list(c("AC","CC","CG"),c("AK","AM","DV")), 
                          colSubset=list(c(2,6),c(2,6,7),c(2,6,7),c(5,10),
                                         c(9,11,61),c(9,11,61)))
    
    testStepName(silent, paste("Explict Representation for feature subset -",
                               "normalized"))
    
    speck <- spectrumKernel(k=2, normalized=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent,
                          features=list(c("AC","CC","CG"),c("AK","AM","DV")), 
                          colSubset=list(c(2,6),c(2,6,7),c(2,6,7),c(5,10),
                                         c(9,11,61),c(9,11,61)))
    
    
    testStepName(silent, paste("Explict Rep with/without zero features against",
                               "kernel matrix - unnormalized"))

    speck <- spectrumKernel(k=2, normalized=FALSE)

    testERAgainstKernelMatrixAll(dc, speck, silent)

    
    testStepName(silent, paste("Explict Rep with/without zero features against",
                               "kernel matrix - normalized"))

    speck <- spectrumKernel(k=2, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, speck, silent)
    
    testStepName(silent, paste("Explict Rep Sparse Subsetting against Explicit",
                               "Rep Dense"))
    
    speck <- spectrumKernel(k=2, normalized=FALSE)
    
    testSubsettingOfSparseERAgainstDenseERAll(dc, speck, silent)
}

testExplicitRepAnnotSpecific <- function(silent)
{
    ## UniProt IDs
    namesX <- c("P53554", "Q50EK3", "P15150", "Q64408", "P15538", "P97720", 
                "Q29527", "Q29552", "P15393", "P51663")
    
    ## first 10 AAs of Cytochrome P450
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
    
    ## double sequences for each annotation letter
    dc <- paste(dc, dc, sep="")
    ac <- paste(ac, ac, sep="")
    names(ac) <- namesX
    names(dc) <- namesX
    
    dna <- DNAStringSet(dc)
    aa <- AAStringSet(ac)
    taa <- translate(dna)
    
    if (!all(taa == aa))
        stop("Error in DNA -> AA translation in Biostrings")
    
    testName(silent, paste("Explicit Representation Annotation Specific for",
                           "Spectrum Kernel"))
    
    testStepName(silent, paste("Explict Rep Ann Spec with/without zero",
                               "features - unnormalized"))
    
    erd1u <- defMat(matrix(
         c(3,0,3,2,2,6,6,2,2,7,1,7,0,0,0,0,1,0,1,6,6,0,0,0,0,1,1,1,1,0,0,    
           4,1,4,4,4,0,0,3,3,2,0,2,1,1,1,1,1,0,1,0,0,1,1,3,3,4,4,3,3,2,2,
           5,1,5,1,1,4,4,2,2,3,0,3,0,0,1,1,2,0,2,3,3,3,3,1,1,2,2,2,2,0,0,
           5,1,5,2,2,3,3,2,2,1,0,1,2,2,1,1,3,0,3,2,2,1,1,1,1,3,3,1,1,2,2,
           5,0,5,1,1,5,5,2,2,3,1,3,0,0,1,1,2,0,2,4,4,1,1,1,1,2,2,2,2,0,0,
           4,0,4,2,2,4,4,2,2,3,0,3,1,1,1,1,2,1,2,3,3,2,2,1,1,2,2,2,2,0,0,
           5,0,5,1,1,5,5,2,2,3,1,3,0,0,1,1,2,0,2,4,4,1,1,1,1,2,2,2,2,0,0,
           5,0,5,0,0,3,3,4,4,4,0,4,0,0,0,0,1,1,1,4,4,4,4,0,0,1,1,3,3,0,0,
           2,0,2,3,3,4,4,2,2,3,0,3,1,1,1,1,2,1,2,2,2,2,2,2,2,3,3,2,2,0,0,
           4,0,4,1,1,3,3,3,3,3,0,3,0,0,1,1,1,1,1,3,3,4,4,1,1,2,2,3,3,0,0),
           nrow=10, byrow=TRUE), namesX)

    rownames(erd1u) <- namesX
    colnames(erd1u) <- 
         c("AAaa","AAab","AAbb","ACaa","ACbb","AGaa","AGbb","ATaa","ATbb",
           "CAaa","CAab","CAbb","CGaa","CGbb","CTaa","CTbb","GAaa","GAab",
           "GAbb","GCaa","GCbb","GGaa","GGbb","GTaa","GTbb","TAaa","TAbb",
           "TGaa","TGbb","TTaa","TTbb")
    
    erd2u <- defMat(matrix(
         c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,
           0,0,0,0,1,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,1,
           0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,1,
           0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
           1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
           1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,
           0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,
           1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
           1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,0,0,1,1,1,1,
           0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,
           0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,0,
           0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,
           0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
           0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
           1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
           0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,
           0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,1,1,1,1,1), nrow=10, byrow=TRUE), namesX)

    rownames(erd2u) <- namesX
    colnames(erd2u) <- 
         c("ADaa","ADbb","AEaa","AEbb","AFaa","AFbb","AIaa","AIbb","AKaa",
           "AKbb","ALaa","ALbb","ARaa","ARbb","ASaa","ASbb","AWaa","AWbb",
           "CMab","DVaa","DVbb","EAaa","EAbb","EVaa","EVbb","FRaa","FRbb",
           "FVaa","FVbb","IAaa","IAbb","ILaa","ILbb","IWaa","IWbb","KAaa",
           "KAbb","KSaa","KSbb","LKaa","LKbb","LRaa","LRbb","LTaa","LTbb",
           "LWaa","LWbb","MAaa","MAbb","MDaa","MDbb","MFaa","MFbb","MTaa",
           "MTbb","NIaa","NIbb","RAaa","RAbb","RLaa","RLbb","RMab","RVaa",
           "RVbb","SDaa","SDbb","SMab","SSaa","SSbb","STaa","STbb","TAaa",
           "TAbb","TIaa","TIbb","TMaa","TMbb","VCaa","VCbb","VMab","VNaa",
           "VNbb","VRaa","VRbb","VTaa","VTbb","VWaa","VWbb","WAaa","WAbb",
           "WMab")
    
    expectedResUnnormalized <- list(erd1u, erd2u) 
    
    annCharset <- "ab"
    ann1 <-
         c(rep("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", 
               length(dna)))
    ann2 <- c(rep("aaaaaaaaaabbbbbbbbbb", length(aa)))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), speck, expectedResUnnormalized,
                    silent, sel=c(3:12), annCharset=annCharset, 
                    ann1=c(c("eeeiii","iiieee"), ann1, c("eeeiii","iiieee")), 
                    ann2=c(c("ei","ie"), ann2, c("ei","ie")))
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - unnormalized"))
    
    erd1up <- erd1u
    erd1up[which(erd1up > 0)] <- 1
    erd2up <- erd2u
    erd2up[which(erd2up > 0)] <- 1
    expectedResPresenceUnnormalized <- list(erd1up, erd2up)
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE, presence=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResPresenceUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features, presence - normalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE, presence=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResPresenceUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    testStepName(silent, paste("Explict Representation for feature",
                               "subset - unnormalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2, 
                          features=list(c("AAab","AAba","ACbb"),
                                        c("ADbb","AEba","AFaa")), 
                          colSubset=list(c(2,5),c(2,3,8),c(2,3,8),c(2,5),
                                         c(12,15,17),c(12,15,17)))
    
    
    testStepName(silent, paste("Explict Representation for feature",
                               "subset - normalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE)
    
    testGetExplicitRepAll(dc, speck, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2, 
                          features=list(c("AAab","AAba","ACbb"),
                                        c("ADbb","AEba","AFaa")), 
                          colSubset=list(c(2,5),c(2,3,8),c(2,3,8),c(2,5),
                                         c(12,15,17),c(12,15,17)))
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=FALSE)
    
    testERAgainstKernelMatrixAll(dc, speck, silent, annCharset=annCharset,
                                 ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))
    
    speck <- spectrumKernel(k=2, annSpec=TRUE, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, speck, silent, annCharset=annCharset,
                                 ann1=ann1, ann2=ann2)
}

testTooShortSequences <- function(seed="789", silent=TRUE)
{
    set.seed(seed)
    
    x <- DNAStringSet(c("T", "TCAC", "ACG","AAAGT","GC"))
    x
    names(x) <- paste("S", 1:length(x), sep="")
    
    testName(silent, "Test sequences shorter than k for Spectrum Kernel")
    
    testStepName(silent, "Kernel with k=3 - unnormalized")
    
    expectedResult1 <- matrix(c(0,0,0,0,0,0,
                                0,0,0,0,1,1,
                                0,0,1,0,0,0,
                                1,1,0,1,0,0,
                                0,0,0,0,0,0), nrow=5, byrow=TRUE)
    
    rownames(expectedResult1) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult1) <- c("AAA","AAG","ACG","AGT","CAC","TCA") 
    spec_3 <- spectrumKernel(k=3, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, spec_3, expRes=expectedResult1,
                                            silent=silent)
    
    annotationMetadata(x, annCharset=c("ei")) <- 
                                         c("i","eeii","eii","eeiii","ei")
    
    annotationCharset(x)
    annotationMetadata(x)
    
    testStepName(silent, "Annotation specific Kernel with k=3 - unnormalized")
    
    expectedResult2 <- matrix(c(0,0,0,0,0,0,
                                0,0,0,0,1,1,
                                0,0,1,0,0,0,
                                1,1,0,1,0,0,
                                0,0,0,0,0,0), nrow=5, byrow=TRUE)
    
    rownames(expectedResult2) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult2) <- c("AAAeei","AAGeii","ACGeii","AGTiii",
                                   "CACeii","TCAeei") 
    
    spec_3a <- spectrumKernel(k=3, annSpec=TRUE, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, spec_3a, 
                                     expRes=expectedResult2, silent=silent)
}

testKernelMatrix <- function(silent)
{
    testPosIndepKernelKernelMatrix(silent)
    testPosSpecKernelMatrix(silent)
    testAnnSpecKernelMatrix(silent)
    testDistWeightKernelMatrix(silent)
}

test_SpectrumKernel <- function(seed=789, silent=TRUE)
{
    set.seed(seed)
    testKernelMatrix(silent)
    testExplicitRepresentation(silent)
    testExplicitRepAnnotSpecific(silent)
    testTooShortSequences(silent=silent)
}
