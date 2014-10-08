##
## Unit Test via RUnit for gappy pair kernel
##

library(RUnit)

testGetExplicitRepAll <- function(xc, kernel, expectedResult, silent, sel,
                                  features=NULL, colSubset=NULL, 
                                  annCharset=NULL, ann1=NULL, ann2=NULL)
{
    dnaColNames <- as.character(rbind(t(as.matrix(
                       paste(rep(c("A","C","G","T"),each=4),
                             rep(c("A","C","G","T"),4), sep=""))),
                       paste(rep(c("A","C","G","T"),each=4), 
                             rep(".", 16), rep(c("A","C","G","T"),4), sep="")))
    rnaColNames <- as.character(rbind(t(as.matrix(
                             paste(rep(c("A","C","G","U"),each=4),
                             rep(c("A","C","G","U"),4), sep=""))),
                       paste(rep(c("A","C","G","U"),each=4), 
                             rep(".", 16), rep(c("A","C","G","U"),4), sep="")))
    aaColNames <- as.character(rbind(t(as.matrix(
        paste(rep(c("A","C","D","E","F","G","H","I","K","L","M","N","P",
                    "Q","R","S","T","U", "V","W","Y"),each=21),
              rep(c("A","C","D","E","F","G","H","I","K","L","M","N","P",
                    "Q","R","S","T","U", "V","W","Y"),21), sep=""))), 
        paste(rep(c("A","C","D","E","F","G","H","I","K","L","M","N","P",
                    "Q","R","S","T","U", "V","W","Y"),each=21), rep(".", 441),
              rep(c("A","C","D","E","F","G","H","I","K","L","M","N","P",
                    "Q","R","S","T","U", "V","W","Y"),21), sep="")))
    
    if (!is.null(annCharset))
    {
        annotNames <- as.character(rbind(t(as.matrix(
            paste(rep(unlist(strsplit(annCharset, split="")), 
                      each=nchar(annCharset)), 
                  rep(unlist(strsplit(annCharset, split="")), 
                      nchar(annCharset)), sep=""))), 
            paste(rep(unlist(strsplit(annCharset, split="")), 
                      each=nchar(annCharset)), 
                  rep(".", (nchar(annCharset))^2),
                  rep(unlist(strsplit(annCharset, split="")), 
                      nchar(annCharset)), sep="")))

        dnaColNames <- 
            as.character(matrix(dnaColNames, nrow=2)[,
                rep(1:(length(dnaColNames)/2), each=nchar(annCharset)^2)])
        dnaAnnotNames <- rep(annotNames, length(dnaColNames)/2)
        dnaColNames <- sapply(1:length(dnaColNames), 
            function(i){paste(dnaColNames[i], dnaAnnotNames[i], sep="")})
        rnaColNames <- as.character(matrix(rnaColNames, nrow=2)[,
            rep(1:(length(rnaColNames)/2), each=nchar(annCharset)^2)])
        rnaColNames <- sapply(1:length(rnaColNames), 
            function(i){paste(rnaColNames[i], dnaAnnotNames[i], sep="")})
        aaAnnotNames <- rep(annotNames, length(aaColNames)/2)
        aaColNames <- as.character(matrix(aaColNames, nrow=2)[,
            rep(1:(length(aaColNames)/2), each=nchar(annCharset)^2)])
        aaColNames <- sapply(1:length(aaColNames), 
            function(i){paste(aaColNames[i], aaAnnotNames[i], sep="")})
    }

    testGetExplicitRepU(xc=xc, kernel=kernel, silent=silent, 
                        expectedResult=expectedResult, sel=sel,
                        features=features, colSubset=colSubset, 
                        annCharset=annCharset, ann1=ann1, ann2=ann2,
                        dnaColNames=dnaColNames, rnaColNames=rnaColNames,
                        aaColNames=aaColNames)
}

testPosIndepKernelMatrix <- function(silent)
{
    testName(silent, "Gappy Pair Kernel Position Independent")
    
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
    
    symres1u <- defMat(matrix(c(27,0,1,0,18,2,1,2,38), nrow=3), namesX, namesX)

    symres2u <- defMat(matrix(c(12,2,2,18), nrow=2), namesY, namesY)

    asymresu <- defMat(matrix(c(9,0,0,2,0,1), nrow=3), namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE)
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized x - y",
                               "with subsetting"))
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), gappy,
                           expectedResUnnormalized, silent, 
                           subsetX=c(3:5), subsetY=c(3:4))
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - y")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=TRUE)

    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - z")

    symres1u <- defMat(matrix(c(27,0,1,0,18,2,1,2,38), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(39), nrow=1), namesZ, namesZ)
    
    asymresu <- defMat(matrix(c(1,0,3), nrow=3), namesX, namesZ)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE)
    
    testGetKernelMatrixAll(xc, zc, gappy, expectedResUnnormalized, silent)
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - z")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, zc, gappy, expectedResNormalized, silent)
    

    testStepName(silent, "Pos Indep - kmer presence unnormalized")
    
    symres1u <- defMat(matrix(c(27,0,1,0,18,2,1,2,35), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(12,2,2,18), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(9,0,0,2,0,1), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE, presence=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)
    

    testStepName(silent, "Pos Indep - kmer presence normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=TRUE, presence=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized translate",
                               "lower to upper"))
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("GGaacGCcCGctT", "cCGAtacCTA", "GACcGAgTTATactct")
    names(xc) <- namesX
    yc <- c("aCCcgcTT", "ccggtTAaTA")
    names(yc) <- namesY
    zc <- "ACACgGAGaagACCATT"
    names(zc) <- namesZ
    
    symres1u <- defMat(matrix(c(27,0,1,0,18,2,1,2,38), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(12,2,2,18), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(9,0,0,2,0,1), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE, ignoreLower=FALSE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)

    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "ignore lowercase"))
    
    symres1u <- defMat(matrix(c(0,0,0,0,0,0,0,0,1), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(0,0,0,0), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(0,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=2, normalized=FALSE, ignoreLower=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent,
                           biovectorsOnly=TRUE )
}

testPosSpecKernelMatrix <- function(silent)
{
    testName(silent, "Gappy Pair Kernel Position Specific")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AACCCA", "AAACCA", "ACCCAA")
    names(xc) <- namesX
    yc <- c("AAACCACCA", "AACCCACC")
    names(yc) <- namesY
    
    
    testStepName(silent, "Pos Spec - unnormalized")
    
    symres1u <- defMat(matrix(c(5,1,0,1,5,0,0,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(11,4,4,9), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,5,0,5,1,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)

    
    testStepName(silent, "Pos Spec - unnormalized with subsetting")
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=FALSE)

    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), gappy,
                        expectedResUnnormalized, silent, 
                        subsetX=c(3:5), subsetY=c(3:4))
    
    
    testStepName(silent, "Pos Spec - normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    testStepName(silent, "Pos Spec - attr pos on x")

    symres1u <- defMat(matrix(c(5,1,3,1,5,0,3,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(11,4,4,9), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(0,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=FALSE)
    
    posX <- c(1,2,0)
    posY <- NULL
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized,
                           silent=silent, posX=posX, posY=posY)
    
    testStepName(silent, "Pos Spec - attr pos on y")

    symres1u <- defMat(matrix(c(5,1,0,1,5,0,0,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(11,0,0,9), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=FALSE)
    
    posX <- NULL
    posY <- c(1,2)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized,
                           silent=silent, posX=posX, posY=posY)
    
    testStepName(silent, "Pos Spec - attr pos on x and y")
    
    symres1u <- defMat(matrix(c(5,1,3,1,5,0,3,0,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(11,0,0,9), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,0,0,0,1,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=1, normalized=FALSE)
    
    posX <- c(1,2,0)
    posY <- c(1,2)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, 
                           silent=silent, posX=posX, posY=posY)
}

testAnnSpecKernelMatrix <- function(silent)
{
    testName(silent, "Gappy Pair Kernel Annotation Specific")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    
    xc <- c("GGAACGCCCGCTT", "CCGATACCTA", "GACCGAGTTATACTCT")
    names(xc) <- namesX
    yc <- c("ACCCGCTT", "CCGGTTAATA")
    names(yc) <- namesY
    
    
    testStepName(silent, "Ann Spec - unnormalized - annotation 1")
    
    symres1u <- defMat(matrix(c(19,0,1,0,13,2,1,2,25), nrow=3), 
                       namesX, namesX)    
    symres2u <- defMat(matrix(c(9,0,0,13), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(7,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=FALSE)
    
    annX1 <- c("eeeeeeeeeiiii", "eeeeiiiiii", "eeeeeeeeeeiiiiii")
    annY1 <- c("eeeeiiii", "eeeeiiiiii")
    annCharset1 <- "ei"
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent,
                           annCharset=annCharset1, annX=annX1, annY=annY1)
    
    
    testStepName(silent, "Ann Spec - normalized - annotation 1")
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent,
                           annCharset=annCharset1, annX=annX1, annY=annY1)
    
    testStepName(silent, "Ann Spec - normalized - annotation 1 with subsetting")
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    annX2 <- c("eeeiii", "iiieee")
    
    gappy <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), gappy,
                           expectedResUnnormalized, silent,
                           annCharset=annCharset1, annX=c(annX2,annX1,annX2),
                           annY=c(annX2,annY1,annX2), subsetX=c(3:5), 
                           subsetY=c(3:4))

    
    testStepName(silent, "Ann Spec - unnormalized - annotation 2")
    
    symres1u <- defMat(matrix(c(19,0,0,0,13,2,0,2,25), nrow=3), namesX, namesX)
    symres2u <- defMat(matrix(c(9,0,0,13), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(3,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=FALSE)
    
    annX2 <- c("eeiiiiiiiiiii", "eeeeiiiiii", "eeeeeeeeeeiiiiii")
    annY2 <- c("eeeiiiii", "eeeeiiiiii")
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent,
                           annCharset=annCharset1, annX=annX2, annY=annY2)
    
    
    testStepName(silent, "Ann Spec - normalized - annotation 2")
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent,
                           annCharset=annCharset1, annX=annX2, annY=annY2)
    
}

testDistWeightKernelMatrix <- function(silent)
{
    testName(silent, "Gappy Pair Kernel with Distance Weights")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AACCCA", "AAACCA", "ACCCAA")
    names(xc) <- namesX
    yc <- c("AAACCACCA", "AACCCACC")
    names(yc) <- namesY
    

    testStepName(silent, "Dist Weight - Explicit weight vector")
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=c(1,1,1,1,1), 
                             normalized=FALSE)
    gappy1 <- gappyPairKernel(k=2, m=1, normalized=FALSE)
    
    compareEqualKernelsAll(xc, yc, gappy, gappy1, silent)
    
    
    testStepName(silent, "Dist Weight - Predefined weight functions")
    
    symres1u <- defMat(matrix(c(5.0, 1.8, 2.4,
                                1.8, 5.0, 0.0,
                                2.4, 0.0, 5.0), nrow=3, byrow=TRUE), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(11.8,4.8,4.8,9.0), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(1.8,5.4,0.0,5.0,1.8,2.4), nrow=3),
                       namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=linWeight(sigma=5),
                             normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    gappy <- gappyPairKernel(k=2, m=1, distWeight=linWeight(sigma=5),
                             normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    symres1u <- defMat(matrix(c(5.000000, 1.818731, 2.456192,
                                1.818731, 5.000000, 0.000000,
                                2.456192, 0.000000, 5.000000), nrow=3,
                              byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(12.097623, 4.818731, 4.818731, 9.000000),
                              nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1.818731, 5.000000,
                                5.548812, 1.818731,
                                0.000000, 2.456192), nrow=3, byrow=TRUE),
                       namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    gappy <- gappyPairKernel(k=2, m=1, distWeight=expWeight(sigma=5),
                             normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=expWeight(sigma=5),
                             normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    symres1u <- defMat(matrix(c(5.000000, 1.960789, 2.882368,
                                1.960789, 5.000000, 0.000000,
                                2.882368, 0.000000, 5.000000), nrow=3,
                              byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(12.395353,4.960789,4.960789,9.000000), nrow=2),
                       namesY, namesY)
    
    asymresu <- defMat(matrix(c(1.960789, 5.000000,
                                5.697676, 1.960789,
                                0.000000, 2.882368), nrow=3, byrow=TRUE),
                       namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    gappy <- gappyPairKernel(k=2, m=1, distWeight=gaussWeight(sigma=5),
                             normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResUnnormalized, silent)
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    gappy <- gappyPairKernel(k=2, m=1, distWeight=gaussWeight(sigma=5),
                             normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, gappy, expectedResNormalized, silent)
    
    
    testStepName(silent, "Dist Weight - User defined weight functions")
    
    gappy <- gappyPairKernel(k=2, m=1, 
        distWeight=function(d, sigma=5){exp(-abs(d)/sigma)}, normalized=FALSE)
    gappy1 <- gappyPairKernel(k=2, m=1, distWeight=expWeight(sigma=5),
        normalized=FALSE)
    
    compareEqualKernelsAll(xc, yc, gappy, gappy1, silent)
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
    
    testName(silent, "Explicit Representation for Gappy Pair Kernel")
    
    testStepName(silent, paste("Explict Representation with/without zero",
                 "features - unnormalized"))

    namesY <- c("AA","A.A","AC","A.C","AG","A.G","AT","A.T","CA","C.A","C.C",
                "CG","C.G","CT","C.T","GA","G.A","GC","G.C","GG","G.G","GT",
                "G.T","TA","T.A","T.C","TG","T.G","TT","T.T")
    
    erd1u <- defMat(matrix(
        c(3,3,2,6,6,3,2,1,7,3,1,0,3,0,0,1,5,6,1,0,0,0,0,1,1,0,1,1,0,0,
          4,4,4,2,0,3,3,2,2,2,0,1,0,1,2,1,3,0,1,1,0,3,1,4,2,1,3,2,2,3,
          5,4,1,1,4,5,2,2,3,3,1,0,0,1,0,2,4,3,2,3,2,1,0,2,1,0,2,2,0,1,
          5,5,2,1,3,4,2,2,1,2,0,2,0,1,2,3,2,2,3,1,1,1,0,3,3,0,1,2,2,1,
          5,4,1,2,5,5,2,2,3,2,1,0,1,1,0,2,5,4,1,1,1,1,0,2,1,1,2,1,0,1,
          4,4,2,2,4,5,2,1,3,2,1,1,1,1,1,2,4,3,2,2,1,1,0,2,1,0,2,2,0,1,
          5,4,1,2,5,5,2,2,3,2,1,0,1,1,0,2,5,4,1,1,1,1,0,2,1,1,2,1,0,1,
          5,4,0,2,3,5,4,1,4,2,0,0,1,0,1,1,5,4,2,4,1,0,0,1,0,0,3,3,0,1,
          2,3,3,2,4,4,2,2,3,1,1,1,2,1,1,2,4,2,2,2,1,2,0,3,2,0,2,2,0,1,
          4,3,1,1,3,5,3,2,3,3,1,0,0,1,0,1,4,3,2,4,2,1,0,2,0,0,3,3,0,2
          ), nrow=10, byrow=TRUE), namesX, namesY)

    namesY2 <- 
        c("A.A","AD","AE","AF","AI","AK","AL","AR","A.R","AS","A.S","A.V",
          "AW","A.W","D.N","D.R","DV","D.W","EA","E.C","EV","E.W","F.L","FR",
          "FV","IA","I.A","IL","I.S","I.T","IW","KA","K.D","K.E","K.R","KS",
          "L.A","LK","L.M","LR","L.S","LT","L.V","LW","MA","MD","MF","M.F",
          "M.I","M.L","MT","M.V","NI","N.L","RA","R.K","RL","R.R","R.T","RV",
          "R.W","S.A","SD","SS","ST","S.T","S.V","TA","T.A","T.D","T.F","TI",
          "TM","T.S","V.A","VC","V.I","VN","VR","VT","VW","WA","W.K")
    
    erd2u <- defMat(matrix(
        c(0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,
          0,1,0,2,1,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,
          0,1,0,0,1,0,0,0,1,1,0,0,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,1,0,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,
          1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,
          0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,
          0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
          0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,
          0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,
          1,1,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,
          0,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,
          0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,
          0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,1,0,0,0,0,0,0,0,2,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,
          0,0,1,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,1,1,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,
          0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,
          0,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
          0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1), 
          nrow=10, byrow=TRUE), namesX, namesY2)
    
    expectedResUnnormalized <- list(erd1u, erd2u) 
    gappy <- gappyPairKernel(k=1, m=1, normalized=FALSE)

    testGetExplicitRepAll(dc, gappy, expectedResUnnormalized, silent)

    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    gappy <- gappyPairKernel(k=1, m=1, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), gappy, expectedResUnnormalized, 
                       silent, sel=c(3:12))
    

    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    gappy <- gappyPairKernel(k=1, m=1, normalized=TRUE)
    
    testGetExplicitRepAll(dc, gappy, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))

    gappy <- gappyPairKernel(k=1, m=1, normalized=FALSE)

    testERAgainstKernelMatrixAll(dc, gappy, silent)

    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))

    gappy <- gappyPairKernel(k=1, m=1, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, gappy, silent)
    
    testStepName(silent, paste("Explict Rep Sparse Subsetting against",
                               "Explicit Rep Dense"))
    
    gappy <- gappyPairKernel(k=1, m=1, normalized=FALSE)
    
    testSubsettingOfSparseERAgainstDenseERAll(dc, gappy, silent)
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
    
    testName(silent, paste("Explicit Representation Annotation Specific",
                           "for Gappy Pair Kernel"))
    
    testStepName(silent, paste("Explict Rep Ann Spec with/without zero",
                               "features - unnormalized"))
    
    erd1u <- defMat(matrix(
        c(3,3,0,3,3,2,6,2,6,6,3,6,3,2,1,0,2,1,7,3,1,7,3,1,1,0,3,0,3,0,0,
          1,0,0,1,5,0,1,1,5,6,1,6,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,1,
          1,1,0,0,0,0,4,4,1,4,4,4,2,4,2,0,3,0,3,3,2,1,3,2,2,2,0,2,2,0,0,
          1,0,1,0,1,2,0,1,2,1,3,0,0,1,3,0,1,0,1,1,0,1,0,3,1,0,3,1,4,2,1,
          4,2,1,1,3,2,3,2,2,3,2,3,5,4,1,5,4,1,1,1,1,4,5,4,5,2,2,1,2,2,3,
          3,0,3,3,1,1,0,0,0,0,1,0,0,1,0,2,4,0,1,2,4,3,2,3,2,3,2,3,2,1,0,
          0,1,0,2,1,0,2,1,0,0,2,2,2,2,0,1,0,1,5,5,1,5,5,2,1,2,1,3,4,3,4,
          2,2,1,2,2,1,2,0,1,2,0,0,2,0,2,0,1,2,0,1,2,3,2,0,1,3,2,2,3,2,3,
          1,1,1,1,1,0,0,1,0,3,3,0,3,3,0,0,1,2,1,2,2,1,2,1,5,4,0,5,4,1,2,
          1,2,5,5,5,5,2,2,0,2,2,3,2,1,3,2,1,1,0,1,0,1,1,0,1,1,0,2,5,0,1,
          2,5,4,1,4,1,1,1,1,1,1,0,0,1,0,2,1,0,2,1,1,1,2,1,2,1,0,1,0,1,4,
          4,0,4,4,2,2,2,2,4,5,4,5,2,1,0,2,1,3,2,0,3,2,1,1,1,1,1,1,1,1,0,
          1,1,2,4,1,1,2,4,3,2,3,2,2,1,2,1,1,0,1,1,0,2,1,0,2,1,0,0,2,2,2,
          2,0,1,0,1,5,4,0,5,4,1,2,1,2,5,5,5,5,2,2,0,2,2,3,2,1,3,2,1,1,0,
          1,0,1,1,0,1,1,0,2,5,0,1,2,5,4,1,4,1,1,1,1,1,1,0,0,1,0,2,1,0,2,
          1,1,1,2,1,2,1,0,1,0,1,5,4,0,5,4,0,2,0,2,3,5,3,5,4,1,0,4,1,4,2,
          0,4,2,0,0,0,1,0,1,0,1,0,0,1,1,5,1,1,1,5,4,2,4,2,4,1,4,1,0,0,1,
          0,0,1,0,0,1,0,0,0,3,3,3,3,0,1,0,1,2,3,0,2,3,3,2,3,2,4,4,4,4,2,
          2,0,2,2,3,1,0,3,1,1,1,1,2,1,2,1,1,0,1,1,2,4,1,1,2,4,2,2,2,2,2,
          1,2,1,2,0,1,2,0,3,2,0,3,2,0,0,2,2,2,2,0,1,0,1,4,3,0,4,3,1,1,1,
          1,3,5,3,5,3,2,0,3,2,3,3,0,3,3,1,1,0,0,0,0,1,0,0,1,0,1,4,1,1,1,
          4,3,2,3,2,4,2,4,2,1,0,1,1,0,2,0,0,2,0,0,0,3,3,3,3,0,2,0,2
          ), nrow=10, byrow=TRUE))
    
    rownames(erd1u) <- namesX
    colnames(erd1u) <- 
        c("AAaa","A.Aa.a","AAab","AAbb","A.Ab.b","ACaa","A.Ca.a","ACbb",
          "A.Cb.b","AGaa","A.Ga.a","AGbb","A.Gb.b","ATaa","A.Ta.a","A.Ta.b",
          "ATbb","A.Tb.b","CAaa","C.Aa.a","CAab","CAbb","C.Ab.b","C.Ca.a",
          "C.Cb.b","CGaa","C.Ga.a","CGbb","C.Gb.b","CTaa","C.Ta.a","C.Ta.b",
          "CTbb","C.Tb.b","GAaa","G.Aa.a","GAab","G.Aa.b","GAbb","G.Ab.b",
          "GCaa","G.Ca.a","GCbb","G.Cb.b","GGaa","G.Ga.a","GGbb","G.Gb.b",
          "GTaa","G.Ta.a","G.Ta.b","GTbb","G.Tb.b","TAaa","T.Aa.a","T.Aa.b",
          "TAbb","T.Ab.b","T.Ca.a","T.Cb.b","TGaa","T.Ga.a","TGbb","T.Gb.b",
          "TTaa","T.Ta.a","TTbb","T.Tb.b")
    
    erd2u <- defMat(matrix(
        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,1,1,0,0,1,1,2,2,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,0,0,
          1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,0,1,
          0,0,0,0,0,0,0,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,
          1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,
          1,1,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,
          0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,
          0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
          0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,
          1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,
          0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,
          0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
          0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,
          1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,
          0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,
          1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,
          0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,
          0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
          1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,
          1,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,
          0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
          1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
          0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
          1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,1,1,0,0,1,1,1,
          1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,
          0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,
          0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,
          1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,
          1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,
          0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,1,
          1,1), nrow=10, byrow=TRUE))
    
    rownames(erd2u) <- namesX
    colnames(erd2u) <- 
        c("A.Aa.a","A.Ab.b","ADaa","ADbb","AEaa","AEbb","AFaa","AFbb","AIaa",
          "AIbb","AKaa","AKbb","ALaa","ALbb","A.Ma.b","ARaa","A.Ra.a","ARbb",
          "A.Rb.b","ASaa","A.Sa.a","ASbb","A.Sb.b","A.Va.a","A.Vb.b","AWaa",
          "A.Wa.a","AWbb","A.Wb.b","C.Aa.b","CMab","D.Na.a","D.Nb.b","D.Ra.a",
          "D.Rb.b","DVaa","DVbb","D.Wa.a","D.Wb.b","EAaa","EAbb","E.Ca.a",
          "E.Cb.b","EVaa","EVbb","E.Wa.a","E.Wb.b","F.La.a","F.Lb.b","F.Ma.b",
          "FRaa","FRbb","FVaa","FVbb","IAaa","I.Aa.a","IAbb","I.Ab.b","ILaa",
          "ILbb","I.Sa.a","I.Sb.b","I.Ta.a","I.Tb.b","IWaa","IWbb","KAaa",
          "KAbb","K.Da.a","K.Db.b","K.Ea.a","K.Eb.b","K.Ra.a","K.Rb.b","KSaa",
          "KSbb","L.Aa.a","L.Ab.b","LKaa","LKbb","L.Ma.a","L.Mb.b","LRaa",
          "LRbb","L.Sa.a","L.Sb.b","LTaa","LTbb","L.Va.a","L.Vb.b","LWaa",
          "LWbb","MAaa","MAbb","MDaa","MDbb","MFaa","M.Fa.a","MFbb","M.Fb.b",
          "M.Ia.a","M.Ib.b","M.La.a","M.Lb.b","MTaa","MTbb","M.Va.a","M.Vb.b",
          "NIaa","NIbb","N.La.a","N.Lb.b","RAaa","R.Aa.b","RAbb","R.Ka.a",
          "R.Kb.b","RLaa","RLbb","RMab","R.Ra.a","R.Rb.b","R.Ta.a","R.Tb.b",
          "RVaa","RVbb","R.Wa.a","R.Wb.b","S.Aa.a","S.Ab.b","SDaa","SDbb",
          "SMab","S.Ma.b","SSaa","SSbb","STaa","S.Ta.a","S.Ta.b","STbb",
          "S.Tb.b","S.Va.a","S.Vb.b","TAaa","T.Aa.a","TAbb","T.Ab.b","T.Da.a",
          "T.Db.b","T.Fa.a","T.Fb.b","TIaa","TIbb","TMaa","TMbb","T.Sa.a",
          "T.Sb.b","V.Aa.a","V.Ab.b","VCaa","VCbb","V.Da.b","V.Ia.a","V.Ib.b",
          "VMab","V.Ma.b","VNaa","VNbb","VRaa","VRbb","VTaa","VTbb","VWaa",
          "VWbb","WAaa","W.Aa.b","WAbb","W.Ka.a","W.Kb.b","WMab")
    
    expectedResUnnormalized <- list(erd1u, erd2u) 
    
    annCharset <- "ab"
    ann1 <- 
        c(rep("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", 
              length(dna)))
    ann2 <- c(rep("aaaaaaaaaabbbbbbbbbb", length(aa)))
    
    gappy <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(dc, gappy, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    gappy <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), gappy, expectedResUnnormalized, 
                          silent, sel=c(3:12), annCharset=annCharset, 
                          ann1=c(c("eeeiii","iiieee"), ann1, 
                                 c("eeeiii","iiieee")), 
                          ann2=c(c("ei","ie"), ann2, c("ei","ie")))
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    gappy <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=TRUE)
    
    testGetExplicitRepAll(dc, gappy, expectedResUnnormalized, silent,
                          annCharset=annCharset, ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))
    
    gappy <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=FALSE)
    
    testERAgainstKernelMatrixAll(dc, gappy, silent, annCharset=annCharset,
                                 ann1=ann1, ann2=ann2)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))
    
    gappy <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, gappy, silent, annCharset=annCharset,
                                 ann1=ann1, ann2=ann2)
}

testTooShortSequences <- function(seed="789", silent=TRUE)
{
    set.seed(seed)
    
    x <- DNAStringSet(c("ACCGAT", "ACGT","AAAGTA","GC"))
    x
    names(x) <- paste("S", 1:length(x), sep="")
    
    testName(silent, "Test sequences shorter than 2*k+m for Gappy Pair Kernel")
    
    testStepName(silent, "Kernel with k=1, m=1 - unnormalized")
    
    expectedResult1 <- matrix(
                              c(0,0,1,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,
                                0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,0,
                                2,1,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0,1,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0), 
                               nrow=4, byrow=TRUE)
    
    rownames(expectedResult1) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult1) <- 
        c("AA","A.A","AC","A.C","AG","A.G","AT","A.T","C.A","CC",
          "CG","C.G","C.T","GA","G.A","GC","GT","G.T","TA") 
    gappy1_1 <- gappyPairKernel(k=1, m=1, normalized=FALSE)

    testDenseAndSparseERAgainstKernelMatrix(x, gappy1_1,
                             expRes=expectedResult1, silent=silent)

    testStepName(silent, "Kernel with k=1, m=3 - unnormalized")

    expectedResult2 <- matrix(
        c(0,0,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0,
          0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,
          2,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0), 
         nrow=4, byrow=TRUE)
    
    rownames(expectedResult2) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult2) <- 
    c("AA","A.A","A..A","A...A","AC","A.C","AG","A.G","A..G",
      "AT","A.T","A..T","A...T","C.A","C..A","CC","CG","C.G",
      "C.T","C..T","C...T","GA","G.A","GC","GT","G.T","TA") 
    
    gappy1_3 <- gappyPairKernel(k=1, m=3, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, gappy1_3,
                                expRes=expectedResult2, silent=silent)
    
    testStepName(silent, "Kernel with k=2, m=1 - unnormalized")

    expectedResult3 <- matrix(
                              c(0,0,0,0,1,1,0,0,1,1,1,
                                0,0,0,0,0,0,1,0,0,0,0,
                                1,1,1,1,0,0,0,1,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0), nrow=4, byrow=TRUE)
    
    rownames(expectedResult3) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult3) <- 
        c("AAAG","AAGT","AA.GT","AA.TA","ACCG","AC.GA","ACGT","AGTA",
          "CC.AT","CCGA","CGAT") 
    
    gappy2_1 <- gappyPairKernel(k=2, m=1, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, gappy2_1,
                              expRes=expectedResult3, silent=silent)
    
    x <- DNAStringSet(paste(as.character(x), as.character(x), sep=""))
    names(x) <- paste("S", 1:length(x), sep="")
    annotationMetadata(x, annCharset=c("ei")) <- 
         c("eeeeeeiiiiii","eeeeiiii","eeeeeeiiiiii","eeii")
    
    annotationCharset(x)
    annotationMetadata(x)
    
    testStepName(silent, paste("Annotation specific Kernel with k=1,",
                               "m=1 - unnormalized"))
    
    expectedResult1 <- matrix(
        c(0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,0,1,1,0,0,1,0,
          0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,
          0,0,0,0,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,1,2,1,
          1,1,2,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
          1,0,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0), 
         nrow=4, byrow=TRUE)
    
    rownames(expectedResult1) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult1) <- 
        c("AAee","A.Ae.e","AAei","A.Ae.i","AAii","A.Ai.i","ACee","A.Ce.e",
          "ACii","A.Ci.i","AGee","A.Ge.e","AGii","A.Gi.i","ATee","A.Te.e",
          "ATii","A.Ti.i","C.Ae.e","C.Ai.i","CCee","C.Ce.i","CCii","CGee",
          "C.Ge.e","CGei","CGii","C.Gi.i","C.Te.e","C.Ti.i","GAee","G.Ae.e",
          "G.Ae.i","GAii","G.Ai.i","GCee","GCii","G.Ge.i","GTee","G.Te.e",
          "GTii","G.Ti.i","TAee","TAei","T.Ae.i","TAii","T.Ce.i")
    
    gappy1_1 <- gappyPairKernel(k=1, m=1, annSpec=TRUE, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, gappy1_1,
                                 expRes=expectedResult1, silent=silent)
    
    testStepName(silent, paste("Annotation specific Kernel with k=2,",
                               "m=1 - unnormalized"))
    
    expectedResult2 <- matrix(
        c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,
          1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
          0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,
          1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0), nrow=4, byrow=TRUE)
    
    rownames(expectedResult2) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult2) <- 
        c("AAAAeiii","AAAGeeee","AA.AGei.ii","AAAGiiii","AAGTeeee",
          "AA.GTee.ee","AAGTiiii","AA.GTii.ii","AA.TAee.ee","AA.TAii.ii",
          "ACCGeeee","ACCGiiii","AC.GAee.ee","AC.GAii.ii","ACGTeeee",
          "ACGTiiii","AC.TAee.ei","AG.AAee.ei","AGTAeeee","AGTAiiii",
          "ATACeeii","AT.CCee.ii","CC.ATee.ee","CC.ATii.ii","CCGAeeee",
          "CCGAiiii","CG.ACee.ii","CGATeeee","CGATiiii","CGTAeeei",
          "CG.TAee.ei","GA.ACee.ii","GATAeeei","GCGCeeii","GTAAeeei",
          "GT.AAee.ii","GTACeeii","GT.CGee.ii","TAAAeeii","TA.AAee.ii",
          "TACCeiii","TACGeiii","TA.CGei.ii","TA.GTei.ii") 
    
    gappy2_1 <- gappyPairKernel(k=2, m=1, annSpec=TRUE, normalized=FALSE)
    testDenseAndSparseERAgainstKernelMatrix(x, gappy2_1,
                                 expRes=expectedResult2, silent=silent)
    
    testStepName(silent, paste("Annotation specific Kernel with k=1,",
                               "m=3 - unnormalized"))
    
    expectedResult3 <- matrix(
        c(0,0,0,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,
          0,1,1,0,1,1,1,1,0,0,1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,1,0,1,0,0,0,1,
          0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,
          0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,
          1,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,
          0,1,0,1,2,1,1,1,1,1,1,1,2,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,
          1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,
          1,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
          0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0
          ), nrow=4, byrow=TRUE)
    
    rownames(expectedResult3) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult3) <- 
        c("AAee","A.Ae.e","A..Ae..e","A...Ae...e","AAei","A.Ae.i","A..Ae..i",
          "A...Ae...i","AAii","A.Ai.i","A..Ai..i","A...Ai...i","ACee","A.Ce.e",
          "A..Ce..i","A...Ce...i","ACii","A.Ci.i","AGee","A.Ge.e","A..Ge..e",
          "A...Ge...i","AGii","A.Gi.i","A..Gi..i","ATee","A.Te.e","A..Te..e",
          "A...Te...e","ATii","A.Ti.i","A..Ti..i","A...Ti...i","C.Ae.e",
          "C..Ae..e","C..Ae..i","C...Ae...i","C.Ai.i","C..Ai..i","CCee",
          "C.Ce.i","C...Ce...i","CCii","CGee","C.Ge.e","CGei","CGii","C.Gi.i",
          "C.Te.e","C..Te..e","C...Te...e","C.Ti.i","C..Ti..i","C...Ti...i",
          "GAee","G.Ae.e","G.Ae.i","G..Ae..i","G...Ae...i","GAii","G.Ai.i",
          "GCee","G..Ce..i","G...Ce...i","GCii","G.Ge.i","G...Ge...i","GTee",
          "G.Te.e","GTii","G.Ti.i","TAee","TAei","T.Ae.i","T..Ae..i",
          "T...Ae...i","TAii","T.Ce.i","T..Ce..i","T..Ge..i","T...Ge...i",
          "T...Te...i") 
    gappy1_3 <- gappyPairKernel(k=1, m=3, annSpec=TRUE, normalized=FALSE)
    
    testDenseAndSparseERAgainstKernelMatrix(x, gappy1_3,
                              expRes=expectedResult3, silent=silent)
    
    testStepName(silent, paste("Annotation specific Kernel with k=4,",
                               "m=1 - unnormalized"))
    
    expectedResult4 <- matrix(
        c(0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
          0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0), nrow=4, byrow=TRUE)
    
    rownames(expectedResult4) <- paste("S", 1:length(x), sep="")
    colnames(expectedResult4) <- 
        c("AAAG.AAAAeeee.eiii","AAAGTAAAeeeeeeii","AAGTAAAAeeeeeiii",
          "AAGT.AAAGeeee.iiii","ACCGATACeeeeeeii","ACCG.TACCeeee.eiii",
          "ACGTACGTeeeeiiii","AGTAAAAGeeeeiiii","AGTA.AAGTeeee.iiii",
          "ATACCGATeeiiiiii","CCGA.ACCGeeee.iiii","CCGATACCeeeeeiii",
          "CGATACCGeeeeiiii","CGAT.CCGAeeee.iiii","GATACCGAeeeiiiii",
          "GATA.CGATeeei.iiii","GTAAAAGTeeeiiiii","GTAA.AGTAeeei.iiii",
          "TAAAAGTAeeiiiiii") 
    gappy4_1 <- gappyPairKernel(k=4, m=1, annSpec=TRUE, normalized=FALSE)
    
    testDenseAndSparseERAgainstKernelMatrix(x, gappy4_1,
                                expRes=expectedResult4, silent=silent)
}

testKernelMatrix <- function(silent)
{
    testPosIndepKernelMatrix(silent)
    testPosSpecKernelMatrix(silent)
#    testAnnSpecKernelMatrix(silent)
    testDistWeightKernelMatrix(silent)
}

test_GappyPairKernel <- function(seed=789, silent=TRUE)
{
    set.seed(seed)
    testKernelMatrix(silent=silent)
    testExplicitRepresentation(silent=silent)
    testExplicitRepAnnotSpecific(silent=silent)
    testTooShortSequences(silent=silent)
}

