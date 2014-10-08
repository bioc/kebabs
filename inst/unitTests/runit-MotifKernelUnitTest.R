##
## Unit Test via RUnit for motif kernel
##

library(RUnit)

testGetExplicitRepAll <- function(xc, kernelDNA, kernelRNA, expectedResult,
                                  silent, sel, annCharset=NULL, ann1=NULL,
                                  ann2=NULL)
{
    dnaColNames <- kernelParameters(kernelDNA)$motifs
    rnaColNames <- kernelParameters(kernelRNA)$motifs
    aaColNames <- kernelParameters(kernelDNA)$motifs
    
    if (!is.null(annCharset))
    {
        annChars <- unlist(strsplit(annCharset, split=""))
        annotNames <- paste(rep(annChars, each=length(annChars)^2), 
                            rep(annChars, each=length(annChars)), 
                            annChars, sep="")
        dnaColNames <- paste(rep(dnaColNames, each=length(annotNames)),
                             annotNames, sep="")
        rnaColNames <- paste(rep(rnaColNames, each=length(annotNames)),
                             annotNames, sep="")
        aaColNames <- paste(rep(aaColNames, each=length(annotNames)),
                            annotNames, sep="")
    }

    testGetExplicitRepU(xc=xc, kernel=kernelDNA, silent=silent,
                        expectedResult=expectedResult, sel=sel,
                        colSubset=colSubset, annCharset=annCharset, 
                        ann1=ann1, ann2=ann2, dnaColNames=dnaColNames,
                        rnaColNames=rnaColNames, aaColNames=aaColNames,
                        kernel2=kernelRNA)
}

testPosIndepKernelMatrix <- function(silent)
{
    testName(silent, "Motif Kernel - Kernel Matrix Position Independent")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC")
    names(xc) <- namesX
    yc <- c("AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    names(yc) <- namesY
    zc <- "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGCATT"
    names(zc) <- namesZ
    
    motifs <- c("AG.", "GG[CT]", "T[^T]T")
    
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - y")
    
    symres1u <- defMat(matrix(c(26,17,13,17,17,15,13,15,14), nrow=3), 
                       namesX, namesX)
    symres2u <- defMat(matrix(c(38,33,33,34), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(31,22,17,27,23,20), nrow=3), namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized x - y",
                               "with subsetting"))
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), mot1,
                           expectedResUnnormalized, silent, 
                           subsetX=c(3:5), subsetY=c(3:4), kernel2=mot2)
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - y")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    
    
    testStepName(silent, "Pos Indep - kmer counts unnormalized x - z")

    symres1u <- defMat(matrix(c(26,17,13,17,17,15,13,15,14), nrow=3), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(19), nrow=1), namesZ, namesZ)
    
    asymresu <- defMat(matrix(c(18,17,14), nrow=3), namesX, namesZ)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)
    
    testGetKernelMatrixAll(xc, zc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)
    
    
    testStepName(silent, "Pos Indep - kmer counts normalized x - z")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=TRUE)
    
    testGetKernelMatrixAll(xc, zc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    

    testStepName(silent, "Pos Indep - kmer presence unnormalized")
    
    symres1u <- defMat(matrix(c(3,3,3,3,3,3,3,3,3), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(3,3,3,3), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(3,3,3,3,3,3), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, normalized=FALSE, presence=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE,
                        presence=TRUE)

    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)
    

    testStepName(silent, "Pos Indep - kmer presence normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, normalized=TRUE, presence=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=TRUE,
                        presence=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    
    
    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "translate lower to upper"))
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AGActtaaggGACCtggtcACCACGCTCGgtgAGgGGGACGGGgtgTAGG", 
            "CCCGgctaATTTTTTTGTattTTTAGtAAaGACGGGGTTTCACCATGTTG", 
            "cacaTGTGCCCTCTgggCCTGGTCACCCCACCACCTGCCCCcagtGTGAC")
    names(xc) <- namesX
    yc <- c("AAGCATGGTGGGATtggcaCAGGTGCACAcagAGAGATcaatGTATACAA", 
            "TTTAGAGaacTGGGTCTTGCTATGTGGCCCAGGTTGGccTCAAACTCCTG")
    names(yc) <- namesY
    zc <- "ATAAAGGTTGCAGacATCATGTCCTTtttgTCCCTAATtatTTCAGCATT"
    names(zc) <- namesZ
    
    symres1u <- defMat(matrix(c(26,17,13,17,17,15,13,15,14), nrow=3), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(38,33,33,34), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(31,22,17,27,23,20), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, normalized=FALSE, ignoreLower=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE,
                        ignoreLower=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)

    testStepName(silent, paste("Pos Indep - kmer counts unnormalized",
                               "ignore lowercase"))
    
    symres1u <- defMat(matrix(c(4,0,0,0,5,5,0,5,5), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(21,17,17,22), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(8,4,4,4,9,9), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, normalized=FALSE, ignoreLower=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE,
                        ignoreLower=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2, biovectorsOnly=TRUE )
}

testPosSpecKernelMatrix <- function(silent)
{
    testName(silent, "Motif Kernel - Kernel Matrix Position Specific")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    namesZ <- c("S6")
    
    xc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC")
    names(xc) <- namesX
    yc <- c("AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    names(yc) <- namesY
    zc <- "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGCATT"
    names(zc) <- namesZ
    
    motifs <- c("AG.", "GG[CT]", "T[^T]T")
    
    testStepName(silent, "Pos Spec - unnormalized")
    
    symres1u <- defMat(matrix(c(8,0,2,0,7,0,2,0,6), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(10,1,1,10), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(2,0,1,0,1,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, distWeight=1, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)

    
    testStepName(silent, "Pos Spec - unnormalized with subsetting")
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    
    mot1 <- motifKernel(motifs, distWeight=1, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=FALSE)

    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), mot1,
                        expectedResUnnormalized, silent, 
                        subsetX=c(3:5), subsetY=c(3:4), kernel2=mot2)
    
    testStepName(silent, "Pos Spec - normalized")

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    mot1 <- motifKernel(motifs, distWeight=1, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    
    testStepName(silent, "Pos Spec - attr pos on x")

    symres1u <- defMat(matrix(c(8,1,2,1,7,1,2,1,6), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(10,1,1,10), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,1,1,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, distWeight=1, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=FALSE)
    
    posX <- c(1,2,1)
    posY <- NULL
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized,
                           silent=silent, posX=posX, posY=posY, 
                           kernel2=mot2)
    
    testStepName(silent, "Pos Spec - attr pos on y")

    symres1u <- defMat(matrix(c(8,0,2,0,7,0,2,0,6), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(10,3,3,10), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,0,0,0,0,0), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, distWeight=1, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=FALSE)
    
    posX <- NULL
    posY <- c(2,12)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized,
                           silent=silent, posX=posX, posY=posY,
                           kernel2=mot2)
    
    testStepName(silent, "Pos Spec - attr pos on x and y")
    
    symres1u <- defMat(matrix(c(8,1,2,1,7,1,2,1,6), nrow=3), namesX, namesX)
    
    symres2u <- defMat(matrix(c(10,3,3,10), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(1,0,1,0,0,2), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, distWeight=1, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=1,
                        normalized=FALSE)
    
    posX <- c(1,2,1)
    posY <- c(2,12)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized,
                           silent=silent, posX=posX, posY=posY,
                           kernel2=mot2)
}

testAnnSpecKernelMatrix <- function(silent)
{
    testName(silent, "Motif Kernel - Kernel Matrix Annotation Specific")
    
    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    
    xc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC")
    names(xc) <- namesX
    yc <- c("AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    names(yc) <- namesY
    
    motifs <- c("AG.", "GG[CT]", "T[^T]T")
    
    
    testStepName(silent, "Ann Spec - unnormalized - annotation 1")
    
    symres1u <- defMat(matrix(c(14,6,5,6,9,8,5,8,10), nrow=3), namesX, namesX)
    symres2u <- defMat(matrix(c(26,10,10,24), nrow=2), namesY, namesY)
    asymresu <- defMat(matrix(c(15,8,11,13,11,9), nrow=3), namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, annSpec=TRUE, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), annSpec=TRUE,
                        normalized=FALSE)
    
    annX1 <- c(rep("aaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbb",3))
    annY1 <- c(rep("aaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbb",2))
    annCharset1 <- "ab"
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, annCharset=annCharset1, annX=annX1,
                           annY=annY1, kernel2=mot2)
    
    
    testStepName(silent, "Ann Spec - normalized - annotation 1")
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, annSpec=TRUE, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), annSpec=TRUE,
                        normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, silent,
                           annCharset=annCharset1, annX=annX1, 
                           annY=annY1, kernel2=mot2)
    
    testStepName(silent, paste("Ann Spec - normalized - annotation 1 with",
                               "subsetting"))
    
    xc1 <- c("ACCAAC",
             "GACATA")
    names(xc1) <- c("NONAME1", "NONAME2")
    annX2 <- c("aaabbb", "aaabbb")
    
    mot1 <- motifKernel(motifs, annSpec=TRUE, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), annSpec=TRUE,
                        normalized=FALSE)
    
    testGetKernelMatrixAll(c(xc1,xc,xc1), c(xc1,yc,xc1), mot1,
                           expectedResUnnormalized, silent, 
                           annCharset=annCharset1, 
                           annX=c(annX2,annX1,annX2), 
                           annY=c(annX2,annY1,annX2), 
                           subsetX=c(3:5), subsetY=c(3:4), kernel2=mot2)
}

testDistWeightKernelMatrix <- function(silent)
{
    testName(silent, "Motif Kernel - Kernel Matrix with Distance Weights")

    namesX <- c("S1", "S2", "S3")
    namesY <- c("S4", "S5")
    
    xc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC")
    names(xc) <- namesX
    yc <- c("AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    names(yc) <- namesY
    
    motifs <- c("AG.", "GG[CT]", "T[^T]T")
    

    testStepName(silent, "Dist Weight - Explicit weight vector")
    
    mot1_1 <- motifKernel(motifs, distWeight=rep(1,50), normalized=FALSE)
    mot1_2 <- motifKernel(chartr("Tt", "Uu", motifs), distWeight=rep(1,50),
                          normalized=FALSE)
    mot2_1 <- motifKernel(motifs, normalized=FALSE)
    mot2_2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)

    compareEqualKernelsAll(xc, yc, mot1_1, mot2_1, silent, kernel3=mot1_2,
                           kernel4=mot2_2)
    
    
    testStepName(silent, "Dist Weight - Predefined weight functions")
    
    symres1u <- defMat(matrix(c(8.0, 1.2, 2.0,
                                1.2, 8.2, 1.0,
                                2.0, 1.0, 6.0), nrow=3, byrow=TRUE), 
                       namesX, namesX)
    
    symres2u <- defMat(matrix(c(14.0,3.2,3.2,12.8), nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(5.2,2.6,3.0,3.0,3.6,0.8), nrow=3), 
                       namesX, namesY)

    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))
    
    mot1 <- motifKernel(motifs, distWeight=linWeight(sigma=5), 
                        normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=linWeight(sigma=5), normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))

    mot1 <- motifKernel(motifs, distWeight=linWeight(sigma=5), 
                        normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=linWeight(sigma=5), normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, silent,
                           kernel2=mot2)
    
    symres1u <- defMat(matrix(c(8.889604, 2.280642, 3.166891,
                                2.280642, 9.092075, 2.057054,
                                3.166891, 2.057054, 7.232344), 
                              nrow=3, byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(16.584441, 5.441963, 5.441963, 15.695080),
                              nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(7.021956, 5.304929,
                                4.893842, 5.866218,
                                4.527119, 2.616651), nrow=3, byrow=TRUE),
                       namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    mot1 <- motifKernel(motifs, distWeight=expWeight(sigma=5), 
                        normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=expWeight(sigma=5), normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)

    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, distWeight=expWeight(sigma=5), 
                        normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=expWeight(sigma=5), normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    
    symres1u <- defMat(matrix(c(8.285069, 2.020650, 2.832538,
                                2.020650, 9.440046, 1.738522,
                                2.832538, 1.738522, 7.017476), 
                              nrow=3, byrow=TRUE), namesX, namesX)
    
    symres2u <- defMat(matrix(c(16.763600, 5.315443, 5.315443, 15.604431),
                              nrow=2), namesY, namesY)
    
    asymresu <- defMat(matrix(c(6.984183, 5.698195,
                                5.006009, 6.015249,
                                4.362391, 1.935495), nrow=3, byrow=TRUE),
                       namesX, namesY)
    
    expectedResUnnormalized <- list(as.KernelMatrix(symres1u), 
                                    as.KernelMatrix(symres2u), 
                                    as.KernelMatrix(asymresu))

    mot1 <- motifKernel(motifs, distWeight=gaussWeight(sigma=5),
                        normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=gaussWeight(sigma=5), normalized=FALSE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResUnnormalized, 
                           silent, kernel2=mot2)
    
    expectedResNormalized <- list(normalizeMatrix(symres1u), 
                                  normalizeMatrix(symres2u), 
                                  normalizeMatrix(asymresu, diag(symres1u), 
                                                            diag(symres2u)))
    
    mot1 <- motifKernel(motifs, distWeight=gaussWeight(sigma=5), 
                        normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), 
                        distWeight=gaussWeight(sigma=5), normalized=TRUE)
    
    testGetKernelMatrixAll(xc, yc, mot1, expectedResNormalized, 
                           silent, kernel2=mot2)
    
    
    testStepName(silent, "Dist Weight - User defined weight functions")
    
    mot1_1 <- motifKernel(motifs, 
        distWeight=function(d, sigma=5){exp(-abs(d)/sigma)}, normalized=TRUE)
    mot1_2 <- motifKernel(chartr("Tt", "Uu", motifs), 
        distWeight=function(d, sigma=5){exp(-abs(d)/sigma)}, normalized=TRUE)
    mot2_1 <- motifKernel(motifs, 
        distWeight=expWeight(sigma=5), normalized=TRUE)
    mot2_2 <- motifKernel(chartr("Tt", "Uu", motifs), 
        distWeight=expWeight(sigma=5), normalized=TRUE)
    
    compareEqualKernelsAll(xc, yc, mot1_1, mot2_1, silent, kernel3=mot1_2,
                           kernel4=mot2_2)
}

testExplicitRepresentation <- function(silent)
{
    namesX <- c("S1", "S2", "S3", "S4", "S5")
    
    dc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC",
            "AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    
    names(dc) <- namesX
    
    motifs <- c("AG.", "GG[CT]", "T[^T]T")

    
    testName(silent, "Motif Kernel - Explicit Representation")
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - unnormalized"))

    erd1u <- defMat(matrix(c(4, 3, 1,
                             2, 2, 3,
                             1, 2, 3,
                             5, 3, 2,
                             3, 4, 3), nrow=5, ncol=3, byrow=TRUE), 
                    namesX, motifs)
    
    erd2u <- defMat(matrix(0, nrow=5, ncol=3, byrow=TRUE), namesX, motifs)
    
    expectedResUnnormalized <- list(erd1u, erd1u) 
    
    
    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)

    testGetExplicitRepAll(dc, mot1, mot2, expectedResUnnormalized, silent)

    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), mot1, mot2, 
                          expectedResUnnormalized, silent, sel=c(3:7))
    

    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    mot1 <- motifKernel(motifs, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=TRUE)
    
    testGetExplicitRepAll(dc, mot1, mot2, expectedResUnnormalized, silent)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))

    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)

    testERAgainstKernelMatrixAll(dc, mot1, silent, kernel2=mot2)

    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))

    mot1 <- motifKernel(motifs, normalized=TRUE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, mot1, silent, kernel2=mot2)
    
    testStepName(silent, paste("Explict Rep Sparse Subsetting against",
                               "Explicit Rep Dense"))
    
    motifs <- c("AG.", "C[AG]C", "A.[GT]", "[^AC]GG", "GG[CT]", "T[^T]T")

    mot1 <- motifKernel(motifs, normalized=FALSE)
    mot2 <- motifKernel(chartr("Tt", "Uu", motifs), normalized=FALSE)
    
    testSubsettingOfSparseERAgainstDenseERAll(dc, mot1, silent, kernel2=mot2)
}

testExplicitRepAnnotSpecific <- function(silent)
{
    namesX <- c("S1", "S2", "S3", "S4", "S5")
    
    dc <- c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGTAGG", 
            "CCCGGCTAATTTTTTTGTATTTTTAGTAAAGACGGGGTTTCACCATGTTG", 
            "CACATGTGCCCTCTGGGCCTGGTCACCCCACCACCTGCCCCCAGTGTGAC",
            "AAGCATGGTGGGATTGGCACAGGTGCACACAGAGAGATCAATGTATACAA", 
            "TTTAGAGAACTGGGTCTTGCTATGTGGCCCAGGTTGGCCTCAAACTCCTG")
    
    names(dc) <- namesX
    
    motifsDNA <- c("AG.", "A.[GT]", "C[AG]C", "GG[CT]", "T[^T]T","[^AC]GG")
    motifsRNA <- c("AG.", "A.[GU]", "C[AG]C", "GG[CU]", "U[^U]U", "[^AC]GG")
        
    testName(silent, paste("Motif Kernel - Explicit Representation",
                           "Annotation Specific"))
    
    testStepName(silent, paste("Explict Rep Ann Spec with/without zero",
                               "features - unnormalized"))
    
    erd1u <- defMat(matrix(c(2,0,2,4,0,3,2,1,0,1,2,0,1,2,0,5,
                             0,1,1,3,1,3,0,0,1,1,1,2,1,0,0,2,
                             0,0,1,1,0,1,1,1,2,2,0,2,1,3,0,0,
                             2,0,3,4,0,2,1,0,2,3,0,0,2,4,0,0,
                             2,0,1,2,0,2,0,0,0,1,3,3,0,2,1,1
                             ), nrow=5, byrow=TRUE))
    
    rownames(erd1u) <- namesX
    colnames(erd1u) <- c("AG.aaa","AG.abb","AG.bbb","A.[GT]aaa",
                         "A.[GT]abb","A.[GT]bbb","C[AG]Caaa","C[AG]Caab",
                         "C[AG]Cbbb","GG[CT]aaa","GG[CT]bbb","T[^T]Taaa",
                         "T[^T]Tbbb","[^AC]GGaaa","[^AC]GGabb","[^AC]GGbbb")
    
    expectedResUnnormalized <- list(erd1u, erd1u) 
    
    annCharset <- "ab"
    annX <- c(rep("aaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbb", 
                  length(dc)))
    
    motDNA <- motifKernel(motifsDNA, annSpec=TRUE, normalized=FALSE)
    motRNA <- motifKernel(motifsRNA, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(dc, motDNA, motRNA, expectedResUnnormalized, 
                          silent, annCharset=annCharset, ann1=annX)
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features with subsetting - unnormalized"))
    
    dc1 <- c("ACCAAC",
             "GACATA")
    
    names(dc1) <- c("NONAME1", "NONAME2")
    
    
    motDNA <- motifKernel(motifsDNA, annSpec=TRUE, normalized=FALSE)
    motRNA <- motifKernel(motifsRNA, annSpec=TRUE, normalized=FALSE)
    
    testGetExplicitRepAll(c(dc1, dc, dc1), motDNA, motRNA, 
                          expectedResUnnormalized, silent, sel=c(3:7), 
                          annCharset=annCharset, 
                          ann1=c(c("aaabbb","bbbaaa"), 
                                 annX, c("aaabbb","bbbaaa")))
    
    
    testStepName(silent, paste("Explict Representation with/without zero",
                               "features - normalized"))
    
    motDNA <- motifKernel(motifsDNA, annSpec=TRUE, normalized=TRUE)
    motRNA <- motifKernel(motifsRNA, annSpec=TRUE, normalized=TRUE)
    
    testGetExplicitRepAll(dc, motDNA, motRNA, expectedResUnnormalized, 
                          silent, annCharset=annCharset, ann1=annX)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - unnormalized"))
    
    motDNA <- motifKernel(motifsDNA, annSpec=TRUE, normalized=FALSE)
    motRNA <- motifKernel(motifsRNA, annSpec=TRUE, normalized=FALSE)
    
    testERAgainstKernelMatrixAll(dc, motDNA, silent, annCharset=annCharset,
                                 ann1=annX, kernel2=motRNA)
    
    
    testStepName(silent, paste("Explict Rep with/without zero features",
                               "against kernel matrix - normalized"))
    
    motDNA <- motifKernel(motifsDNA, annSpec=TRUE, normalized=TRUE)
    motRNA <- motifKernel(motifsRNA, annSpec=TRUE, normalized=TRUE)
    
    testERAgainstKernelMatrixAll(dc, motDNA, silent, annCharset=annCharset,
                                 ann1=annX, kernel2=motRNA)
}

testKernelMatrix <- function(silent)
{
    testPosIndepKernelMatrix(silent)
    testPosSpecKernelMatrix(silent)
    testAnnSpecKernelMatrix(silent)
    testDistWeightKernelMatrix(silent)
}

test_MotifKernel <- function(seed=789, silent=TRUE)
{
    set.seed(seed)
    testKernelMatrix(silent)
    testExplicitRepresentation(silent)
    testExplicitRepAnnotSpecific(silent)
}

