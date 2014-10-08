##
## common parts of kernel test
##

getNonzeroCols <- function(er)
{
    mat <- as(er, "matrix")
    which(colSums(mat) > 0)
}

testGetKernelMatrix <- function(x, y, kernel, expectedResult, silent, 
                                subsetX=NULL, subsetY=NULL, posX=NULL,
                                posY=NULL, annCharset=NULL, annX=NULL,
                                annY=NULL)
{
    if (!silent) cat("\nClass: ", class(x), "\n")
    
    if (!is.null(posX))
        positionMetadata(x) <- posX

    if (!is.null(annX))
        annotationMetadata(x, annCharset=annCharset) <- annX

    if (!is.null(posY))
        positionMetadata(y) <- posY

    if (!is.null(annY))
        annotationMetadata(y, annCharset=annCharset) <- annY

    km <- NULL
    km1 <- NULL
    km2 <- NULL
    km3 <- NULL

    if (is.null(subsetX) && is.null(subsetY))
    {
        km <- kernel(x)
        km1 <- kernel(y)
        km2 <- kernel(x,y)
        km3 <- kernel(y,x)
    }
    else
    {
        km <- kernel(x, selx=subsetX)
        km1 <- kernel(y, selx=subsetY)
        km2 <- kernel(x,y, selx=subsetX, sely=subsetY)
        km3 <- kernel(y,x, selx=subsetY, sely=subsetX)
    }
    if (!silent)
    {
        print(km)
        print(km1)
        print(km2)
        print(km3)
    }
    checkEquals(km, expectedResult[[1]], tolerance=5e-7)
    checkEquals(km, t(km))
    checkEquals(km1, expectedResult[[2]], tolerance=5e-7)
    checkEquals(km1, t(km1))
    checkEquals(km2, expectedResult[[3]], tolerance=5e-7)
    checkEquals(km3, t(km2))
}

testGetKernelMatrixAll <- function(xc, yc, kernel, expectedResult, silent, 
                                   subsetX=NULL, subsetY=NULL, posX=NULL,
                                   posY=NULL, annCharset=NULL, annX=NULL,
                                   annY=NULL, biovectorsOnly=FALSE, 
                                   kernel2=NULL)
{
    if (biovectorsOnly)
        range <- c(2,4,6)
    else
        range <- 1:6
    
    for (i in range)
    {
        currKernel <- kernel
        switch(i,
               "1" = {x <- DNAStringSet(chartr("Uu", "Tt", xc)); 
                      names(x) <- names(xc);
                      y <- DNAStringSet(chartr("Uu", "Tt", yc));
                      names(y) <- names(yc)},
               "2" = {x <- DNAVector(chartr("Uu", "Tt", xc)); 
                      names(x) <- names(xc);
                      y <- DNAVector(chartr("Uu", "Tt", yc));
                      names(y) <- names(yc)},
               "3" = {x <- RNAStringSet(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc); 
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      y <- RNAStringSet(chartr("Tt", "Uu", yc));
                      names(y) <- names(yc)},
               "4" = {x <- RNAVector(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc);
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      y <- RNAVector(chartr("Tt", "Uu", yc));
                      names(y) <- names(yc)},
               "5" = {x <- AAStringSet(xc);
                      y <- AAStringSet(yc)},
               "6" = {x <- AAVector(xc);
                      y <- AAVector(yc)}
               )

        if (is(currKernel, "MismatchKernel") && 
            (class(x) %in% c("AAStringSet","AAVector")))
            expResult <- expectedResult[4:6]
        else
            expResult <- expectedResult[1:3]


        testGetKernelMatrix(x=x, y=y, kernel=currKernel, silent=silent,
                            expectedResult=expResult, subsetX=subsetX,
                            subsetY=subsetY, posX=posX, posY=posY,
                            annCharset=annCharset, annX=annX, annY=annY)
    }
}

compareEqualKernels <- function(x, y, kernel1, kernel2, silent)
{        
    if (!silent) cat("\nClass: ", class(x), "\n")
        
    kml <- list(NA, 8)
    kml[[1]] <- kernel1(x)
    checkEquals(dim(kml[[1]]), c(length(x), length(x)))
    kml[[2]] <- kernel2(x)
    kml[[3]] <- kernel1(y)
    checkEquals(dim(kml[[3]]), c(length(y), length(y)))
    kml[[4]] <- kernel2(y)
    kml[[5]] <- kernel1(x,y)
    checkEquals(dim(kml[[5]]), c(length(x), length(y)))
    kml[[6]] <- kernel2(x,y)
    kml[[7]] <- kernel1(y,x)
    checkEquals(dim(kml[[7]]), c(length(y), length(x)))
    kml[[8]] <- kernel2(y,x)

    for (i in 1:4)
    {
        if (!silent)
        {
            print(kml[[2*i-1]])
            print(kml[[2*i]])
        }
        checkTrue(!any(is.na(kml[[2*i-1]])))
        checkEquals(kml[[2*i-1]], kml[[2*i]])
    }
}

compareEqualKernelsAll <- function(xc, yc, kernel1, kernel2, silent,
                                   kernel3=NULL, kernel4=NULL)
{
    for (i in 1:6)
    {
        currKernel1 <- kernel1
        currKernel2 <- kernel2
 
        switch(i,
               "1" = {x <- DNAStringSet(chartr("Uu", "Tt", xc)); 
                      names(x) <- names(xc);
                      y <- DNAStringSet(chartr("Uu", "Tt", yc));
                      names(y) <- names(yc)},
               "2" = {x <- DNAVector(chartr("Uu", "Tt", xc)); 
                      names(x) <- names(xc);
                      y <- DNAVector(chartr("Uu", "Tt", yc));
                      names(y) <- names(yc)},
               "3" = {x <- RNAStringSet(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc);
                      if (!is.null(kernel3))
                          { currKernel1 <- kernel3; currKernel2 <- kernel4 }; 
                      y <- RNAStringSet(chartr("Tt", "Uu", yc));
                      names(y) <- names(yc)},
               "4" = {x <- RNAVector(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc);
                      if (!is.null(kernel3))
                          { currKernel1 <- kernel3; currKernel2 <- kernel4 }; 
                      y <- RNAVector(chartr("Tt", "Uu", yc));
                      names(y) <- names(yc)},
               "5" = {x <- AAStringSet(xc);
                      y <- AAStringSet(yc)},
               "6" = {x <- AAVector(xc);
                      y <- AAVector(yc)}
               )

        compareEqualKernels(x=x, y=y, kernel1=currKernel1, 
                            kernel2=currKernel2, silent=silent)
    }
}

testGetExplicitRep <- function(x, kernel, expectedResult, silent, sel,
                               features=NULL, colSubset=NULL, colNames, 
                               annCharset=NULL, ann1=NULL, ann2=NULL)
{
    if (!silent) cat("\nClass: ", class(x), "\n")
        
    if (!is.null(ann1))
    {
        if (class(x) %in% c("DNAStringSet", "DNAVector", 
                           "RNAStringSet", "RNAVector") || is.null(ann2))
            annotationMetadata(x, annCharset=annCharset) <- ann1
        else
            annotationMetadata(x, annCharset=annCharset) <- ann2
    }

    erd1 <- NULL
    erd2 <- NULL
    erd3 <- NULL
    ers1 <- NULL
    ers2 <- NULL
    ers3 <- NULL
    erd1f <- NULL
    erd2f <- NULL
    erd3f <- NULL
    ers1f <- NULL
    ers2f <- NULL
    ers3f <- NULL

    if (missing(sel))
    {
        erd1 <- getExRep(x, kernel=kernel, sparse=FALSE, features=features)
        ers1 <- getExRep(x, kernel=kernel, sparse=TRUE, features=features)

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            erd2 <- getExRep(x, kernel=kernel, sparse=FALSE, 
                             zeroFeatures=TRUE, features=features)
            erd3 <- getExRep(x, kernel=kernel, sparse=FALSE, 
                             useRowNames=FALSE, useColNames=FALSE, 
                             zeroFeatures=TRUE, features=features)
            ers2 <- getExRep(x, kernel=kernel, sparse=TRUE, 
                             zeroFeatures=TRUE, features=features)
            ers3 <- getExRep(x, kernel=kernel, sparse=TRUE, 
                             useRowNames=FALSE, useColNames=FALSE, 
                             zeroFeatures=TRUE, features=features)
        }

        if (!is.null(features))
        {
            erdf <- getExRep(x, kernel=kernel, sparse=FALSE, 
                             features=NULL, zeroFeatures=TRUE)
            erd1f <- erdf[, features]
            erd1f <- erd1f[,getNonzeroCols(erd1f)]
            erd2f <- erdf[, features]
            erd3f <- erd2f
            dimnames(erd3f) <- list(NULL, NULL)
            ersf <- getExRep(x, kernel=kernel, sparse=TRUE, 
                             features=NULL, zeroFeatures=TRUE)
            ers1f <- ersf[, features]
            ers1f <- ers1f[,getNonzeroCols(ers1f)]
            ers2f <- ersf[, features]
            ers3f <- ers2f
            dimnames(ers3f) <- list(NULL, NULL)
        }
    }
    else
    {
        erd1 <- getExRep(x, kernel=kernel, sparse=FALSE, 
                         selx=sel, features=features)
        ers1 <- getExRep(x, kernel=kernel, sparse=TRUE, selx=sel,
                         features=features)

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            erd2 <- getExRep(x, kernel=kernel, sparse=FALSE, 
                             zeroFeatures=TRUE, selx=sel, features=features)
            erd3 <- getExRep(x, kernel=kernel, sparse=FALSE, 
                             useRowNames=FALSE, useColNames=FALSE, 
                             zeroFeatures=TRUE, selx=sel, features=features)
            ers2 <- getExRep(x, kernel=kernel, sparse=TRUE, zeroFeatures=TRUE,
                             selx=sel, features=features)
            ers3 <- getExRep(x, kernel=kernel, sparse=TRUE, useRowNames=FALSE,
                             useColNames=FALSE, zeroFeatures=TRUE, selx=sel,
                             features=features)
        }

        if (!is.null(features))
        {
            erdf <- getExRep(x, kernel=kernel, sparse=FALSE, selx=sel,
                             features=NULL, zeroFeatures=TRUE)
            erd1f <- erdf[, features]
            erd1f <- erd1f[,getNonzeroCols(erd1f)]
            erd2f <- erdf[, features]
            erd3f <- erd2f
            dimnames(erd3f) <- list(NULL, NULL)
            ersf <- getExRep(x, kernel=kernel, sparse=TRUE, selx=sel,
                             features=NULL, zeroFeatures=TRUE)
            ers1f <- ersf[, features]
            ers1f <- ers1f[,getNonzeroCols(ers1f)]
            ers2f <- ersf[, features]
            ers3f <- ers2f
            dimnames(ers3f) <- list(NULL, NULL)
        }            
    }

    if (!silent)
    {
        print(erd1)
        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            print(erd2)
            print(erd3)
        }
        print(ers1)
        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            print(ers2)
            print(ers3)
        }
    }
    
    normalizeERD <- kernelParameters(kernel)$normalized

    if (class(x) %in% c("DNAStringSet", "DNAVector", 
                       "RNAStringSet", "RNAVector"))
    {
        expRes1 <- expectedResult[[1]]

        if (class(x) %in% c("RNAStringSet", "RNAVector"))
            colnames(expRes1) <- chartr("Tt", "Uu", colnames(expRes1))

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            ## expand result with zero features
            expRes2 <- matrix(0,nrow(expectedResult[[1]]), length(colNames))
            colnames(expRes2) <- colNames
            rownames(expRes2) <- rownames(expectedResult[[1]])

            if (class(x) %in% c("RNAStringSet", "RNAVector"))
            {
                expRes2[ , chartr("Tt", "Uu", colnames(expectedResult[[1]]))] <-
                    expectedResult[[1]]
            }
            else
                expRes2[ , colnames(expectedResult[[1]])] <- expectedResult[[1]]

            expRes3 <- expRes2
            ## delete names for expRes3 below
        }
    }
    else
    {
        expRes1 <- expectedResult[[2]]

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            expRes2 <- matrix(0,nrow(expectedResult[[2]]), length(colNames))
            rownames(expRes2) <- rownames(expectedResult[[2]])
            colnames(expRes2) <- colNames
            expRes2[ , colnames(expectedResult[[2]])] <- expectedResult[[2]]
            expRes3 <- expRes2
            ## delete names for expRes3 below
        }
    }
        
    if (normalizeERD)
    {
        expRes1 <- normalizeERD(expRes1)

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            expRes2 <- normalizeERD(expRes2)
            expRes3 <- normalizeERD(expRes3)
        }
    }

    if (!is.null(colSubset))
    {
        expRes1 <- expRes1[, colSubset[[1]]]

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            expRes2 <- expRes2[, colSubset[[2]]]
            expRes3 <- expRes3[, colSubset[[3]]]
        }
    }

    expRes1 <- new("ExplicitRepresentationDense", expRes1)
    expRes1@usedKernel <- kernel

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        expRes2 <- new("ExplicitRepresentationDense", expRes2)
        expRes2@usedKernel <- kernel
        expRes3 <- new("ExplicitRepresentationDense", expRes3)
        expRes3@usedKernel <- kernel
    }

    expResS1 <- as(expRes1, "ExplicitRepresentationSparse")
    expResS1@usedKernel <- kernel

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        expResS2 <- as(expRes2, "ExplicitRepresentationSparse")
        expResS2@usedKernel <- kernel
        expResS3 <- as(expRes3, "ExplicitRepresentationSparse")
        expResS3@usedKernel <- kernel

        ## remove names from third variant
        dimnames(expRes3) <- list(NULL, NULL)
        dimnames(expResS3) <- list(NULL, NULL)
    }

    if (!silent)
    {
        print(str(expRes1))

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            print(str(expRes2))
            print(str(expRes3))
        }

        print(str(expResS1))

        if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
        {
            print(str(expResS2))
            print(str(expResS3))
        }
    }

    checkEquals(erd1, expRes1, tolerance=1e-7)

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        checkEquals(erd2, expRes2, tolerance=1e-7)
        checkEquals(erd3, expRes3, tolerance=1e-7)
    }

    checkEquals(ers1, expResS1, tolerance=1e-7)

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        checkEquals(ers2, expResS2, tolerance=1e-7)
        checkEquals(ers3, expResS3, tolerance=1e-7)
    }

    if (!is.null(features))
    {
        checkEquals(erd1, erd1f, tolerance=1e-7)
        checkEquals(erd2, erd2f, tolerance=1e-7)
        checkEquals(erd3, erd3f, tolerance=1e-7)
        checkEquals(ers1, ers1f, tolerance=1e-7)
        checkEquals(ers2, ers2f, tolerance=1e-7)
        checkEquals(ers3, ers3f, tolerance=1e-7)
    }
}

testGetExplicitRepU <- function(xc, kernel, expectedResult, silent, sel,
                                features=NULL, colSubset=NULL, 
                                annCharset=NULL, ann1=NULL, ann2=NULL, 
                                dnaColNames, rnaColNames, aaColNames,
                                kernel2=NULL)
{    
    dss <- DNAStringSet(chartr("Uu", "Tt", xc));
    names(dss) <- names(xc)

    if (is(kernel, "MotifKernel"))
        ass <- dss
    else
        ass <- translate(dss)

    feat <- NULL
    colsub <- NULL
    
    for (i in 1:6)
    {
        currKernel <- kernel
        switch(i,
               "1" = {x <- dss; if (!is.null(features)) feat <- features[[1]]; 
                      if (!is.null(features)) colsub <- colSubset[1:3]; 
                      colNames <- dnaColNames},
               "2" = {x <- DNAVector(chartr("Uu", "Tt", xc))
                      if (!is.null(features)) feat <- features[[1]];
                      if (!is.null(features)) colsub <- colSubset[1:3]; 
                      colNames <- dnaColNames},
               "3" = {x <- RNAStringSet(chartr("Tt", "Uu", xc));
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      if (!is.null(features)) 
                          feat <- chartr("Tt", "Uu", features[[1]]);
                      if (!is.null(features)) colsub <- colSubset[1:3]; 
                      colNames <- rnaColNames},
               "4" = {x <- RNAVector(chartr("Tt", "Uu", xc));
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      if (!is.null(features)) 
                          feat <- chartr("Tt", "Uu", features[[1]]); 
                      if (!is.null(features)) colsub <- colSubset[1:3]; 
                      colNames <- rnaColNames},
               "5" = {x <- ass; names(x) <- names(xc);
                      if (!is.null(features)) feat <- features[[2]];
                      if (!is.null(features)) colsub <- colSubset[4:6]; 
                      colNames <- aaColNames},
               "6" = {x <- AAVector(as.character(ass)); names(x) <- names(xc);
                      if (!is.null(features)) feat <- features[[2]];
                      if (!is.null(features)) colsub <- colSubset[4:6]; 
                      colNames <- aaColNames}
               )
        
        testGetExplicitRep(x=x, kernel=currKernel, silent=silent, 
                           expectedResult=expectedResult, sel=sel,
                           features=feat, colSubset=colsub, 
                           colNames=colNames, annCharset=annCharset,
                           ann1=ann1, ann2=ann2)
    }    
}

testERAgainstKernelMatrix <- function(x, kernel, silent, annCharset,
                                      ann1, ann2)
{        
    if (!silent) cat("\nClass: ", class(x), "\n")
    
    if (!is.null(ann1))
    {
        if (class(x) %in% c("DNAStringSet", "DNAVector", 
                            "RNAStringSet", "RNAVector") || is.null(ann2))
            annotationMetadata(x, annCharset=annCharset) <- ann1
        else
            annotationMetadata(x, annCharset=annCharset) <- ann2
    }
    
    erd1 <- NULL
    erd2 <- NULL
    erd3 <- NULL
    erd1 <- getExRep(x, kernel=kernel, sparse=FALSE)
    kmerd1 = linearKernel(x=erd1)

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        erd2 <- getExRep(x, kernel=kernel, sparse=FALSE, zeroFeatures=TRUE)
        erd3 <- getExRep(x, kernel=kernel, sparse=FALSE, useRowNames=FALSE,
                         useColNames=FALSE, zeroFeatures=TRUE)
        kmerd2 = linearKernel(x=erd2)
        kmerd3 = linearKernel(x=erd3)
    }
    
    km1 <- NULL
    km2 <- NULL
    km1 <- kernel(x)
    km2 <- km1
    dimnames(km2) <- NULL
    
    if(!silent)
    {
        print(attributes(kmerd1))
        print(attributes(km1))
    }
    
    checkEquals(kmerd1, km1, tolerance=1e-7)

    if (!(is(kernel, "MotifKernel") && !is.null(ann1)))
    {
        checkEquals(kmerd2, km1, tolerance=1e-7)
        checkEquals(kmerd3, km2, tolerance=1e-7)
    }
}

testERAgainstKernelMatrixAll <- function(xc, kernel, silent, annCharset=NULL,
                                         ann1=NULL, ann2=NULL, kernel2=NULL)
{
    dss <- DNAStringSet(chartr("Uu", "Tt", xc));
    names(dss) <- names(xc)

    if (is(kernel, "MotifKernel"))
        ass <- dss
    else
        ass <- translate(dss)
    
    for (i in 1:6)
    {
        currKernel <- kernel
        switch(i,
               "1" = {x <- dss},
               "2" = {x <- DNAVector(chartr("Uu", "Tt", xc));
                      names(x) <- names(xc)},
               "3" = {x <- RNAStringSet(chartr("Tt", "Uu", xc));
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      names(x) <- names(xc)},
               "4" = {x <- RNAVector(chartr("Tt", "Uu", xc));
                      if (!is.null(kernel2)) currKernel <- kernel2;
                      names(x) <- names(xc)},
               "5" = {x <- ass; names(x) <- names(xc)},
               "6" = {x <- AAVector(as.character(ass));
                      names(x) <- names(xc)}
               )
        
        testERAgainstKernelMatrix(x=x, kernel=currKernel, silent=silent,
                                  annCharset=annCharset, ann1=ann1, 
                                  ann2=ann2)
    }
}

testDenseAndSparseERAgainstKernelMatrix <- function(seqs, kernel, expRes,
                                                    silent)
{
    FEATURE_SPACE_LIMIT <- 1000000
    erd <- getExRep(seqs, kernel, sparse=FALSE)
    ers <- getExRep(seqs, kernel, sparse=TRUE)
    
    if (!silent)
    {
        print(paste("FeatureSpaceDimension:", 
                    getFeatureSpaceDimension(kernel, seqs)))
        cat("\n")
        print(erd)
        print(ers)
    }
    
    checkEquals(erd@.Data, expRes)
    checkEquals(erd@.Data, as.matrix(ers))
    
    if (getFeatureSpaceDimension(kernel, seqs) <= FEATURE_SPACE_LIMIT)
    {
        erdz <- getExRep(seqs, kernel, sparse=FALSE, zeroFeatures=TRUE)
        ersz <- getExRep(seqs, kernel, sparse=TRUE, zeroFeatures=TRUE)
        
        if (!silent)
        {
            print(erdz)
            print(ersz)
        }
        
        checkEquals(erd, erdz[,which(colSums(erdz) > 0)])
        checkEquals(erdz@.Data, as.matrix(ersz))
    }
    
    km1 <- linearKernel(erd)
    km2 <- linearKernel(ers)
    km3 <- kernel(seqs)
    
    if (!silent)
    {
        print(km1)
        print(km2)
        print(km3)
    }    
    
    checkEquals(km1, km2)
    checkEquals(km1, km3)
}

testSubsettingOfSparseERAgainstDenseER <- function(x, kernel, silent,
                              annCharset=annCharset, ann1=ann1, ann2=ann2)
{
    if (!silent) cat("\nClass: ", class(x), "\n")
    
    if (!is.null(ann1))
    {
        if (class(x) %in% c("DNAStringSet", "DNAVector", 
                            "RNAStringSet", "RNAVector"))
            annotationMetadata(x, annCharset=annCharset) <- ann1
        else
            annotationMetadata(x, annCharset=annCharset) <- ann2
    }
    
    erd1 <- NULL
    erd1 <- getExRep(x, kernel=kernel, sparse=FALSE, zeroFeatures=TRUE)

    ers1 <- NULL
    ers1 <- getExRep(x, kernel=kernel, sparse=TRUE, zeroFeatures=TRUE)
    
    ## row subsetting
    erd1_1 <- NULL
    erd1_2 <- NULL
    erd1_3 <- NULL
    erd1_4 <- NULL
    erd1_5 <- NULL
    erd1_1 <- erd1[1,]
    erd1_2 <- erd1[2:3,]
    erd1_3 <- erd1[c(1,3,5),]
    erd1_4 <- erd1[-c(1,3,5),]
    erd1_5 <- erd1[c(names(x)[2], names(x[4])),]

    ers1_1 <- NULL
    ers1_2 <- NULL
    ers1_3 <- NULL
    ers1_4 <- NULL
    ers1_5 <- NULL
    ers1_1 <- ers1[1,]
    ers1_2 <- ers1[2:3,]
    ers1_3 <- ers1[c(1,3,5),]
    ers1_4 <- ers1[-c(1,3,5),]
    ers1_5 <- ers1[c(names(x)[2], names(x[4])),]
        
    if(!silent)
    {
        print(erd1_1)
        print(ers1_1)
        print(erd1_2)
        print(ers1_2)
        print(erd1_3)
        print(ers1_3)
        print(erd1_4)
        print(ers1_4)
        print(erd1_5)
        print(ers1_5)
    }

    checkEquals(as(ers1_1, "matrix"), as(erd1_1, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_2, "matrix"), as(erd1_2, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_3, "matrix"), as(erd1_3, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_4, "matrix"), as(erd1_4, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_5, "matrix"), as(erd1_5, "matrix"), tolerance=1e-7)

    ## col subsetting
    erd1_1 <- NULL
    erd1_2 <- NULL
    erd1_3 <- NULL
    erd1_4 <- NULL
    erd1_5 <- NULL
    erd1_1 <- erd1[,1]
    erd1_2 <- erd1[,2:3]
    erd1_3 <- erd1[,c(1,3,5)]
    erd1_4 <- erd1[,-c(1,3,5)]
    erd1_5 <- erd1[,c(colnames(erd1)[2], colnames(erd1)[4])]

    ers1_1 <- NULL
    ers1_2 <- NULL
    ers1_3 <- NULL
    ers1_4 <- NULL
    ers1_5 <- NULL
    ers1_1 <- ers1[,1]
    ers1_2 <- ers1[,2:3]
    ers1_3 <- ers1[,c(1,3,5)]
    ers1_4 <- ers1[,-c(1,3,5)]
    ers1_5 <- ers1[,c(colnames(erd1)[2], colnames(erd1)[4])]

    if(!silent)
    {
        print(erd1_1)
        print(ers1_1)
        print(erd1_2)
        print(ers1_2)
        print(erd1_3)
        print(ers1_3)
        print(erd1_4)
        print(ers1_4)
        print(erd1_5)
        print(ers1_5)
    }

    checkEquals(as(ers1_1, "matrix"), as(erd1_1, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_2, "matrix"), as(erd1_2, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_3, "matrix"), as(erd1_3, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_4, "matrix"), as(erd1_4, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_5, "matrix"), as(erd1_5, "matrix"), tolerance=1e-7)

    ## row and col subsetting
    erd1_1 <- NULL
    erd1_2 <- NULL
    erd1_3 <- NULL
    erd1_4 <- NULL
    erd1_5 <- NULL
    erd1_1 <- erd1[2:3,1]
    erd1_2 <- erd1[1,2:3]
    erd1_3 <- erd1[2:3,c(1,3,5)]
    erd1_4 <- erd1[c(1,3),-c(1,3,5)]
    erd1_5 <- erd1[c(names(x)[2], names(x[4])),c(colnames(erd1)[2], 
                                                 colnames(erd1)[4])]

    ers1_1 <- NULL
    ers1_2 <- NULL
    ers1_3 <- NULL
    ers1_4 <- NULL
    ers1_5 <- NULL
    ers1_1 <- ers1[2:3,1]
    ers1_2 <- ers1[1,2:3]
    ers1_3 <- ers1[2:3,c(1,3,5)]
    ers1_4 <- ers1[c(1,3),-c(1,3,5)]
    ers1_5 <- ers1[c(names(x)[2], names(x[4])),c(colnames(erd1)[2], 
                                                 colnames(erd1)[4])]

    if(!silent)
    {
        print(erd1_1)
        print(ers1_1)
        print(erd1_2)
        print(ers1_2)
        print(erd1_3)
        print(ers1_3)
        print(erd1_4)
        print(ers1_4)
        print(erd1_5)
        print(ers1_5)
    }

    checkEquals(as(ers1_1, "matrix"), as(erd1_1, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_2, "matrix"), as(erd1_2, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_3, "matrix"), as(erd1_3, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_4, "matrix"), as(erd1_4, "matrix"), tolerance=1e-7)
    checkEquals(as(ers1_5, "matrix"), as(erd1_5, "matrix"), tolerance=1e-7)
}

testSubsettingOfSparseERAgainstDenseERAll <- function(xc, kernel, silent,
                                                      annCharset=NULL,
                                                      ann1=NULL, ann2=NULL,
                                                      kernel2=NULL)
{
    dss <- DNAStringSet(chartr("Uu", "Tt", xc));
    names(dss) <- names(xc)

    if (is(kernel, "MotifKernel"))
        ass <- dss
    else
        ass <- translate(dss)
    
    for (i in 1:6)
    {
        currKernel <- kernel
        switch(i,
               "1" = {x <- dss},
               "2" = {x <- DNAVector(chartr("Uu", "Tt", xc)); 
                      names(x) <- names(xc)},
               "3" = {x <- RNAStringSet(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc);
                      if (!is.null(kernel2)) currKernel <- kernel2},
               "4" = {x <- RNAVector(chartr("Tt", "Uu", xc));
                      names(x) <- names(xc);
                      if (!is.null(kernel2)) currKernel <- kernel2},
               "5" = {x <- ass; names(x) <- names(xc)},
               "6" = {x <- AAVector(as.character(ass));
                      names(x) <- names(xc)}
               )

        testSubsettingOfSparseERAgainstDenseER(x=x, kernel=currKernel,
            silent=silent, annCharset=annCharset, ann1=ann1, ann2=ann2)
    }
}

