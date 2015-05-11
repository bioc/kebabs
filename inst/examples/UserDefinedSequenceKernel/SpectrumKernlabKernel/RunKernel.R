library(kebabs)

## load user defined sequence kernel
source(paste(path.package("kebabs"),
             "/examples/UserDefinedSequenceKernel/SpectrumKernlabKernel/SpectrumKernlabKernel.R",
             sep=""))

## load data
data(TFBS)

## redefine sequences as BioVector because kernlab does not
## support Biostrings classes
## for kernels supporting the XStringSet classes this
## conversion is not necessary
seqs <- DNAVector(as.character(enhancerFB))

## assign short names to sequences
names(seqs) <- paste("S", 1:length(seqs), sep="_")

## show sequences
seqs

## define kernel object of user defined sequence kernel which
## in this example is the kernlab stringdot spectrum kernel
skk <- spectrumKernlabKernel(k=5)

## show user-defined kernel object
skk

## computation of kernel matrix with user defined sequence kernel
#################################################################

## generate kernel matrix for 100 sequences
km <- getKernelMatrix(skk, seqs[1:100])

km[1:6,1:6]

## training and prediction with user defined sequence kernel
## train only on 100 sequences to shorten time for
## computation of kernel matrix
#################################################################
train <- sample(1:length(seqs), 100)
test <- c(1:length(seqs))[-train]

model <- kbsvm(seqs[train], yFB[train], skk, pkg="kernlab",
               svm="C-svc", cost=20, featureWeights="no")
pred <- predict(model, seqs[test])
evaluatePrediction(pred, yFB[test], unique(yFB))

## cross validation
#################################################################
modelc <- kbsvm(seqs[train], yFB[train], skk, pkg="kernlab",
                svm="C-svc", cost=20, cross=5, featureWeights="no")
cvResult(modelc)

## grid search
#################################################################
modelg <- kbsvm(seqs[train], yFB[train], skk, pkg="kernlab",
                svm="C-svc", cost=c(1,2,3,5,10,20,50,100),
                cross=10, noCross=5, featureWeights="no",
                explicit="no", showProgress=TRUE)
modelSelResult(modelg)

## model selection
#################################################################
modelm <- kbsvm(seqs[train], yFB[train], skk, pkg="kernlab",
                svm="C-svc", cost=c(2,10,20,50,100),
                cross=10, nestedCross=5, featureWeights="no",
                explicit="no", showProgress=TRUE)
modelSelResult(modelm)
cvResult(modelm)

## no support for explicit representation, feature weights
## and prediction profiles for user-defined kernel
## user-defined kernel can not be used as single
## instance kernel in symmetric pair kernel
