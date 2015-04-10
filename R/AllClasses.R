##345678901234567890123456789012345678901234567890123456789012345678901234567890

###################################################
##
## Sequence kernel classes
##
###################################################

#' @rdname SequenceKernel-class
#' @name SequenceKernel-class
#' @aliases
#' class:SequenceKernel
#' SequenceKernel-class
#' SequenceKernel
#' @title Sequence Kernel Class
#'
#' @details
#' This class represents the parent class for all sequence kernels. It is an
#' abstract class and must not be instantiated.
#'
#' @slot .Data the kernel function is stored in this slot. It is executed
#'       when the kernel matrix is created through invoking the kernel object.
#' @slot .userDefKernel indicates whether kernel is user defined or not.
#        The boolean is set by default and reset by the \code{kebabs} kernels.
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("SequenceKernel",
    representation = representation
    (
        .Data           = "function",
        .userDefKernel  = "logical"
    ),
    prototype = prototype
    (
        .userDefKernel  = TRUE
    )
)

#' @rdname SpectrumKernel-class
#' @name SpectrumKernel-class
#' @aliases
#' class:SpectrumKernel
#' SpectrumKernel-class
#' SpectrumKernel
#' @title Spectrum Kernel Class
#'
#' @details
#' Instances of this class represent a kernel object for the spectrum
#' kernel. The class is derived from \link{SequenceKernel}.
#'
#' @slot k length of the substrings considered by the kernel
#' @slot r exponent (for details see \link{spectrumKernel})
#' @slot annSpec when set the kernel evaluates annotation information
#' @slot distWeight distance weighting function or vector
#' @slot normalized data generated with this kernel object is normalized
#' @slot exact use exact character set for evaluation
#' @slot ignoreLower ignore lower case characters in the sequence
#' @slot presence consider only the presence of kmers not their counts
#' @slot revComplement consider a kmer and its reverse complement
#'       as the same feature
#' @slot mixCoef mixing coefficients for mixture kernel
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("SpectrumKernel",
    representation = representation
    (
        k               = "numeric",
        r               = "numeric",
        annSpec         = "logical",
        distWeight      = "ANY",
        normalized      = "logical",
        exact           = "logical",
        ignoreLower     = "logical",
        presence        = "logical",
        revComplement   = "logical",
        mixCoef         = "numeric"
    ),
    prototype = prototype
    (
        k               = 3,
        r               = 1,
        annSpec         = FALSE,
        distWeight      = numeric(0),
        normalized      = TRUE,
        exact           = TRUE,
        ignoreLower     = TRUE,
        presence        = FALSE,
        revComplement   = FALSE,
        mixCoef         = numeric(0)
    ),
    contains="SequenceKernel"
)

#' @rdname MismatchKernel-class
#' @name MismatchKernel-class
#' @aliases
#' class:MismatchKernel
#' MismatchKernel-class
#' MismatchKernel
#' @title Mismatch Kernel Class
#'
#' @details
#' Instances of this class represent a kernel object for the mismatch
#' kernel. The class is derived from \link{SequenceKernel}.
#'
#' @slot k length of the substrings considered by the kernel
#' @slot m maximum number of mismatches
#' @slot r exponent (for details see \link{mismatchKernel})
#' @slot annSpec not used for mismatch kernel
#' @slot distWeight not used for mismatch kernel
#' @slot normalized data generated with this kernel object is normalized
#' @slot exact use exact character set for evaluation
#' @slot ignoreLower ignore lower case characters in the sequence
#' @slot presence consider only the presence of kmers not their counts
#' @slot revComplement not used for mismatch kernel
#' @slot mixCoef not used for mismatch kernel
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("MismatchKernel",
    representation=representation
    (
        k               = "numeric",
        m               = "numeric",
        r               = "numeric",
        annSpec         = "logical",
        distWeight      = "ANY",
        normalized      = "logical",
        exact           = "logical",
        ignoreLower     = "logical",
        presence        = "logical",
        revComplement   = "logical"
    ),
    prototype = prototype
    (
        k               = 3,
        m               = 1,
        r               = 1,
        annSpec         = FALSE,
        distWeight      = numeric(0),
        normalized      = TRUE,
        exact           = TRUE,
        ignoreLower     = TRUE,
        presence        = FALSE,
        revComplement   = FALSE
    ),
    contains="SequenceKernel"
)

#' @rdname MotifKernel-class
#' @name MotifKernel-class
#' @aliases
#' class:MotifKernel
#' MotifKernel-class
#' MotifKernel
#' @title Motif Kernel Class
#'
#' @details
#' Instances of this class represent a kernel object for the motif
#' kernel. The class is derived from \link{SequenceKernel}. The motif
#' character vector is not stored in the kernel object.
#'
#' @slot r exponent (for details see \link{motifKernel})
#' @slot annSpec when set the kernel evaluates annotation information
#' @slot distWeight distance weighting function or vector
#' @slot normalized data generated with this kernel object is normalized
#' @slot exact use exact character set for evaluation
#' @slot ignoreLower ignore lower case characters in the sequence
#' @slot presence consider only the presence of motifs not their counts
#' @slot revComplement consider a kmer and its reverse complement
#'       as the same feature
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("MotifKernel",
    representation = representation
    (
        ## motifs stored only in environment
        r               = "numeric",
        annSpec         = "logical",
        distWeight      = "ANY",
        normalized      = "logical",
        exact           = "logical",
        ignoreLower     = "logical",
        presence        = "logical",
        revComplement   = "logical"
    ),
    prototype = prototype
    (
        r               = 1,
        annSpec         = FALSE,
        distWeight      = numeric(0),
        normalized      = TRUE,
        exact           = TRUE,
        ignoreLower     = TRUE,
        presence        = FALSE,
        revComplement   = FALSE
    ),
    contains="SequenceKernel"
)

#' @rdname GappyPairKernel-class
#' @name GappyPairKernel-class
#' @aliases
#' class:GappyPairKernel
#' GappyPairKernel-class
#' GappyPairKernel
#' @title Gappy Pair Kernel Class
#'
#' @details
#' Instances of this class represent a kernel object for the gappy pair
#' kernel. The kernel considers adjacent pairs of kmers with up to m
#' irrelevant characters between the pair. The class is derived from
#' \link{SequenceKernel}.
#'
#' @slot k length of the substrings considered by the kernel
#' @slot m maximum number of irrelevant character between two kmers
#' @slot r exponent (for details see \link{gappyPairKernel})
#' @slot annSpec when set the kernel evaluates annotation information
#' @slot distWeight distance weighting function or vector
#' @slot normalized data generated with this kernel object is normalized
#' @slot exact use exact character set for evaluation
#' @slot ignoreLower ignore lower case characters in the sequence
#' @slot presence consider only the presence of kmers not their counts
#' @slot revComplement consider a kmer and its reverse complement
#'       as the same feature
#' @slot mixCoef mixing coefficients for mixture kernel
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("GappyPairKernel",
    representation = representation
    (
        k               = "numeric",
        m               = "numeric",
        r               = "numeric",
        annSpec         = "logical",
        distWeight      = "ANY",
        normalized      = "logical",
        exact           = "logical",
        ignoreLower     = "logical",
        presence        = "logical",
        revComplement   = "logical",
        mixCoef         = "numeric"
    ),
    prototype = prototype
    (
        k               = 3,
        m               = 1,
        r               = 1,
        annSpec         = FALSE,
        distWeight      = numeric(0),
        normalized      = TRUE,
        exact           = TRUE,
        ignoreLower     = TRUE,
        presence        = FALSE,
        revComplement   = FALSE,
        mixCoef         = numeric(0)
    ),
    contains="SequenceKernel"
)

#' @rdname SymmetricPairKernel-class
#' @name SymmetricPairKernel-class
#' @aliases
#' class:SymmetricPairKernel
#' SymmetricPairKernel-class
#' SymmetricPairKernel
#' @title Symmetric Pair Kernel Class
#'
#' @details
#' Instances of this class represent a kernel object for the symmetric pair
#' kernel. The kernel does not compute similarity between single samples but
#' between two pairs of samples based on a regular sequence kernel for single
#' samples. The class is derived from \link{SequenceKernel}.
#'
#' @slot siKernel single instance kernel
#' @slot kernelType type of pair kernel
#' @slot r exponent (for details see \link{gappyPairKernel})
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("SymmetricPairKernel",
    representation = representation
    (
        siKernel        = "ANY",
        kernelType      = "character",
        r               = "numeric"
    ),
    prototype = prototype
    (
        siKernel        = NULL,
        kernelType      = "arithMean",
        r               = 1
    ),
    contains="SequenceKernel"
)


###################################################
##
## BioVector kernel classes
##
###################################################

## imported classes from IRanges
## setClassUnion("DataTableORNULL", c("NULL", "DataTable"))
## setClassUnion("characterORNULL", c("NULL", "character"))

#' @rdname BioVector-class
#' @title BioVector, DNAVector, RNAVector and AAVector Classes
#' @name BioVector-class
#' @aliases
#' class:BioVector
#' BioVector-class
#'
#' @details
#' This class is the parent class for representing sets of biological
#' sequences with support of lowercase characters. The derived classes
#' \code{\linkS4class{DNAVector}}, \code{\linkS4class{RNAVector}} and
#' \code{\linkS4class{AAVector}} hold DNA-, RNA- or AA-sequences which can
#' contain also lowercase characters. In many cases repeat regions are coded
#' as lowercase characters and with the \code{\linkS4class{BioVector}} based
#' classes sequence analysis with and without repeat regions can be performed
#' from the same sequence set. Whenever lowercase is not needed please use the
#' \code{\linkS4class{XStringSet}} based classes as they provide much richer
#' functionality. The class \code{BioVector }is derived from "character" and
#' holds the sequence information as character vector. Interfaces for the
#' small set of functions needed in KeBABS are designed consistent with
#' \code{\linkS4class{XStringSet}}.
#'
#' @slot NAMES sequence names
#' @slot elementMetadata element metadata, which is applicable per element
#'       and holds a DataTable with one entry per sequence in each column.
#'       KeBABS uses the column names "annotation" and "offset".
#' @slot metadata metadata applicable for the entire sequence set as list.
#'       KeBABS stores the annotation character set as list element
#'       named "annotationCharset".
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("BioVector",
    representation=representation
    (
        NAMES           = "characterORNULL",
        elementMetadata = "DataTableORNULL",
        metadata        = "list"
    ),
    prototype = prototype
    (
        NAMES           = NULL,
        elementMetadata = NULL,
        metadata        = list()
    ),
    contains="character"
)

#' @rdname BioVector-class
#' @name DNAVector-class
#' @aliases
#' class:DNAVector
#' DNAVector-class
#'
#' @details
#' Instances of the \code{DNAVector} class are used for representing sets of
#' DNA sequences.
#'

setClass("DNAVector",
    contains="BioVector"
)

#' @rdname BioVector-class
#' @name RNAVector-class
#' @aliases
#' class:RNAVector
#' RNAVector-class
#'
#' @details
#' Instances of the \code{RNAVector} class are used for representing sets of
#' RNA sequences.
#'

setClass("RNAVector",
    contains="BioVector"
)

#' @rdname BioVector-class
#' @name AAVector-class
#' @aliases
#' class:AAVector
#' AAVector-class
#'
#' @details
#' Instances of the \code{AAVector} class are used for representing sets of
#' amino acid sequences.
#'

setClass("AAVector",
    contains="BioVector"
)


###################################################
##
##  Kernel matrix class
##
###################################################

#' @rdname KernelMatrix-class
#' @name KernelMatrix-class
#' @aliases
#' class:KernelMatrix
#' KernelMatrix-class
#' KernelMatrix
#' @title Kernel Matrix Class
#'
#' @details
#' Instances of this class are used in KeBABS for storing a kernel matrix.
#' The hidden data part ".Data" contains the matrix.
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("KernelMatrix",representation("matrix"),
         prototype=structure(.Data=matrix()))

###################################################
##
##  Explicit representation classes
##
###################################################

#' @rdname ExplicitRepresentation-class
#' @name ExplicitRepresentation
#' @aliases
#' class:ExplicitRepresentation
#' ExplicitRepresentation-class
#' ExplicitRepresentation
#' @title Explicit Representation Dense and Sparse Classes
#'
#' @details
#' In KeBABS this class is the virtual parent class for explicit representations
#' generated from a set of biological sequences for a given kernel. The derived
#' classes \code{\linkS4class{ExplicitRepresentationDense}} and
#' \code{\linkS4class{ExplicitRepresentationSparse}} are meant to hold explicit
#' representations in dense or sparse format. The kernel used to generate the
#' explicit representation is stored together with the data.
#'
#' @slot usedKernel kernel used for generating the explicit representation
#' @slot quadratic boolean indicating a quadratic explicit representation
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("ExplicitRepresentation",
    representation = representation
    (
        usedKernel        = "ANY",
        quadratic          = "logical",
        "VIRTUAL"
    ),
    prototype = prototype
    (
        usedKernel        = NULL,
        quadratic          = FALSE
    )
)

#' @rdname ExplicitRepresentation-class
#' @aliases
#' class:ExplicitRepresentationDense
#' ExplicitRepresentationDense-class
#' ExplicitRepresentationDense
#'
#' @details
#' Instances of this class are used for storing explicit representations
#' in dense matrix format. This class is derived from
#' \code{\linkS4class{ExplicitRepresentation}}.


setClass("ExplicitRepresentationDense",
         contains = c("ExplicitRepresentation", "matrix"))

#' @rdname ExplicitRepresentation-class
#' @aliases
#' class:ExplicitRepresentationSparse
#' ExplicitRepresentationSparse-class
#' ExplicitRepresentationSparse
#'
#' @details
#' Instances of this class are used for storing explicit representations
#' in sparse \code{\linkS4class{dgRMatrix}} format. This class is derived from
#' \code{\linkS4class{ExplicitRepresentation}}.


setClass("ExplicitRepresentationSparse",
         contains = c("ExplicitRepresentation", "dgRMatrix"))


###################################################
##
##  Prediction profile class
##
###################################################

#' @rdname PredictionProfile-class
#' @name PredictionProfile-class
#' @aliases
#' class:PredictionProfile
#' PredictionProfile-class
#' PredictionProfile
#' @title Prediction Profile Class
#'
#' @details
#' This class stores prediction profiles generated for a set of biological
#' sequences from a trained model. Prediction profiles show the relevance
#' of individual sequence positions for the prediction result.
#'
#' @slot sequences sequence information for the samples with profiles
#' @slot baselines baselines generated from the offset in the model spread\cr
#'       to all sequence positions
#' @slot profiles prediction profile information stored as dense matrix with\cr
#'       the rows as samples and the columns as positions in the sample
#' @slot kernel kernel used for training the model on which these prediction\cr
#'       profiles are based
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("PredictionProfile",
    representation=representation
    (
        sequences         = "ANY",
        baselines         = "numeric",
        profiles         = "ANY",
        kernel           = "ANY"
    ),
    prototype=prototype
    (
        sequences         = NULL,
        baselines         = 0,
        profiles         = NULL,
        kernel           = NULL
    )
)


###################################################
##
##  ROCData class
##
###################################################

#' @rdname ROCData-class
#' @name ROCData-class
#' @aliases
#' class:ROCData
#' ROCData-class
#' ROCData
#' @title ROC Data Class
#'
#' @details
#' This class stores receiver operating characteristics (ROC) data.
#'
#' @slot AUC area under ROC curve
#' @slot TPR true positive rate for varying threshold
#' @slot FPR false positive rate for varying threshold
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("ROCData",
    representation=representation
    (
        AUC              = "numeric",
        TPR              = "numeric",
        FPR              = "numeric"
    ),
    prototype=prototype
    (
        AUC              = numeric(0),
        TPR              = numeric(0),
        FPR              = numeric(0)
    )
)


###################################################
##
##  kebabs model classes
##
###################################################

#' @rdname SVMInformation-class
#' @name SVMInformation-class
#' @aliases
#' class:SVMInformation
#' SVMInformation-class
#' SVMInformation
#' @title SVM Information Class
#'
#' @details
#' Instances of this class store SVM related information.
#'
#' @slot availPackages installed SVM packages
#' @slot reqSVM user requested SVM implementation
#' @slot reqPackage user requested package
#' @slot reqSVMPar user requested SVM parameters
#' @slot reqKernel user requested kernel
#' @slot reqExplicit user requested indictor of expl. rep. processing
#' @slot reqExplicitType user requested expl. rep. type
#' @slot reqFeatureType user requested feature type
#' @slot selSVM selected SVM implementation
#' @slot selPackage selected package
#' @slot selSVMPar selected SVM parameters
#' @slot selKernel selected kernel
#' @slot selExplicit selected indictor of expl. rep. processing
#' @slot explicitKernel kernel for explicit representation
#' @slot featureWeights indicator for feature weights
#' @slot weightLimit cutoff value for feature weights
#' @slot probModel indicator for probability model
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("SVMInformation",
    representation = representation
    (
        availPackages   = "character",
        reqSVM          = "character",
        reqPackage      = "character",
        reqSVMPar       = "list",
        reqKernel       = "ANY",
        reqExplicit     = "character",
        reqExplicitType = "character",
        reqFeatureType  = "character",
        selSVM          = "character",
        selPackage      = "character",
        selSVMPar       = "list",
        selKernel       = "ANY",
        selExplicit     = "logical",
        explicitKernel  = "character",
        featureWeights  = "character",
        weightLimit     = "numeric",
        probModel       = "logical"
    ),
    prototype = prototype
    (
        availPackages   = character(0),
        reqSVM          = character(0),
        reqPackage      = character(0),
        reqSVMPar       = list(),
        reqKernel       = NULL,
        reqExplicit     = character(0),
        reqExplicitType = character(0),
        reqFeatureType  = "linear",
        selSVM          = character(0),
        selPackage      = character(0),
        selSVMPar       = list(),
        selKernel       = NULL,
        selExplicit     = FALSE,
        explicitKernel  = character(0),
        featureWeights  = "auto",
        weightLimit     = 0,
        probModel       = FALSE
    ),
)

#' @rdname ControlInformation-class
#' @name ControlInformation-class
#' @aliases
#' class:ControlInformation
#' ControlInformation-class
#' ControlInformation
#' @title KeBABS Control Information Class
#'
#' @details
#' Instances of this class store control information for the
#' KeBABS meta-SVM.
#'
#' @slot classification indicator for classification task
#' @slot multiclassType type of multiclass SVM
#' @slot featureWeights feature weights control information
#' @slot selMethod selected processing method
#' @slot onlyDense indicator that only dense processing can be performed
#' @slot sparse indicator for sparse processing
#' @slot runtimeWarning indicator for runtime warning
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("ControlInformation",
    representation = representation
    (
        classification   = "logical",
        multiclassType   = "character",
        featureWeights   = "logical",
        selMethod        = "character",
        onlyDense        = "logical",
        sparse           = "logical",
        runtimeWarning   = "logical"
    ),
    prototype = prototype
    (
        classification   = FALSE,
        multiclassType   = character(0),
        featureWeights   = FALSE,
        selMethod        = character(0),
        onlyDense        = FALSE,
        sparse           = FALSE,
        runtimeWarning   = TRUE
    ),
)

#' @rdname CrossValidationResult-class
#' @name CrossValidationResult-class
#' @aliases
#' class:CrossValidationResult
#' CrossValidationResult-class
#' CrossValidationResult
#' @title Cross Validation Result Class
#'
#' @details
#' Instances of this class store the result of cross validation.
#'
#' @slot cross number of folds for cross validation
#' @slot noCross number of CV runs
#' @slot groupBy group assignment of samples
#' @slot perfParameters collected performance parameters
#' @slot outerCV flag indicating outer CV
#' @slot folds folds used in CV
#' @slot cvError cross validation error
#' @slot foldErrors fold errors
#' @slot noSV number of support vectors
#' @slot ACC cross validation accuracy
#' @slot BACC cross validation balanced accuracy
#' @slot MCC cross validation Matthews correlation coefficient
#' @slot AUC cross validation area under the ROC curve
#' @slot foldACC fold accuracy
#' @slot foldBACC fold balanced accuracy
#' @slot foldMCC fold Matthews correlation coefficient
#' @slot foldAUC fold area under the ROC curve
#' @slot sumAlphas sum of alphas
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("CrossValidationResult",
    representation = representation
    (
        cross            = "numeric",
        noCross          = "numeric",
        groupBy          = "ANY",
        perfParameters   = "character",
        outerCV          = "logical",
        folds            = "list",
        cvError          = "numeric",
        foldErrors       = "matrix",
        ACC              = "numeric",
        BACC             = "numeric",
        MCC              = "numeric",
        AUC              = "numeric",
        foldACC          = "matrix",
        foldBACC         = "matrix",
        foldMCC          = "matrix",
        foldAUC          = "matrix",
        noSV             = "numeric",
        sumAlphas        = "numeric"
    ),
    prototype = prototype
    (
        cross            = numeric(0),
        noCross          = numeric(0),
        groupBy          = NULL,
        perfParameters   = character(0),
        outerCV          = FALSE,
        folds            = list(),
        cvError          = numeric(0),
        foldErrors       = matrix(0),
        ACC              = numeric(0),
        BACC             = numeric(0),
        MCC              = numeric(0),
        AUC              = numeric(0),
        foldACC          = matrix(0),
        foldBACC         = matrix(0),
        foldMCC          = matrix(0),
        foldAUC          = matrix(0),
        noSV             = numeric(0),
        sumAlphas        = numeric(0)
    ),
)

#' @rdname ModelSelectionResult-class
#' @name ModelSelectionResult-class
#' @aliases
#' class:ModelSelectionResult
#' ModelSelectionResult-class
#' ModelSelectionResult
#' @title Model Selection Result Class
#'
#' @details
#' Instances of this class store the result of grid search or
#' model selection.
#'
#' @slot cross number of folds for cross validation
#' @slot noCross number of CV runs
#' @slot groupBy group assignment of samples
#' @slot nestedCross number of folds for outer CV
#' @slot noNestedCross number of runs of outer CV
#' @slot perfParameters collected performance parameters
#' @slot perfObjective performance criterion for grid search / model selection
#' @slot gridRows rows in grid search (i.e. kernels)
#' @slot gridCols columns in grid search
#' @slot gridErrors grid errors
#' @slot gridACC grid accuracy
#' @slot gridBACC grid balanced accuracy
#' @slot gridMCC grid Matthews correlation coefficient
#' @slot gridAUC grid area under the ROC curve
#' @slot gridNoSV grid number of support vectors
#' @slot gridSumAlphas grid sum of alphas
#' @slot smallestCVError smallest CV error
#' @slot selGridRow grid row of best result
#' @slot selGridCol grid col of best result
#' @slot fullModel full model for best result
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("ModelSelectionResult",
    representation = representation
    (
        cross            = "numeric",
        noCross          = "numeric",
        groupBy          = "ANY",
        nestedCross      = "numeric",
        noNestedCross    = "numeric",
        perfParameters   = "character",
        perfObjective    = "character",
        gridRows         = "list",
        gridCols         = "ANY",
        gridErrors       = "ANY",
        gridACC          = "ANY",
        gridBACC         = "ANY",
        gridMCC          = "ANY",
        gridAUC          = "ANY",
        gridNoSV         = "ANY",
        gridSumAlphas    = "ANY",
        smallestCVError  = "numeric",
        selGridRow       = "ANY",
        selGridCol       = "ANY",
        fullModel        = "ANY"
    ),
    prototype = prototype
    (
        cross            = numeric(0),
        noCross          = numeric(0),
        groupBy          = NULL,
        nestedCross      = numeric(0),
        noNestedCross    = numeric(0),
        perfParameters   = character(0),
        perfObjective    = character(0),
        gridRows         = list(),
        gridCols         = NULL,
        gridErrors       = NULL,
        gridACC          = NULL,
        gridBACC         = NULL,
        gridMCC          = NULL,
        gridAUC          = NULL,
        gridNoSV         = NULL,
        gridSumAlphas    = NULL,
        smallestCVError  = Inf,
        selGridRow       = NULL,
        selGridCol       = NULL,
        fullModel        = NULL
    ),
)

#' @rdname KBModel-class
#' @name KBModel-class
#' @aliases
#' class:KBModel
#' KBModel-class
#' KBModel
#' @title KeBABS Model Class
#'
#' @details
#' Instances of this class represent a model object for the KeBABS
#' meta-SVM.
#'
#' @slot call invocation string of KeBABS meta-SVM
#' @slot numSequences number of sequences used for training
#' @slot sel index subset of samples used for training
#' @slot y vector of target values
#' @slot levels levels of target
#' @slot numClasses number of classes
#' @slot classNames class labels
#' @slot classWeights class weights
#' @slot SV support vectors
#' @slot svIndex support vector indices
#' @slot alphaIndex list of SVM indices per SVM
#' @slot trainingFeatures feature names used in training
#' @slot featureWeights feature Weights
#' @slot b model offset
#' @slot probA fitted logistic function parameter A
#' @slot probB fitted logistic function parameter A
#' @slot sigma scale of Laplacian fitted to regression residuals
#' @slot cvResult cross validation result of class
#'       \code{\linkS4class{CrossValidationResult}}
#' @slot modelSelResult model selection / grid search result of class
#'       \code{\linkS4class{ModelSelectionResult}}
#' @slot ctlInfo KeBABS control info of class
#'       \code{\linkS4class{ControlInformation}}
#' @slot svmInfo info about requested / used SVM of class
#'       \code{\linkS4class{SVMInformation}}
#' @slot svmModel original model returned from SVM
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setClass("KBModel",
    representation = representation
    (
        call             = "character",       # invocation string
        numSequences     = "numeric",         # number of sequences
        sel              = "numeric",         # selected subset of sequences
        y                = "ANY",             # target values
        levels           = "ANY",             # levels of target
        numClasses       = "numeric",         # number of classes
        classNames       = "ANY",             # class labels
        classWeights     = "numeric",         # class weights
        SV               = "ANY",             # support vectors
        svIndex          = "ANY",             # support vector indices
        alphaIndex       = "ANY",             # list of SV indices per SVM
        trainingFeatures = "characterORNULL", # feature names of training
        featureWeights   = "ANY",             # feature weights
        b                = "numeric",         # b (bias)
        probA            = "numeric",         # fitted logistic function param A
        probB            = "numeric",         # fitted logistic function param B
        sigma            = "numeric",         # for regression only - scale of
                                              # laplacian fitted to residuals
        cvResult         = "CrossValidationResult", # cross validation result
        modelSelResult   = "ANY",             # model selection result
        ctlInfo          = "ControlInformation",    # kb control info
        svmInfo          = "SVMInformation",  # info about requested/used svm
        svmModel         = "ANY"              # model object of ksvm, LiblineaR,
                                              #   libsvm via e1071, psvm,
                                              #   rvm, lssvm, penalizedSVM
    ),
    prototype = prototype
    (
        call             = character(0),
        numSequences     = 0,
        sel              = NULL,
        y                = NULL,
        levels           = NULL,
        numClasses       = 0,
        classNames       = NULL,
        classWeights     = numeric(0),
        SV               = NULL,
        svIndex          = numeric(0),
        alphaIndex       = numeric(0),
        trainingFeatures = character(0),
        featureWeights   = NULL,
        b                = 0,
        probA            = as.numeric(NA),
        probB            = as.numeric(NA),
        sigma            = as.numeric(NA),
        cvResult         = NULL,
        modelSelResult   = NULL,
        ctlInfo          = NULL,
        svmInfo          = NULL,
        svmModel         = NULL
    ),
)

setClass("KebabsInfo",
    representation=representation
    (
        supportedPkgs         = "character",
        allowedSeqClasses     = "character",
        allowedSeqSetClasses  = "character",
        allowedSeqTypes       = "numeric",
        classifierMap         = "matrix",
        svmParMap             = "matrix",
        modelParMap           = "matrix",
        kebabsDebug           = "logical"
    ),
    prototype=prototype
    (
        supportedPkgs         = character(0),
        allowedSeqClasses     = character(0),
        allowedSeqSetClasses  = character(0),
        allowedSeqTypes       = numeric(0),
        classifierMap         = matrix(0),
        svmParMap             = matrix(0),
        modelParMap           = matrix(0),
        kebabsDebug           = FALSE
    ),
)
