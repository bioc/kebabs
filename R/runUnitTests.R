#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname runUnitTests
#' @aliases runit
#' @title RUnit Tests for Kernels
#'
#' @description Start Command for RUnit Tests
#'
#' @param unit just for internal testing purpose.
#' @param echo just for internal testing purpose.
#' @param silent just for internal testing purpose.
#' @param seed just for internal testing purpose.
#' @param ... just for internal testing purpose.
#'
#' @details
#' This command is just for internal testing purpose and is exported only for 
#' simplifying developer testing.
#' @return
#' see above
#' @examples
#' runit("MotifKernel")
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @export

runit <- function(unit=NULL, echo=FALSE, silent=TRUE, seed=789, ...)
{
    source(paste(.libPaths()[1],"/kebabs/unitTests/runUnitTestUtils.R",
                 sep=""), echo=echo, ...)
    source(paste(.libPaths()[1],"/kebabs/unitTests/kernelTestCommon.R",
                 sep=""), echo=echo, ...)

    if (!is.null(unit))
    {
        filePath <- paste(.libPaths()[1],"/kebabs/unitTests/runit-",
                          unit, "UnitTest.R", sep="")

        if (!file.exists(filePath))
            stop("invalid RUnit test specified\n")

        ## run single unit test
        source(paste(.libPaths()[1],"/kebabs/unitTests/runit-",
                     unit, "UnitTest.R", sep=""), echo=echo, ...)
        cat(unit, " Result:\n\n")
        expr <- parse(text=paste("test_", unit, "(seed=", seed,
                                 ", silent=", silent, ")",sep=""))
        startTime <- proc.time()
        eval(expr)
        print(proc.time() - startTime)
        cat("\n")
    }
    else
    {
        ## run all unit tests
        for (nm in list.files(paste(.libPaths()[1],"/kebabs/unitTests/",
                                    sep="")))
        {
            if (substr(nm,1,6) == "runit-")
            {
                cat(nm,":\n")
                source(paste(.libPaths()[1],"/kebabs/unitTests/",
                             nm, sep=""), echo=echo, ...)
                unit <- substr(nm, 7, nchar(nm) - 10)
                cat(unit, " Result:\n\n")
                expr <- parse(text=paste("test_", unit, "(seed=", seed,
                                         ", silent=", silent, ")",sep=""))
                startTime <- proc.time()
                eval(expr)
                print(proc.time() - startTime)
                cat("\n")
            }
        }
    }
}
