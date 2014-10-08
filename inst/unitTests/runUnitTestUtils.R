#2345678901234567890123456789012345678901234567890123456789012345678901234567890
testName <- function(silent, testName)
{
    if (!silent)
    {
        cat(paste("\n\n****************************************",
                  "***********************************\n", sep=""))
        cat("Test: ", testName, "\n")
        cat(paste("\n\n****************************************",
                  "***********************************\n", sep=""))
    }
}

testStepName <- function(silent, teststepName)
{
    if (!silent)
    {
        cat("\nTeststep:", teststepName, "\n")
        cat(paste("------------------------------------------",
                  "---------------------------------\n", sep=""))
    }
}

defMat <- function(mat, namesX, namesY)
{
    if (!missing(namesX))
    {
        checkEquals(length(namesX), dim(mat)[1])
        rownames(mat) <- namesX
    }
    
    if (!missing(namesY))
    {
        checkEquals(length(namesY), dim(mat)[2])
        colnames(mat) <- namesY
    }
    
    mat
}

normalizeMatrix <- function(km, kvX, kvY)
{
    if (missing(kvX) && missing(kvY))
    {  
        if (dim(km)[1] != dim(km)[2])
        {
            print(km)
            checkTrue(FALSE)
        }
        x1 <- sqrt(diag(km))
        y1 <- sqrt(diag(km))
    }
    else
    {
        if (length(kvX) != dim(km)[1] || length(kvY) != dim(km)[2])
        {
            print(km)
            checkTrue(FALSE)
        }
        x1 <- sqrt(kvX)
        y1 <- sqrt(kvY)
    }
    
    km <- km / outer(x1, y1)
    
    as.KernelMatrix(km)
}

normalizeERD <- function(erd)
{
    erdn <- t(apply(erd,1,function(x){x/sqrt(sum(x^2))}))
    
    if (length(rownames(erd)) > 0)
        rownames(erdn) <- rownames(erd)
    
    if (length(colnames(erd)) > 0)
        colnames(erdn) <- colnames(erd)
    
    erdn
}
