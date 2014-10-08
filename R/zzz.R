.onLoad <- function(lib, pkg, ...)
{
    #initKebabs(lib, pkg)
}

.onUnload <- function(lib)
{
    cleanupKebabs(lib)
}
