datacache <- new.env(hash=TRUE, parent=emptyenv())

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    genelocate <- system.file("extdata", "genelocate.txt", package=pkgname, lib.loc=libname)
    genelocate <- read.table(genelocate,header=TRUE,sep="\t",stringsAsFactors=FALSE)
    assign("genelocate", genelocate, envir=datacache)
    chromLength <- system.file("extdata", "chromLength.txt", package=pkgname, lib.loc=libname)
    chromLength <- read.table(chromLength,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    assign("chromLength", chromLength, envir=datacache)
}