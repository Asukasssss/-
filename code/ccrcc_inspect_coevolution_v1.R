args <- commandArgs(TRUE)
e <- new.env()
load(args[1], envir=e)
for (nm in ls(e)) {
 x <- e[[nm]]
 cat('OBJECT',nm,'CLASS',class(x),'DIM',dim(x),'\n')
 if (is.data.frame(x) || is.matrix(x)) cat('COLUMNS',paste(head(colnames(x),30),collapse=' | '),'\n')
 if (is.list(x) && !is.data.frame(x)) cat('NAMES',paste(head(names(x),30),collapse=' | '),'\n')
}
