## S4 methods so we can use non-package based annotation databases (e.g., from AnnotationHub)
## as if they were installed packages

setGeneric("isValidKey", function(ids, pkg) standardGeneric("isValidKey"))

setGeneric("allValidKeys", function(pkg) standardGeneric("allValidKeys"))
