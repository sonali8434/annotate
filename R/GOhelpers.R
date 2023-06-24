##Copyright R. Gentleman, 2004
##simple functions to get Evidence codes

.isMissingGOEntry <- function(x) (length(x) == 1L && is.na(x))

##get then GO term names for a particular (sub)ontology
getOntology = function(inlist, ontology=c("MF", "BP", "CC")) {
   which = match.arg(ontology)
   onts = sapply(inlist, function(z) {
       if (!.isMissingGOEntry(z))
         z$Ontology
       else
         z
       })
   onts = onts[!is.na(onts)]
   unique(names(inlist[onts %in% which]))
}


##get GO evidence codes
getEvidence = function(inlist) {
    ans <- sapply(inlist, function(z) {
         if (!.isMissingGOEntry(z))
           z$Evidence
         else
           z
     })
    ans[!is.na(ans)]
}


##drop a specified set of evidence codes
dropECode = function(inlist, code = "IEA") {
    hasCode = sapply(inlist, function(z) {
        if (!.isMissingGOEntry(z))
          z$Evidence
        else
          z
        })
    hasCode <- hasCode[!is.na(hasCode)]
    badVals = hasCode %in% code
    inlist[!badVals]
}


## helper function, determines if there is a GO annotation for the
## desired mode
hasGOannote <- function(x, which="MF") {
    if (is(x, "GOTerms")) {
        cat <- Ontology(x)
        if (!is.na(cat) && cat == which)
          return(TRUE) else return(FALSE)
    }
    if (is.list(x)) {
        gT <- sapply(x, function(y) is(y, "GOTerms"))
        if (any(gT)) {
            if (all(gT)) {
                cats <- sapply(x, Ontology)
                return(cats == which)
            }
            else
              stop("mixed arguments not allowed")
        }
    }
    if (!is.character(x))
      stop("wrong argument")
    tm <- getGOOntology(x)
    return(tm == which)
}


##three functions to get all the GO information for a set of GO terms
##FIXME: these need to be renovated - probably removed even..
 getGOOntology <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return( character(0))
     loadNamespace("GO.db")
     wh <- mget(x, envir=GO.db::GOTERM, ifnotfound=NA)
     return( sapply(wh, Ontology) )
 }

 getGOParents <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     loadNamespace("GO.db")
     MF_parents <- mget(x, envir=GO.db::GOMFPARENTS, ifnotfound=NA)
     BP_parents <- mget(x, envir=GO.db::GOBPPARENTS, ifnotfound=NA)
     CC_parents <- mget(x, envir=GO.db::GOCCPARENTS, ifnotfound=NA)
     lapply(setNames(seq_along(x), x),
         function(i) {
             xi_parents <- MF_parents[[i]]
             if (!identical(xi_parents, NA))
                 return(list(Ontology="MF", Parents=xi_parents))
             xi_parents <- BP_parents[[i]]
             if (!identical(xi_parents, NA))
                 return(list(Ontology="BP", Parents=xi_parents))
             xi_parents <- CC_parents[[i]]
             if (!identical(xi_parents, NA))
                 return(list(Ontology="CC", Parents=xi_parents))
             stop(paste(x[[i]], "is not a member of any ontology"))
         }
     )
 }

 getGOChildren <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     loadNamespace("GO.db")
     MF_children <- mget(x, envir=GO.db::GOMFCHILDREN, ifnotfound=NA)
     BP_children <- mget(x, envir=GO.db::GOBPCHILDREN, ifnotfound=NA)
     CC_children <- mget(x, envir=GO.db::GOCCCHILDREN, ifnotfound=NA)
     lapply(setNames(seq_along(x), x),
         function(i) {
             xi_children <- MF_children[[i]]
             if (!identical(xi_children, NA))
                 return(list(Ontology="MF", Children=xi_children))
             xi_children <- BP_children[[i]]
             if (!identical(xi_children, NA))
                 return(list(Ontology="BP", Children=xi_children))
             xi_children <- CC_children[[i]]
             if (!identical(xi_children, NA))
                 return(list(Ontology="CC", Children=xi_children))
             list()  # not an error (unlike for getGOParents() above)
         }
     )
 }

 getGOTerm <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     loadNamespace("GO.db")
     terms <- mget(x, envir=GO.db::GOTERM, ifnotfound=NA)
     isNA = sapply(terms,function(x) !(isS4(x) && is(x, "GOTerms")))
     if( any(isNA) )
         terms = terms[!isNA]

     ontology <- sapply(terms, Ontology)
     terms = sapply(terms, Term)
     return(split(terms, ontology))
 }


filterGOByOntology <- function(goids, ontology=c("BP", "CC", "MF")) {
    ontology <- match.arg(ontology)
    eName <- switch(ontology,
                    BP="GOBPPARENTS",
                    CC="GOCCPARENTS",
                    MF="GOMFPARENTS",
                    stop("invalid ontology ", ontology))
    e <- get(eName)
    goids %in% ls(e)
}

aqListGOIDs <- function(ont) {
    ## Return all GO IDs in the specified ontologies
    ont <- unique(ont)
    knownOnts <- c("BP", "CC", "MF")
    badOnt <- ont[!(ont %in% knownOnts)]
    if (length(badOnt))
      stop("Unknown ontology codes: ", paste(badOnt, collapse=", "),
           "\nvalid codes are: ", paste(knownOnts, collapse=", "))
    ## determine size
    lens <- integer(length(ont))
    for (i in seq(along=ont))
      lens[i] <- length(getAnnMap(paste(ont[i], "PARENTS", sep=""),
                                  chip="GO"))
    ## retrieve IDs
    ans <- character(sum(lens))
    lens <- c(0L, lens)
    for (i in seq(along=ont)) {
        ans[lens[i]+1:lens[i+1]] <- ls(getAnnMap(paste(ont[i], "PARENTS", sep=""),
                                               chip="GO"))
    }
    ans
}
