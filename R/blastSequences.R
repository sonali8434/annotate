.blastSequencesToDNAMultipleAlignment <- function(xml) {
   loadNamespace("Biostrings")
   loadNamespace("IRanges")
   qseq <- xpathSApply(xml, "//Hsp_qseq", xmlValue)
   hseq <- xpathSApply(xml, "//Hsp_hseq", xmlValue)
   res <- vector("list", length(qseq))
   for(i in seq_along(qseq)){
     res[[i]] <- Biostrings::DNAMultipleAlignment(
         c(hseq[[i]],qseq[[i]]),
         rowmask=as(IRanges::IRanges(), "NormalIRanges"),
         colmask=as(IRanges::IRanges(), "NormalIRanges"))
   }
   res
}

.blastSequencesToDataFrame <- function(xml) {
    if (xpathSApply(xml, "count(//Hit)") == 0L) {
        message("'blastSequences' returned 0 matches")
        return(data.frame())
    }

    iter <- xml["//Iteration"]
    iterlen <- sapply(iter, xpathSApply, "count(.//Hsp)")
    iterdf <- xmlToDataFrame(iter, stringsAsFactors=FALSE)

    hit <- xml["//Hit"]
    hitlen <- sapply(hit, xpathSApply, "count(.//Hsp)")
    hitdf <- xmlToDataFrame(hit, stringsAsFactors=FALSE)
    hitdf <- hitdf[, names(hitdf) != "Hit_hsps", drop=FALSE]

    hsp <- xmlToDataFrame(xml["//Hsp"] , stringsAsFactors=FALSE)

    df <- cbind(
        iterdf[rep(seq_len(nrow(iterdf)), iterlen),, drop=FALSE],
        hitdf[rep(seq_len(nrow(hitdf)), hitlen),, drop=FALSE],
        hsp)
    rownames(df) <- NULL
    df
}

.tryParseResult <- function(baseUrl, rid, rtoe, timeout) {
    message("estimated response time ", rtoe, " seconds")
    start <- Sys.time()
    end <- Sys.time() + timeout
    url <- sprintf("%s?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s",
                   baseUrl, rid)
    Sys.sleep(min(rtoe, timeout))
    repeat {
        elapsed <- as.double(Sys.time() - start, units="secs")
        ## RCurl::getURL(url, followlocation=TRUE) has issues.
        ## See getURL2() in R/query.R
        result <- as(htmlParse(getURL2(url),
                               error = xmlErrorCumulator(immediate=FALSE)),
                     "character")

        if (grepl("Status=FAILED", result))
            stop("BLAST search failed")
        else if  (grepl("Status=UNKNOWN", result))
            stop("BLAST search expired")
        else if (grepl("Status=READY", result)) {
            url <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
            ## RCurl::getURL(url, followlocation=TRUE) has issues.
            ## See getURL2() in R/query.R
            result <- xmlParse(getURL2(url),
                               error = xmlErrorCumulator(immediate=FALSE))
            return(result)
        } else if (grepl("Status=WAITING", result)) {
            message(sprintf("elapsed time %.0f seconds", elapsed))
            if (Sys.time() > end && interactive()) {
                msg <- sprintf("wait another %d seconds? [y/n] ", timeout)
                repeat {
                    ans <- substr(trimws(tolower(readline(msg))), 1, 1)
                    if (ans %in% c("y", "n"))
                        break
                }
                if (ans == "n")
                    break
                end <- Sys.time() + timeout
            }
            Sys.sleep(10)
        } else
            stop("BLAST search unknown response") 
    }
    msg <- sprintf("'blastSequences' timeout after %.0f seconds",
                   elapsed)
    stop(msg, call.=FALSE)
}

## Using the REST-ish API described at
## http://www.ncbi.nlm.nih.gov/blast/Doc/node2.html
blastSequences <- function(x, database="nr",
                           hitListSize="10",
                           filter="L",
                           expect="10",
                           program="blastn",
                           timeout=40,
                           as=c("DNAMultipleAlignment", "data.frame", "XML"))
{
    PARSE <- switch(match.arg(as),
                    DNAMultipleAlignment=.blastSequencesToDNAMultipleAlignment,
                    data.frame=.blastSequencesToDataFrame,
                    XML=identity)
    ## TODO: lots of argument checking and testing.  Also,
    ## depending on which program string is used we need to make the correct
    ## kind of object at the end (so blastn means DNAMultipleAlignment, and
    ## blastp means AAMultipleAlignment etc.

    ## So:
    ## 1) get online values these parameters can be
    ## 2) document those
    ## 3) restrict their vals in the code here.
    ## 4) for program, use this to determine what object is returned.
    
    ## assemble the query
    baseUrl <- "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
    query <- paste("QUERY=", URLencode(as.character(x)), "&DATABASE=",database,
                   "&HITLIST_SIZE=",hitListSize,"&FILTER=",filter,
                   "&EXPECT=",expect,"&PROGRAM=",program, sep="")
    url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
    ## RCurl::getURL(url, followlocation=TRUE) has issues.
    ## See getURL2() in R/query.R
    post <- htmlParse(getURL2(url0))
    
    x <- post[['string(//comment()[contains(., "QBlastInfoBegin")])']]
    rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
    rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
    result <- .tryParseResult(baseUrl, rid, rtoe, timeout)
    PARSE(result)
}

## took 11.5 minutes to do a blast...  (ugh)
