\name{blastSequences}
\alias{blastSequences}

\title{
  Run a blast query to NCBI for either a string or an entrez gene ID and
  then return a series of MultipleAlignment objects.
}

\description{
  This function sends a query to NCBI as a string of sequence or an
  entrez gene ID and then returns a series of MultipleAlignment objects.
}

\usage{
  blastSequences(x, database, hitListSize, filter, expect, program,
      timeout=40, as=c("DNAMultipleAlignment", "data.frame", "XML"))
}

\arguments{
  \item{x}{
    A sequence as a character vector or an integer corresponding to an
    entrez gene ID. Submit multiple sequences as a length-1 character
    vector, \code{x = ">ID-1\nACATGCTA\n>ID-2\nAAACCACTT"}.
  }
  \item{database}{
    Which NCBI database to use. If not \dQuote{blastn}, then set
    \code{as="XML"}
  }
  \item{hitListSize}{
    Number of hits to keep.
  }
  \item{filter}{
    Sequence filter; \dQuote{L} for Low Complexity, \dQuote{R} for Human Repeats,
    \dQuote{m} for Mask lookup
  }
  \item{expect}{
    The BLAST \sQuote{expect} value above which matches will be
    returned.
  }
  \item{program}{
    Which program do you want to use for blast.
  }
  \item{timeout}{
    Approximate maximum length of time, in seconds, to wait for a result.
  }
  \item{as}{
    character(1) indicating whether the result from the NCBI server
    should be parsed to a list of \code{DNAMultipleAlignment} instances,
    represented as a \code{data.frame}, or returned as XML.
  }

}

\details{
  Right now the function only works for "blastn".

  The NCBI URL api used by this function is documented at
  \url{https://www.ncbi.nlm.nih.gov/blast/Doc/urlapi.html}
}

\value{

  By default, a series of \code{DNAMultipleAlignment} (see
  \code{\link[Biostrings]{MultipleAlignment-class}}
  objects. Alternatively, a \code{data.frame} or XML document returned
  from the NCBI server. The \code{data.frame} is a \sQuote{long form}
  representation of the \sQuote{Iteration}, \sQuote{Hit} and
  \sQuote{Hsp} results returned from the server. The XML document is the
  result of the \code{xmlParse} function of the XML library, and follows
  the format described by
  \url{https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd} and
  \url{https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd}. }

\author{
  M. Carlson
}

\examples{

## x can be an entrez gene ID
blastSequences(17702, timeout=40, as="data.frame")

if (interactive()) {

    ## or x can be a sequence
    blastSequences(x = "GGCCTTCATTTACCCAAAATG")

    ## hitListSize does not promise that you will get the number of
    ## matches you want..  It will just try to get that many.
    blastSequences(x = "GGCCTTCATTTACCCAAAATG", hitListSize="20")

}
}
