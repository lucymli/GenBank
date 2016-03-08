

#' download.fasta
#'
#' @return
#' @export
#'
#' @examples
download.fasta <- function(accessions, database="nucleotide", file.name=NULL,
                           return.seq=ifelse(is.null(file.name), TRUE, FALSE)) {
  search.url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=",
                       database, "&rettype=fasta&retmode=text&id=",
                       paste(accessions, collapse=","))
  if (is.null(file.name)) file.name <- tempfile()
  download.file(search.url, file.name)
  if (return.seq) {
    lines <- readLines(file.name)
    titles.pos <- grep(">", lines)
    sequence.starts <- titles.pos+1
    sequence.ends <- c(titles.pos[-1]-1, length(lines))
    sequences <- lapply(seq_along(titles.pos), function (i) {
      paste(lines[sequence.starts[i]:sequence.ends[i]], collapse="")
    })
    names(sequences) <- lines[titles.pos]
    return(sequences)
  } else {
    invisible(NULL)
  }
}

#' download.gb
#'
#' @return
#' @export
#'
#' @examples
download.gb <- function(accessions, database="nucleotide", file.name=NULL,
                        return.gb=ifelse(is.null(file.name), TRUE, FALSE)) {
  search.url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=",
                       database, "&rettype=gb&retmode=text&id=",
                       paste(accessions, collapse=","))
  if (is.null(file.name)) file.name <- tempfile()
  download.file(search.url, file.name)
  if (return.gb) {
    gb.lines <- readLines(file.name)
    gb <- list()
    gb$locus <- lapply(strsplit(grep("LOCUS", gb.lines, ignore.case=FALSE, value=TRUE), " "), function (x) {
      y <- x[nchar(x)>0]
      c(accession=y[2], length=y[3], seqtype=y[5], genometype=y[6], kingdom=y[7],
        submission.date=as.Date(y[8], format="%d-%b-%Y"))
    })
    gb$definition <- gsub("DEFINITION[ ]*", "", grep("DEFINITION", gb.lines, ignore.case=FALSE, value=TRUE))
    gb$accession <- gsub("ACCESSION[ ]*", "", grep("ACCESSION", gb.lines, ignore.case=FALSE, value=TRUE))
    gb$version <- lapply(strsplit(gsub("VERSION[ ]*", "", grep("VERSION", gb.lines, ignore.case=FALSE, value=TRUE)), " "), function (x) {
      c(accession=x[1], GI=gsub("GI:", "", tail(x, 1)))
    })
    gb$keywords <- gsub("KEYWORDS[ ]*", "", grep("KEYWORDS", gb.lines, ignore.case=FALSE, value=TRUE))
    gb$source <- gsub("SOURCE[ ]*", "", grep("SOURCE", gb.lines, ignore.case=FALSE, value=TRUE))
    gb$organism <- gsub("ORGANISM[ ]*", "", grep("ORGANISM", gb.lines, ignore.case=FALSE, value=TRUE))
    source.pos <- grep("source", gb.lines, ignore.case=FALSE)
    CDS.pos <- grep("CDS", gb.lines, ignore.case=FALSE)
    ORIGIN.pos <- grep("ORIGIN", gb.lines, ignore.case=FALSE)
    gb$features.source <- lapply(seq_along(source.pos), function (i) {
      raw.source.lines <- sub("^\\s+", "", gb.lines[(source.pos[i]+1):(CDS.pos[i]-1)])
      split.source.lines <- strsplit(raw.source.lines, "=")
      unlist(lapply(split.source.lines, function (x) {
        NAME <- substr(x[1], 2, nchar(x[1]))
        value <- c(substr(x[2], 2, nchar(x[2])-1))
        names(value) <- NAME
        return(value)
      }))
    })
    gb$features.CDS <- lapply(seq_along(CDS.pos), function (i) {
      raw.source.lines <- sub("^\\s+", "", gb.lines[(CDS.pos[i]+1):(ORIGIN.pos[i]-1)])
      split.source.lines <- strsplit(raw.source.lines, "=")
      unlist(lapply(split.source.lines, function (x) {
        NAME <- substr(x[1], 2, nchar(x[1]))
        value <- c(substr(x[2], 2, nchar(x[2])-1))
        names(value) <- NAME
        return(value)
      }))
    })
    return(gb)
  } else {
    invisible(NULL)
  }
}
