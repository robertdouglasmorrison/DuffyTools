# table2html.R
#


`table2html` <- function( x, fileout, title="", maxRows=NULL, 
			extraHTMLtextFile=NULL, extraHTMLtext=NULL, 
			linkColumnNames="GENE_ID", linkPaths=".", linkExtensions=".png", 
			linkTargetColumnNames=linkColumnNames) {

	options("useUnquotedTagnames"=TRUE)
	if ( ! is.null( linkColumnNames)) {
		linkPaths <- rep( linkPaths, length.out=length(linkColumnNames))
		linkExtensions <- rep( linkExtensions, length.out=length(linkColumnNames))
		linkTargetColumnNames <- rep( linkTargetColumnNames, length.out=length(linkColumnNames))
	}

	# force upper limit on what gets written out
	if ( !is.null( maxRows) && nrow(x) > maxRows) {
		x <- x[1:maxRows, ]
	}

	# build the final result
	if (nrow(x)) {
		options( digits=5)
		out <- format( x)

		if ( ! (is.null( linkColumnNames) || any( is.na( linkColumnNames)))) {
	  		for ( j in 1:length( linkColumnNames)) {
				gTextName <- x[[ linkColumnNames[j] ]]
				gTargetName <- x[[ linkTargetColumnNames[j] ]]
				# quietly ignore if those names are not in the file
				if ( ! (linkColumnNames[j] %in% colnames(x))) next
				# we are making it easier to have a link say/show one thing (Name) 
				# but hyperlink contentID (Target) is from another column
				plotLinks <- makeGenePlotLinks( gTargetName, path=linkPaths[j], type="html", 
							extension=linkExtensions[j], gText=gTextName)
				out[[ linkColumnNames[j] ]] <- plotLinks
	  		}
		}
	} else {
		out <- x
	}

	writeAsWebPage( out, fileout, title=title, extraHTMLtextFile=extraHTMLtextFile,
			extraHTMLtext=extraHTMLtext)
	# cat("\nWrote HTML file:  ", fileout, "\n")
}


`htmlTitle` <- function( label="", sampleID=NULL, groupID=NULL, speciesID=NULL) {
	
	HTML_BR <- "<br>"
	HTML_SP <- "&nbsp;"

	out <- ""
	if ( ! is.null( sampleID)) {
		out <- paste( out, "SampleID =", HTML_SP, sampleID, HTML_BR)
	}
	if ( ! is.null( groupID)) {
		out <- paste( out, "GroupID =", HTML_SP, groupID, HTML_BR)
	}
	if ( ! is.null( speciesID)) {
		out <- paste( out, "SpeciesID =", HTML_SP, speciesID, HTML_BR)
	}
	out <- paste( out, label)
	out
}


`format4html` <- function( x, format=NULL, cleanColumns=NULL) {

	out <- x

	# format the values by column
	if ( ! is.null( format)) {
		for ( i in 1:min( ncol(x), length(format))) {
			fmt <- format[i]
			if ( is.null(fmt) || is.na(fmt) || fmt == "" || fmt == "c") next
			v <- x[[i]]
			comma <- if (regexpr( ',', fmt, fixed=T) > 0) "," else ""
			digits <- as.integer( substr( fmt, 2,2))
			if ( is.na(digits)) digits <- NULL
			fmt <- substr( tolower(fmt), 1,1)

			if (fmt %in% c("f","e")) {
				v <- as.numeric(v)
				v <- formatC( v, format=fmt, digits=digits, big.mark=comma)

			} else if ( fmt == "d") {
				v <- as.numeric(v)
				v <- formatC( round(v), format="d", big.mark=comma)
			}
			out[[i]] <- v
		}
	}

	# clean the column names
	if ( ! is.null( cleanColumns)) {
		toClean <- intersect( 1:ncol(x), cleanColumns)
		if ( length( toClean) > 0) {
			nams <- colnames(x)[ toClean]
			nams <- gsub( "_", " ", nams, fixed=T)
			nams <- gsub( ".", " ", nams, fixed=T)
			colnames(out)[ toClean] <- nams
		}
	}

	out
}


`addExcelLinksToTable` <- function( x, linkColumnNames="GENE_ID", linkPaths=".", linkExtensions=".png") {

	# options("useUnquotedTagnames"=TRUE)

	# build the final result
	options( digits=6)
	out <- format( x)

	if ( ! (is.null( linkColumnNames) || any( is.na( linkColumnNames)))) {
	  for ( j in 1:length( linkColumnNames)) {
		gNames <- x[[ linkColumnNames[j] ]]
		plotLinks <- makeGenePlotLinks( gNames, path=linkPaths[j], type="excel", extension=linkExtensions[j])
		out[[ linkColumnNames[j] ]] <- plotLinks
	  }
	}

	return( out)
}


makeGenePlotLinks <- function(gTarget, path="./png", type="excel", prefix="", suffix="", extension=".png",
				gText=gTarget) {

	# create text that excel will treat as a hyperlink to our gene plot
	fSet <- base::paste( prefix, gTarget, suffix, extension, sep="")

	# don't let bad OS filename characters ruin the links
	fSet <- file.cleanSpecialCharactersFromFileName( fSet)

	urlSet <- file.path( path, fSet)
	if ( type == "excel") {
		urlSet <- gsub( "/", "\\", urlSet, fixed=T)
		txt <- base::paste("=HYPERLINK(\"", urlSet, "\",\"", gText, "\")", sep="")
	} else {
		txt <- as.link( urlSet, gText)
	}
	return( txt)
}


writeAsWebPage <- function( df, file, title="FactDB Results:  ", highlightRows=NULL, 
			extraHTMLtextFile=NULL, extraHTMLtext=NULL) {

	# catch standard out into that file
	sink( file=file)
	on.exit( sink())
	HTMLheader( title=file)
	# stuff the file name in as a heading line...
	for ( text in title) {
		tag( "h3"); cat( text); untag( "h3")
	}
	tag( "i")
	cat( "Updated: ", date())
	untag( "i")
	br(1)
	HTMLincludeText( extraHTMLtext)
	HTMLincludeFile( extraHTMLtextFile)
	if ( nrow( df) > 0) {
		HTMLtable( df, rowNumbers=TRUE, cellpadding=3, highlightRows=highlightRows)
	} else {
		br(1)
		tag( "h4"); cat( "Empty Table with zero rows..."); untag( "h4")
	}
	HTMLfooter()
}

#
# CGIwithRtools.R  - modified from the "CGIwithR" package
# Bob Morrison   4/08
#	7/29	- added a "table" function

# needed by some of the CGIwithRtools.R routines, so set it by hand now...
options("useUnquotedTagnames"=TRUE)


CGIparse <- function(string=Sys.getenv("QUERY_STRING"), collapse = TRUE)
{
    boundary <- getMultiPartBoundary()  
    if(length(boundary)) 
       return(parseMultiPartFormData(boundary, string, splitLines = TRUE))
      
    the.split.string <- lapply(strsplit(string, "&"), 
                             function(string)strsplit(string, "="))[[1]]
    arglist <- lapply(the.split.string, function(pair) pair[2])
    names(arglist) <- sapply(the.split.string, function(pair) pair[1])
    ans <- lapply(arglist, hexDecode)

    if(collapse) {
      w <- which(duplicated(names(ans)))
      for(i in w) {
         first <- base::match(names(ans)[i], names(ans))
         ans[first] <- c(ans[first], ans[i])
      }
    }

    ans
}
      

hexDecode <- function(string){
    string <- base::gsub("\\+", " ", string)
    string <- base::gsub("%09", "\t", string)
    string <- base::gsub("%0D%0A", "\n", string)
    string <- base::gsub("%0D|%0A", "\n", string)
    pieces <- strsplit(string, "%")[[1]]
    dehex <- function(string){
              hex <- base::substr(string, 1, 2)
              base::paste(ascii[hex], substring(string, 3), sep = "")
              }
    pieces <- c(pieces[1], sapply(pieces[-1], dehex))
    base::paste(pieces, collapse = "")
}


HTMLheader <-
function(title = character(0), css = character(0), hasBody=TRUE)
{
    cat("<!doctype html public \"-//W3C/DTD HTML 4.0/EN\">\n")
    lf(2)

    cat("<HTML>\n<HEAD>\n")
    if(length(title))
        cat("<title>", base::paste(title, collapse = " "), "</title>\n", sep = "")

    if(length(css)) {
        cssNames = names(css)
        base::sapply(seq(along = css),
               function(i) {
                       cat('<LINK rel="StyleSheet" HREF="', css[i], '"',
                           ifelse(length(cssNames) && cssNames[i] != "", base::paste(' MEDIA="', cssNames[i], '" ', sep = ""), ""),
                           'TYPE="text/css">\n', sep = "")
               })
    }

    
    cat("</HEAD>\n")
    if (hasBody) cat("\n<BODY>\n")
}


HTMLfooter =
  #
  # Output the closing </BODY></HTML>
  #
function( hasBody=TRUE)
{
  if (hasBody) cat("\n</BODY>")
  cat("\n</HTML>\n")
}


HTMLincludeText <- function( text) {
	if (is.null( text)) return()
	if (length( text) > 0) {
		cat( base::paste( text, collapse="\n"))
		br(2)
	}
}


HTMLincludeFile <- function( file) {
	if (is.null( file)) return()
	if ( file.exists( file)) {
		txt <- readLines( con=file)
		if (length( txt) > 0) {
			cat( base::paste( txt, collapse="\n"))
			br(2)
		}
	}
}


HTMLtable <- function( df, rowNumbers=FALSE, cellpadding=4, border=1, highlightRows=NULL) {
	# if row numbers wanted, do it explicitly
	if (rowNumbers) {
		dfNames <- colnames( df)
		dfRownames <- rownames(df)
		if ( is.null(dfRownames)) dfRownames <- as.character(1:nrow(df))
		if ( ! all( as.character( as.integer( dfRownames)) == dfRownames)) {
			dfRownames <- 1:nrow(df)
		}
		df <- data.frame( dfRownames, df, stringsAsFactors=FALSE)
		colnames(df) <- c("Row", dfNames)
	}
	cat( "<table cellpadding=", cellpadding, " border=", border, "> \n", sep="")
	# bold the column headings, turning any 'spaces' that got turned to preiods back to spaces...
	HTMLtableRow( base::sub( ".", " ", colnames(df), fixed=TRUE), bold=TRUE, underline=FALSE)
	for ( i in 1:nrow(df)) {
		oneRow <- as.character( df[ i, ])
		oneRow[ oneRow %in% c("", " ", "\t")] <- " &nbsp; "
		highlight <- NULL
		if ( !is.null(highlightRows) && !is.na( highlightRows[i])) highlight <- highlightRows[i]
		HTMLtableRow( oneRow, bgcolor=highlight)
	}
	cat( "</table>", "\n")
}


HTMLtableRow <- function( text, bold=FALSE, underline=FALSE, bgcolor=NULL) {
	if (underline) text <- base::paste( "<u>", text, "</u>")
	if (bold) text <- base::paste( "<b>", text, "</b>")
	tr <- "<tr>"
	if ( !is.null( bgcolor) && !is.na( bgcolor)) tr <- base::paste( "<tr  bgcolor=\"", bgcolor, "\">", sep="")
	cat( tr, base::paste( "<td>", text, "</td>", collapse=" "), "</tr>", "\n")
}


as.link <- function( url, label, target=NULL, style="") {
	if ( is.null( target)) {
		result <- base::paste( "<a ", style, " href=\"", url, "\"> ", label, " </a>", sep="")
	} else {
		result <- base::paste( "<a ", style, " href=\"", url, "\"  target=\"", target, "\"> ", label, " </a>", sep="")
	}
	return( result)
}


br <- function(n = 1){
    cat(base::paste(rep("<BR>", n), collapse = ""), "\n")
    }


tag <- function(tagname, ...) 
{
    if (getOption("useUnquotedTagnames")) {
        result <- as.character(substitute(tagname))
    }
    dots <- list(...)
    if (length(dots) > 0) {
        dotnames <- names(dots)
        dots <- base::paste(dotnames, base::paste(dots, "\"", sep = ""),
                      sep = "=\"")
        dots <- base::paste(dots, collapse = " ")
        result <- base::paste(result, dots, sep = " ")
    }
    cat(base::paste("<", result, ">", sep = ""))
}
  

untag <- function(tagname){
    if (getOption("useUnquotedTagnames")) 
        tagname <- as.character(substitute(tagname))
    cat("</", tagname, ">", sep = "")
}
  

lf <- function(n = 1) cat(base::paste(rep("\n", n), sep = ""))


comments <- function(text) cat("<!--", text, "-->")


mailto <- function(text, address){
    cat("<a href=\"mailto:", address, "\">", text, "</a>", sep="")
}
    

linkto <- function(text, URL){
    cat("<A href=\"", URL, "\">", text, "</A>", sep = "")
}

 
webPNG <- function(file, ..., graphDir){
    if (missing(graphDir)) {
      if(!exists("graphDir", envir = globalenv(), inherits = TRUE)) {
        cat("Error in webPNG(): graphDir not set\n\n")
        q("no") ##  abort if graphDir not specified
      }
      graphDir = get("graphDir", envir = globalenv(), inherits = TRUE)
    }

    if (!file.exists(graphDir)) {
        cat("Error in webPNG():", graphDir, "does not exist\n\n")
        q("no") ##  abort if specified graphDir does not exist
    }
    n <- base::nchar(graphDir)
    separator <- if ( base::substr(graphDir, n, n) == "/") "" else "/"

    if(require(GDD, quietly = TRUE)) {
       GDD(file = base::paste(graphDir, file, sep = separator), ...)
    } else
       bitmap(file = base::paste(graphDir, file, sep = separator), ...)
    invisible(NULL)
}

    
img <- function (src, ..., graphURLroot = "") 
{
    result <- src
    if(missing(graphURLroot)) {
       if (exists("graphURLroot", envir = globalenv(), inherits = TRUE))
         graphURLroot = get("graphURLroot", envir = globalenv(), inherits = TRUE)
       else
         graphURLroot = ""
    }
    result <- file.path(graphURLroot, src)
    result <- base::paste("<IMG SRC=\"", result, "\"", sep = "")
    dots <- list(...)
    if (length(dots) > 0) {
        dotnames <- names(dots)
        dots <- base::paste(dotnames, base::paste(dots, "\"", sep = ""),
                      sep = "=\"")
        dots <- base::paste(dots, collapse = " ")
        result <- base::paste(result, " ", dots, ">", sep = "")
    } else result <- base::paste(result, ">", sep = "")
    cat(result)
    invisible(result)
}

"ascii" <-
  structure(c(
  "\x01", "\x02", "\x03", "\x04", "\x05", "\x06", "\x07",
  "\x08", "\x09", "\x0A", "\x0B", "\x0C", "\x0D", "\x0E", "\x0F",
  "\x10", "\x11", "\x12", "\x13", "\x14", "\x15", "\x16", "\x17",
  "\x18", "\x19", "\x1A", "\x1B", "\x1C", "\x1D", "\x1E", "\x1F",
  "\x20", "\x21", "\x22", "\x23", "\x24", "\x25", "\x26", "\x27",
  "\x28", "\x29", "\x2A", "\x2B", "\x2C", "\x2D", "\x2E", "\x2F",
  "\x30", "\x31", "\x32", "\x33", "\x34", "\x35", "\x36", "\x37",
  "\x38", "\x39", "\x3A", "\x3B", "\x3C", "\x3D", "\x3E", "\x3F",
  "\x40", "\x41", "\x42", "\x43", "\x44", "\x45", "\x46", "\x47",
  "\x48", "\x49", "\x4A", "\x4B", "\x4C", "\x4D", "\x4E", "\x4F",
  "\x50", "\x51", "\x52", "\x53", "\x54", "\x55", "\x56", "\x57",
  "\x58", "\x59", "\x5A", "\x5B", "\x5C", "\x5D", "\x5E", "\x5F",
  "\x60", "\x61", "\x62", "\x63", "\x64", "\x65", "\x66", "\x67",
  "\x68", "\x69", "\x6A", "\x6B", "\x6C", "\x6D", "\x6E", "\x6F",
  "\x70", "\x71", "\x72", "\x73", "\x74", "\x75", "\x76", "\x77",
  "\x78", "\x79", "\x7A", "\x7B", "\x7C", "\x7D", "\x7E", "\x7F",
  "\x80", "\x81", "\x82", "\x83", "\x84", "\x85", "\x86", "\x87",
  "\x88", "\x89", "\x8A", "\x8B", "\x8C", "\x8D", "\x8E", "\x8F",
  "\x90", "\x91", "\x92", "\x93", "\x94", "\x95", "\x96", "\x97",
  "\x98", "\x99", "\x9A", "\x9B", "\x9C", "\x9D", "\x9E", "\x9F",
  "\xA0", "\xA1", "\xA2", "\xA3", "\xA4", "\xA5", "\xA6", "\xA7",
  "\xA8", "\xA9", "\xAA", "\xAB", "\xAC", "\xAD", "\xAE", "\xAF",
  "\xB0", "\xB1", "\xB2", "\xB3", "\xB4", "\xB5", "\xB6", "\xB7",
  "\xB8", "\xB9", "\xBA", "\xBB", "\xBC", "\xBD", "\xBE", "\xBF",
  "\xC0", "\xC1", "\xC2", "\xC3", "\xC4", "\xC5", "\xC6", "\xC7",
  "\xC8", "\xC9", "\xCA", "\xCB", "\xCC", "\xCD", "\xCE", "\xCF",
  "\xD0", "\xD1", "\xD2", "\xD3", "\xD4", "\xD5", "\xD6", "\xD7",
  "\xD8", "\xD9", "\xDA", "\xDB", "\xDC", "\xDD", "\xDE", "\xDF",
  "\xE0", "\xE1", "\xE2", "\xE3", "\xE4", "\xE5", "\xE6", "\xE7",
  "\xE8", "\xE9", "\xEA", "\xEB", "\xEC", "\xED", "\xEE", "\xEF",
  "\xF0", "\xF1", "\xF2", "\xF3", "\xF4", "\xF5", "\xF6", "\xF7",
  "\xF8", "\xF9", "\xFA", "\xFB", "\xFC", "\xFD", "\xFE", "\xFF"
),
  .Names = c(
  "01", "02", "03", "04", "05", "06", "07",
  "08", "09", "0A", "0B", "0C", "0D", "0E", "0F",
  "10", "11", "12", "13", "14", "15", "16", "17",
  "18", "19", "1A", "1B", "1C", "1D", "1E", "1F",
  "20", "21", "22", "23", "24", "25", "26", "27",
  "28", "29", "2A", "2B", "2C", "2D", "2E", "2F",
  "30", "31", "32", "33", "34", "35", "36", "37",
  "38", "39", "3A", "3B", "3C", "3D", "3E", "3F",
  "40", "41", "42", "43", "44", "45", "46", "47",
  "48", "49", "4A", "4B", "4C", "4D", "4E", "4F",
  "50", "51", "52", "53", "54", "55", "56", "57",
  "58", "59", "5A", "5B", "5C", "5D", "5E", "5F",
  "60", "61", "62", "63", "64", "65", "66", "67",
  "68", "69", "6A", "6B", "6C", "6D", "6E", "6F",
  "70", "71", "72", "73", "74", "75", "76", "77",
  "78", "79", "7A", "7B", "7C", "7D", "7E", "7F",
  "80", "81", "82", "83", "84", "85", "86", "87",
  "88", "89", "8A", "8B", "8C", "8D", "8E", "8F",
  "90", "91", "92", "93", "94", "95", "96", "97",
  "98", "99", "9A", "9B", "9C", "9D", "9E", "9F",
  "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7",
  "A8", "A9", "AA", "AB", "AC", "AD", "AE", "AF",
  "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7",
  "B8", "B9", "BA", "BB", "BC", "BD", "BE", "BF",
  "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
  "C8", "C9", "CA", "CB", "CC", "CD", "CE", "CF",
  "D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7",
  "D8", "D9", "DA", "DB", "DC", "DD", "DE", "DF",
  "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7",
  "E8", "E9", "EA", "EB", "EC", "ED", "EE", "EF",
  "F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7",
  "F8", "F9", "FA", "FB", "FC", "FD", "FE", "FF"
  )
  )


scanText <- function(string, what = character(0), ...){
    tc <- textConnection(string)
    result <- scan(tc, what = what, quiet = TRUE, ...)
    close(tc)
    return(result)}
    

indentPrint <- function(object, indent = 4, ...){
    tc <- textConnection("zz", "w", local = TRUE)
    sink(tc)
    try(print(object, ...))
    sink()
    close(tc)
    indent <- base::paste(rep(" ", indent), sep = "", collapse = "")
    cat(base::paste(indent, zz, sep = ""), sep = "\n")}
    

getParameters =
  #
  #
  #
function(str = Sys.getenv("QUERY_STRING"))
{
  params = strsplit(strsplit(str, "&")[[1]], "=")
  print( base::paste("\n", params, "\n"))
  names(params) = sapply(params, function(x) x[1])
  params = lapply(params, function(x) x[-1])
}

writeRequestInfo =
function(env = c("PATH_INFO", "HTTP_HOST", "HTTP_USER_AGENT", "HTTP_REFERER", "REMOTE_ADDR"), size = -1)
{
   cat('<center><hr width="50%"></center>\n')
   cat('<font size="', as.integer(size), '">\n', sep = '')
   cat("Generated by CGIwithR, version 1.0 (SBRI mod.) <br/>\n")
   cat("QUERY_STRING = <pre>", Sys.getenv("QUERY_STRING"), "</pre><p/>\n")

   if(is.logical(env) && env)
     env = Sys.getenv()

   if(length(env))
     showEnvironmentVariables(env)

   cat('</font>\n')

   invisible(env)
}


showEnvironmentVariables =
  #
  # Display the values of some or all environment variables
  # for help in debugging or reporting the settings for the
  # actual computations in the CGI script.
  #
function(env = Sys.getenv())
{
 if(is.null(names(env))) 
    env = Sys.getenv(env)

  cat("<table>\n")
  for(i in names(env))
      cat('<tr align="left"><th>', i, '</th><th>', env[i], '</th></tr>\n')
  cat("</table>\n")

  invisible(env)
}






splitDataString =
  #  Split the given line by sep[1] and then each of those elements by sep[2].
  #  and use the names of the name = value pairs generated from sep[2]
  #  as the names for the resulting list.
  #  e.g.  'form-data; name="dataFile"; filename="bootData"'
  # 
function(txt, sep = c("; ", "="))
{
  els = strsplit(txt, sep[1])[[1]]

  if(length(sep) > 1) {
    els = sapply(els, strsplit, sep[2])
    tmp = sapply(els, function(x) x[1])
    els = lapply(els, function(x) if(length(x) == 1) x[1] else x[-1])
    names(els) = tmp
  }

  els
}

getMultiPartBoundary =
  #
  # Determines whether the request is a multipart/form-data
  # and if so, extracts the boundary that identifies the different
  # segments of the request.
  # This returns either the boundary or a zero-length character vector.
  # If an actual string is returned, then the contents of the request
  # were retrieved from the standard input. This is important if
  # we use R directly rather than via cgi.R.
function()
{  
   if(Sys.getenv("REQUEST_METHOD") == "POST") {
       type = Sys.getenv("CONTENT_TYPE")
       if(type != "") {
          els = strsplit(type, "; ")[[1]]
          idx = base::match("multipart/form-data", els)
          if(!is.na(idx)) {
             els = sapply(els[-idx], strsplit, "=")
             tmp = sapply(els, function(x) x[1])
             els = sapply(els, function(x) x[-1])
             names(els) = tmp

             if("boundary" %in% names(els))
               return(els["boundary"])
          }
       }
    }

    return(character(0))
}


parseMultiPartFormData =
  #
  # This is used when we have a form with enctype (encoding type)
  # of multipart/form-data.
  # This is quite similarl to attachments in a mail message.
  # Each "element" is identified by a pair of --boundary lines
  # (and the last one with a --boundary--)
  # Each element is of the formi.e. 
  #   field: text
  #   field: text
  #   <blank line>
  #   body....

  #
  #  This returns an object of class "MultiPartFormData"
  #  which is a list of the processed elements indexed by their
  #  name.   Each element typically has a Content-Disposition field
  #  and within this there is a name="text" string. We use that for the
  #  name of the element in the R.
  #
  #
  # XXX
function(boundary, content = readLines(stdin()), collapseBody = "\n",
          splitLines = FALSE)
{
    # If the content was not specified, read it from the standard input.
    # This would be used if we were calling R directly rather than indirectly
    # via R.cgi.  This may happen soon.
    # If the cotent is specified ans so is sep, then we need to break
    # the string into a vector of lines. So we split first by new line
    # and then by control-feed.  The reason for doing this in two steps
    # is that the contents of a file may be embedded within the input and
    # separated by \n rather than \r\n. So not all lines are identified by the
    # \r\n, only those created by the form mechanism.
  if(!missing(content) && splitLines) {
       # Read the lines from standard input
     content = strsplit(content, "\n")[[1]]
     content = base::gsub("\r$", "", content)    
  }

    # Find the locations of the boundary lines, i.e. --boundary
    # Won't match the last one since that is --boundary--, but we don't need that one.
  starts = grep(base::paste("^--", boundary, "$", sep = ""), content)

    # Loop over each of the positions of the start of a block in the form
    # and process the lines within that block.
 vals = lapply(seq(along = starts),
         function(i) {
              # Get the lines of the block.
            txt = content[seq(starts[i], ifelse(!is.na(starts[i+1]), starts[i+1], length(content)) - 1)]

              # Find the first blank line which separates the header from the body of the block.
            b = base::match("", txt)
              #
            body = txt[-(1:(b))]
            if(!is.na(collapseBody))
              body = base::paste(body, collapse = collapseBody)
              
            header = txt[2:(b-1)]
         
             # Break each line of the header by the :,
             # e.g. Content-Type: value...
            els = strsplit(header, ": ")
            fieldNames = sapply(els, function(x) x[1])
            fields = lapply(els, function(x) {
                                    splitDataString(base::paste(x[-1], collapse = ": "))
                                  })
         
            
            names(fields) = fieldNames
            list(fields = fields, body = body)
  })

    # Now put the names of the actual data elements of the form onto the R form
    # object itself.
  names(vals)  = stripQuotes(sapply(vals, function(x) x$fields[["Content-Disposition"]]["name"]))

   # Handle duplicated name elements by combining them into a vector for that name.
  w = duplicated(names(vals))
  if(any(w)) {
    sapply(which(w),
            function(i) {
               first = base::match(names(vals)[i], names(vals))
               vals[[first]][["body"]] <<- c(vals[[first]]$body, vals[[i]]$body)
             })
    vals = vals[!w]
  }

  class(vals) <- "MultiPartFormData"
  vals
}

"$.MultiPartFormData" <-
  # Method to access the body associated with the given Content-Disposition: name="name" value
  # in the form.
function(x, name)
{
  el = x[[name]]
  el$body
}



stripQuotes =
  #
  # Remove leading and trailing quotes.
  # This is used in dealing with the name="value" pairs in the
  # multipart/form-data fields, e.g. the Content-Disposition: field.
function(txt) {
  base::gsub('^"(.*)"$', "\\1", txt)
}
