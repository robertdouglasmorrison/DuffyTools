# convertHypertext.R  - undo special hypertext substitutions

convertHypertext <- function( txt) {

	hexTable <- c( "%20","%21","%22","%23","%24","%25","%26","%27",
			"%28","%29","%2A","%2B","%2C","%2D","%2E","%2F",
			            "%2a","%2b","%2c","%2d","%2e","%2f",
		       "%3A","%3B","%3C","%3D","%3E","%3F",
		       "%3a","%3b","%3c","%3d","%3e","%3f",
		       "%40","%60",
		       "%5B","%5C","%5D","%5E","%5F",
		       "%5b","%5c","%5d","%5e","%5f",
		       "%7B","%7C","%7D","%7E",
		       "%7b","%7c","%7d","%7e")
	valTable <- c( " ",  "!",  "\"", "#",  "$",  "%",  "&",  "'",  
			"(",  ")",  "*",  "+",  ",",  "-",  ".",  "/", 
			            "*",  "+",  ",",  "-",  ".",  "/",
	               ":",  ";",  "<",  "=",  ">",  "?", 
	               ":",  ";",  "<",  "=",  ">",  "?", 
		       "@",  "`",
		       "[",  "\\", "]",  "^",  "_",  
		       "[",  "\\", "]",  "^",  "_",
		       "{",  "|",  "}",  "~",
		       "{",  "|",  "}",  "~")

	textTable <- c( "&apos;", "&colon;", "&quot;", "&nbsp;")
	textValTable <- c( "'", ":", "'", " ")

	newtxt <- txt

	for ( j in 1:length(txt)) {
	    repeat {
		marks <- gregexpr( "%", newtxt[j], fixed=TRUE)
		if ( marks[[1]][1] < 1) break
		gothit <- FALSE
		for ( from in marks[[1]]) {
			esc <- substr( newtxt[j], from, from+2)
			hit <- match( esc, hexTable, nomatch=0)[1]
			if ( hit > 0) {
				newtxt[j] <- sub( hexTable[hit], valTable[hit], newtxt[j], fixed=TRUE)
				gothit <- TRUE
				break
			}
		}
		# if we made no changes, we're done
		if ( ! gothit) break
	    }
	}

	for ( j in 1:length(txt)) {
	    repeat {
		marks <- gregexpr( "&", newtxt[j], fixed=TRUE)
		if ( marks[[1]][1] < 1) break
		endmarks <- gregexpr( ";", newtxt[j], fixed=TRUE)
		gothit <- FALSE
		for ( k in 1:length(marks[[1]])) {
			from <- marks[[1]][k]
			to <- endmarks[[1]][k]
			amper <- substr( newtxt[j], from, to)
			hit <- match( amper, textTable, nomatch=0)[1]
			if ( hit > 0) {
				newtxt[j] <- sub( textTable[hit], textValTable[hit], newtxt[j], fixed=TRUE)
				gothit <- TRUE
				break
			}
		}
		# if we made no changes, we're done
		if ( ! gothit) break
	    }
	}

	return( newtxt)
}
