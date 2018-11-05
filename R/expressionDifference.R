# expressionDifference.R

# run a matrix of expression intensities  through a differential expression test

`expressionDifference` <- function( x, classes, offset=1, FUN=wilcox.test, ...) {

	# decipher the classes arg to know which columns to use
	if ( length( classes) != ncol(x)) stop( "'classes' arg must equal 'ncol(x)'")

	lvls <- levels( factor( classes))
	NL <- length(lvls)

	NG <- nrow(x)
	gnames <- rownames(x)
	outF <- outP <- rep.int(NA, NG)
	out <- vector( mode="list")
	nout <- 0

	for ( i in 1:NL) {
		lab1 <- lvls[i]
		lab2 <- if (NL == 2) lvls[3-i] else paste( "Not", lvls[i], sep="_")
		myLabel <- paste( lab1, lab2,sep=".v.")
		cat( "\n  ", myLabel)
		cols1 <- which( classes == lab1)
		cols2 <- setdiff( which( !is.na(classes)), cols1)

		for ( j in 1:NG) {
			v1 <- x[ j, cols1] + offset
			v2 <- x[ j, cols2] + offset
			ans <- FUN( v1, v2, ...)
			outP[j] <- ans$p.value
			outF[j] <- log2( mean(v1, na.rm=T) / mean( v2, na.rm=T))
		}
		sml <- data.frame( "GENE_ID"=gnames, "LOG2_FOLD"=outF, "P_VALUE"=outP, stringsAsFactors=F)
		ord <- diffExpressRankOrder( sml$LOG2_FOLD, sml$P_VALUE)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:nrow(sml)

		nout <- nout + 1
		out[[nout]] <- sml
		names(out)[nout] <- myLabel
	}

	return( out)
}


