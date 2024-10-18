# lifeCycleTools.R

`verifyLifeCycleSetup` <- function() {

	isReady <- exists( "VectorSpace", envir=LifeCycleEnv)
	isRightSpecies <- FALSE
	if ( isReady) {
		curSpecies <- get( "Species", envir=LifeCycleEnv)
		if( ! is.null( curSpecies)) isRightSpecies <- ( curSpecies == getCurrentSpecies())
	}
	if ( !isReady || !isRightSpecies) LifeCycleSetup( dataset="RNAseq", unitVectorMode="absolute")
	return()
}


`LifeCycleSetup` <- function( dataset=c( "RNAseq", "YoungWinzeler", "DeRisi", "Custom"), 
				unitVectorMode=c("absolute","relative","none"), 
				stageMode=c("all", "blood.only"),
				custom.file=NULL, preNormalize=TRUE, postNormalize=TRUE, 
				doRMA=FALSE, min.spread=3.0, verbose=FALSE) {

	dataset <- match.arg( dataset)
	unitVectorMode <- match.arg( unitVectorMode)
	stageMode <- match.arg( stageMode)
	ans <- NULL

	cat( "\nSetting up LifeCycle dataset:  ", dataset)

	# get the data from the package, its called 'lifeCycleData'
	if ( dataset == "RNAseq") {
		data( LifeCycleData_RNAseq, package="DuffyTools", envir=environment())

		# use explicit colnames lookup
		sgspzPtrs <- match( c("sgSPZ.3D7","sgSPZ.NF54", "Spz1","Spz2","Spz3","Spz4", 
					"SGspz1","SGspz2","SGspz3","SPZ.SG"), colnames(lifeCycleData))
		lsPtrs <- match( c("Day4.GFP.r1.TR3","Day4.GFP.r2","Day5.GFP.r1","Day5.GFP.r2","Day5.GFP.r4", 
					"Day6.GFP.r1","Day6.GFP.r2","Day6.GFP.r4"), colnames(lifeCycleData))
		merozPtrs <- match( c("OttoHr48","StunHr5","HoeijHr05"), colnames(lifeCycleData))
		eringPtrs <- match( c("OttoHr0","OttoHr8","StunHr10","HoeijHr10"), colnames(lifeCycleData))
		lringPtrs <- match( c("OttoHr16","StunHr15","HoeijHr15","StunHr20","HoeijHr20",
					"Ring1","Ring2"), colnames(lifeCycleData))
		etrophPtrs <- match( c("OttoHr24","StunHr25","HoeijHr25",
					"HoeijHr30","StunHr30"), colnames(lifeCycleData))
		ltrophPtrs <- match( c("OttoHr32","StunHr35","HoeijHr35"), colnames(lifeCycleData))
		schizPtrs <- match( c("OttoHr40","StunHr40","HoeijHr40","SanariaBLOOD_1",
					"SanariaBLOOD_2"), colnames(lifeCycleData))
		gcyte5Ptrs <- match( c("Su.GcyteV","AO_Gcyte","SauerGcyte_Female","SauerGcyte_Male1",
					"SauerGcyte_Male2","Gcyte1","Gcyte2"), colnames(lifeCycleData))
		ookinePtrs <- match( c("Su.Ookinete3D7","Su.Ookinete7G8"), colnames(lifeCycleData))
		oospzPtrs <- match( c("ooSPZ.NF54","OOspz1","OOspz2","OOspz3"), colnames(lifeCycleData))
		allPtrs <- c( sgspzPtrs, lsPtrs, merozPtrs, eringPtrs, lringPtrs, etrophPtrs, ltrophPtrs,
				schizPtrs, gcyte5Ptrs, ookinePtrs, oospzPtrs)
		if ( any ( is.na( allPtrs))) {
			bad <- which( is.na( allPtrs))
			cat( "\nBad/missing Life Cycle Stage column name: ", bad)
			stop()
		}

		# these are column numbers in the original file
		ans <- buildLifeCycleVectorSpace( file=NULL, tbl=lifeCycleData, geneColumn="Gene",
			dimensions=list( "SalivarySPZ"=sgspzPtrs, "LiverStage"=lsPtrs, "Merozoite"=merozPtrs,
			"EarlyRing"=eringPtrs, "LateRing"=lringPtrs, "EarlyTroph"=etrophPtrs, "LateTroph"=ltrophPtrs, 
			"Schizont"=schizPtrs, "Gametocyte"=gcyte5Ptrs, "Ookinete"=ookinePtrs, 
			"OocystSPZ"=oospzPtrs), 
			unitVectorMode=unitVectorMode, min.spread=min.spread, preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, stageMode=stageMode, verbose=verbose)
		ans2 <- buildLifeCycleVectorSpace( file=NULL, tbl=lifeCycleData, geneColumn="Gene",
			dimensions=list( "SalivarySPZ"=sgspzPtrs, "LiverStage"=lsPtrs, "Merozoite"=merozPtrs,
			"EarlyRing"=eringPtrs, "LateRing"=lringPtrs, "EarlyTroph"=etrophPtrs, "LateTroph"=ltrophPtrs, 
			"Schizont"=schizPtrs, "Gametocyte"=gcyte5Ptrs, "Ookinete"=ookinePtrs, 
			"OocystSPZ"=oospzPtrs), 
			unitVectorMode="none", min.spread=NULL, preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, stageMode=stageMode, verbose=FALSE)
	}

	# get the data from the package, its called 'tbl'
	if ( dataset == "YoungWinzeler") {
		data( LifeCycleData_YoungWinzeler, package="DuffyTools", envir=environment())

		# these are column numbers in the original file
		ans <- buildLifeCycleVectorSpace( file=NULL, tbl=tbl, geneColumn="Gene",
			dimensions=list( "Sporozoite"=5, "Merozoite"=c(12,19), "EarlyRing"=c(6,13), 
			"LateRing"=c(7,14), "EarlyTroph"=c(8,15), "LateTroph"=c(9,16), 
			"EarlySchiz"=c(10,17), "LateSchiz"=c(11,18), "EarlyGameto"=c( 24:26), 
			"LateGameto"=c( 27:29)), unitVectorMode=unitVectorMode,
			min.spread=min.spread, preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
		ans2 <- buildLifeCycleVectorSpace( file=NULL, tbl=tbl, geneColumn="Gene",
			dimensions=list( "Sporozoite"=5, "Merozoite"=c(12,19), "EarlyRing"=c(6,13), 
			"LateRing"=c(7,14), "EarlyTroph"=c(8,15), "LateTroph"=c(9,16), 
			"EarlySchiz"=c(10,17), "LateSchiz"=c(11,18), "EarlyGameto"=c( 24:26), 
			"LateGameto"=c( 27:29)), unitVectorMode="none",
			min.spread=NULL, preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=FALSE)
	}

	if ( dataset == "DeRisi") {
		data( LifeCycleData_DeRisi, package="DuffyTools", envir=environment())
		ans <- buildLifeCycleVectorSpace( file=NULL, tbl=tbl, geneColumn="Gene",
			dimensions=list( "EarlyRing"=1:6, "MidRing"=7:12, "LateRing"=13:18,
			"EarlyTroph"=19:25, "LateTroph"=26:31, "EarlySchiz"=32:39, "LateSchiz"=40:48), 
			unitVectorMode=unitVectorMode, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
		ans2 <- buildLifeCycleVectorSpace( file=NULL, tbl=tbl, geneColumn="Gene",
			dimensions=list( "EarlyRing"=1:6, "MidRing"=7:12, "LateRing"=13:18,
			"EarlyTroph"=19:25, "LateTroph"=26:31, "EarlySchiz"=32:39, "LateSchiz"=40:48), 
			unitVectorMode="none", min.spread=NULL, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=FALSE)
	}

	if ( dataset == "Custom") {

		# pre read the file to make the dimensions construct
		tmp <- read.delim(custom.file, as.is=T)
		dimNames <- setdiff( colnames(tmp), "Gene")
		dimensions <- vector( mode="list")
		for (k in 1:length(dimNames)) dimensions[[k]] <- match( dimNames[k], colnames(tmp))
		names(dimensions) <- dimNames

		cat( "\nBuilding custom Life Cycle dataset:   Dimensions:\n")
		print( dimensions)
		cat( "\nFirst genes: \n")
		print( head( tmp))


		ans <- buildLifeCycleVectorSpace( file=custom.file, tbl=NULL, geneColumn="Gene",
			dimensions=dimensions,
			unitVectorMode=unitVectorMode, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
		ans2 <- buildLifeCycleVectorSpace( file=custom.file, tbl=NULL, geneColumn="Gene",
			dimensions=dimensions,
			unitVectorMode="none", min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
	}


	LifeCycleEnv[[ "VectorSpace" ]] <- ans
	LifeCycleEnv[[ "IntensitySpace" ]] <- ans2
	LifeCycleEnv[[ "Species" ]] <- getCurrentSpecies()
	if ( ! is.null( ans)) {
		LifeCycleEnv[[ "N_STAGES" ]] <- ncol(ans) - 2   # geneID, Product come first
		LifeCycleEnv[[ "STAGE_NAMES" ]] <- colnames(ans)[3:ncol(ans)]
	}
	cat( "\nDone.\n")
}


`getLifeCycleMatrix` <- function( mode=c("IntensitySpace", "VectorSpace")) {

	verifyLifeCycleSetup()
	mode <- match.arg( mode)

	tbl <- LifeCycleEnv[[ mode]]
	genes <- tbl$GENE_ID

	# get just the values, not the IDs and Products
	m <- as.matrix( tbl[ ,3:ncol(tbl)])
	rownames(m) <- genes
	return(m)
}


`buildLifeCycleVectorSpace` <- function( file="", tbl=NULL, dimensions, geneColumn="Gene", 
					unitVectorMode=c( "absolute", "relative", "none"), 
					stageMode=c("all", "blood.only"),
					min.value=0.01, min.spread=2, preNormalize=FALSE, 
					postNormalize=TRUE, doRMA=TRUE, verbose=FALSE) {

	unitVectorMode <- match.arg( unitVectorMode)
	stageMode <- match.arg( stageMode)

	# get all the original reference set of gene data
	if ( is.null( tbl)) {
		tbl <- read.delim( file, as.is=TRUE)
	}

	# verify needed columns are all present
	if ( !( geneColumn %in% colnames(tbl))) {
		stop( paste( "gene name column not found:   tried: ", 
			geneColumn, "\tfound: ", colnames(tbl)))
	}

	allDs <- base::unlist( dimensions)
	isNumbers <-  ! any( is.na( colNumbers <- as.integer( allDs)))
	if (isNumbers) {
		if ( max( colNumbers) > ncol(tbl)) {
			stop( paste( "Invalid dimension column:   given: :", base::sort(colNumbers), 
				"\tColumns in table: ", ncol(tbl)))
		}
	} else {
		if ( !( all( allDs %in% colnames(tbl)))) {
			stop( paste( "some dimension name columns not found:   tried: ", 
				allD, "\tfound: ", colnames(tbl)))
		}
		allDs <- base::match( allDs, colnames(tbl))
		for( k in 1:length( dimensions)) dimensions[[k]] <- base::match( dimensions[[k]], colnames(tbl))
	}
	if ( verbose) {
		for( k in 1:length(dimensions)) cat( "\n\nStage: ", names(dimensions)[k], 
			"\n\tColumns: ", colnames(tbl)[ dimensions[[k]] ])
		cat( "\n")
	}

	# get the gene names and resolve to current annotation
	oldGenes <- tbl[[ geneColumn]]

	# allow a mapping by ortholog to the current species...  the Life cycle data is PF3D7 only
	thisSpecies <- getCurrentSpecies()
	if ( toupper(thisSpecies) %in% setdiff( toupper(PARASITE_SPECIES), "PF3D7")) {
		setCurrentSpecies( "Pf3D7")
		newGenes <- alias2Gene( oldGenes)
		orthos <- ortholog( newGenes, from="Pf3D7", to=thisSpecies)
		# drop rows with no ortholog
		keep <- which( orthos != "")
		if ( length(keep) < 1000) {
			cat( "\nProblem:  Not enough Pf3D7 GeneIDs detected.  Verify gene names...")
			stop()
		}
		tbl <- tbl[ keep, ]
		newGenes <- orthos[ keep]
		setCurrentSpecies( thisSpecies)
	} else {
		newGenes <- alias2Gene( oldGenes)
	}

	# extract just the wanted columns, impose a minimun intensity, and perhaps normalize the intensity columns
	dataColumns <- allDs
	mIn <- as.matrix( tbl[ , dataColumns])
	colnames( mIn) <- base::paste( "V", dataColumns, sep="")
	mIn[ is.na(mIn)] <- min.value
	mIn[ mIn < min.value] <- min.value

	# with microarrays, RMA is more appropriate, otherwise just plain quantile normalize
	if ( preNormalize) {
		if ( doRMA) {
			mIn <- duffyRMA( mIn, verbose=verbose)
		} else {
			mIn <- duffyMagnitudeNormalize( mIn, verbose=verbose)
		}
	}

	# there may be duplicate rows because of gene re-annotation or orthollogging...
	# combine those by mean average
	uGenes <- unique.default( newGenes)
	nUnique <- length( uGenes)
	mNew <- matrix( 0, nrow=nUnique, ncol=ncol(mIn))
	colnames(mNew) <- colnames(mIn)
	for ( i in 1:nUnique) {
		who <- which( newGenes == uGenes[i])
		if ( length(who) == 1) {
			mNew[ i, ] <- mIn[ who, ]
		} else {
			mNew[ i, ] <- apply( mIn[ who, ], MARGIN=2, FUN=sqrtmean)
		}
	}

	# use the dimensionality list to combine columns
	N_STAGES <- length( dimensions)
	STAGE_NAMES <- names( dimensions)
	mNewer <- matrix( nrow=nrow(mNew), ncol=N_STAGES)
	for( i in 1:N_STAGES) {
		who <- dimensions[[i]]
		where <- base::match( base::paste( "V", who, sep=""), colnames(mNew))
		if ( length(who) == 1) {
			oneV <- mNew[ , where]
		} else {
			oneV <- apply( mNew[ , where], MARGIN=1, FUN=sqrtmean)
		}
		mNewer[ , i] <- oneV
	}
	colnames( mNewer) <- STAGE_NAMES

	# once the stages have been set, now we can eliminate some if wanted...
	if ( stageMode == "blood.only") {
		keepers <- grep( "merozoite|ring|troph|schizont", tolower(colnames(mNewer)))
		# put a second copy of merozoite at the end, to complete the cycle
		mNewer <- mNewer[ , c(keepers,keepers[1])]
		colnames(mNewer)[ncol(mNewer)] <- "Merozoite2"
		N_STAGES <- ncol( mNewer)
		STAGE_NAMES <- colnames( mNewer)
	}

	# with microarrays, RMA is more appropriate, otherwise just plain quantile normalize
	if ( postNormalize) {
		if ( doRMA) {
			mNewer <- duffyRMA( mNewer, verbose=verbose)
		} else {
			mNewer <- duffyMagnitudeNormalize( mNewer, verbose=verbose)
		}
	}
	# now normalize each row to sum to unit vector (length=1)
	if ( unitVectorMode != "none") {
	    for( i in 1:nrow( mNewer)) {
		oneV <- mNewer[ i, ]
		# subtract the min so the proportions have at least one 0%
		if ( unitVectorMode == "relative") oneV <- oneV - min(oneV)
		oneV <- oneV / sum(oneV)
		mNewer[i, ] <- oneV
	    }
	}

	if ( ! is.null( min.spread)) {
		min.spread <- as.numeric( min.spread)
		mins <- apply( mNewer, MARGIN=1, min, na.rm=T)
		maxs <- apply( mNewer, MARGIN=1, max, na.rm=T)
		spreads <- maxs / mins
		drops <- which( spreads < min.spread)
		if ( length(drops) > 0) {
			mNewer <- mNewer[ -drops, ]
			uGenes <- uGenes[ -drops]
			cat( "\nDropping Genes with too little intensity change between dimensions: ", length(drops))
		}
	}

	ans <- data.frame( "GENE_ID"=uGenes, "PRODUCT"=gene2Product( uGenes), mNewer,
			stringsAsFactors=FALSE)
	rownames(ans) <- 1:nrow(ans)
	return( ans)
}
 

`calcLifeCycleStage` <- function( g, inten, dropGenes=vector()) {

	verifyLifeCycleSetup()

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- LifeCycleEnv[[ "VectorSpace"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	# build up a tally of how much intensity goes to each stage
	# do it for all genes, and the make stage histogram from the marker genes
	genes <- alias2Gene(g)

	# allow omitting some GeneIDs...
	if ( length( dropGenes)) {
		droppers <- alias2Gene( dropGenes)
		toDrop <- which( genes %in% droppers)
		if ( length( toDrop)) {
			g <- g[ -toDrop]
			genes <- genes[ -toDrop]
			inten <- inten[ -toDrop]
		}
	}

	# the Vector Space table is now in the current species, no need to ortholog here!
	where <- base::match( genes, vectorSpace$GENE_ID, nomatch=0)

	# build the matrix of stage fractions for every gene we were given.
	# switching to only use genes that are defined... others will not contribute
	myFractions <- matrix( 0, nrow=length(genes), ncol=N_STAGES)
	for( i in 1:N_STAGES) {
		myFractions[ where > 0, i] <- unitVectors[ where, i]
	}

	# account for background intensity to put microarray and RNA, etc., on equal footing .. now done elsewhere.!
	useInten <- inten

	# spread these intensitys over those stages
	myVectors <- myIntens <- matrix( 0, nrow=length(genes), ncol=N_STAGES)
	rownames(myVectors) <- rownames(myIntens) <- g
	colnames(myVectors) <- colnames(myIntens) <- STAGE_NAMES
	for( i in 1:N_STAGES) {
		myVectors[ , i] <- useInten * myFractions[ , i]
		myIntens[ , i] <- inten * myFractions[ , i]
	}

	# now do a summation by stage, and express as percentages...
	allSums <- apply( myVectors, MARGIN=2, FUN=sum, na.rm=T)
	bigSum <- sum( allSums)
	ans <- allSums * 100 / bigSum

	out <- list( "Stage"=ans, "IntensityVectors"=myIntens)
	return( out)
}


`calcLifeCycleStageFromFile` <- function( f, geneColumn="GENE_ID", intensityColumn="INTENSITY", 
					dropGenes=vector(), sep="\t") {

	verifyLifeCycleSetup()

	# open that file and find the needed columns
	tmp <- read.delim( f, as.is=T, sep=sep)
	if ( nrow( tmp) < 1) return(NULL)

	gset <- tmp[[ geneColumn]]
	if ( is.null(gset)) {
		cat( "calcLifeCycleStage:  gene name column not found.\nFile: ",f,
			"\nTried: ", geneColumn)
		return( NULL)
	}

	inten <- tmp[[ intensityColumn]]
	if ( is.null(inten)) {
		cat( "calcLifeCycleStage:  gene intensity column not found.\nFile: ",f,
			"\nTried: ", intensityColumn)
		return( NULL)
	}

	return( calcLifeCycleStage( gset, inten, dropGenes=dropGenes))
}


`plotLifeCycleStageFromFileSet` <- function( fnames, fids, fcolors=NULL, geneColumn="GENE_ID", 
		intensityColumn="RPKM_M", yMax=NULL, barSplines=FALSE, bar.lwd=5, legend.cex=1, max.labels=20,
		dropGenes=vector(), label="your label goes here...", sep="\t") {

	verifyLifeCycleSetup()

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nLifeStage:  Some transcript files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}

	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	# build the storage
	nFiles <- length( fnames)
	m <- matrix( nrow=nFiles, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# load each file in turn
	cat( "\nLoading:  ")
	for( i in 1:nFiles) {
		cat( " ", basename(fnames[i]))
		ans <- calcLifeCycleStageFromFile( fnames[i], geneColumn=geneColumn,
				intensityColumn=intensityColumn, dropGenes=dropGenes, sep=sep)
		m[ i, ] <- ans$Stage
	}
	cat( "\n")

	# plot it
	splineAns <- plotLifeCycleStages(m, col=fcolors, label=label, yMax=yMax, barSplines=barSplines,
			bar.lwd=bar.lwd, legend.cex=legend.cex, max.labels=max.labels)

	if ( is.null( splineAns)) {
		return( m)
	} else {
		return( list( "stage.vectors"=m, "spline.maximums"=splineAns))
	}
}


`plotLifeCycleStageFromMatrix` <- function( geneSet, intenMatrix, fids=colnames(intenMatrix), 
			fcolors=NULL, yMax=NULL, barSplines=FALSE, bar.lwd=5, legend.cex=1, max.labels=20,
			dropGenes=vector(), label="your label goes here...") {

	verifyLifeCycleSetup()

	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	# build the storage
	nColumns <- ncol( intenMatrix)
	m <- matrix( nrow=nColumns, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# load each dataset in turn
	for( i in 1:nColumns) {
		ans <- calcLifeCycleStage( geneSet, intenMatrix[ ,i], dropGenes=dropGenes)
		m[ i, ] <- ans$Stage
	}

	# plot it
	splineAns <- plotLifeCycleStages(m, col=fcolors, label=label, yMax=yMax, barSplines=barSplines,
			bar.lwd=bar.lwd, legend.cex=legend.cex, max.labels=max.labels)

	if ( is.null( splineAns)) {
		return( m)
	} else {
		return( list( "stage.vectors"=m, "spline.maximums"=splineAns))
	}
}


`plotLifeCycleStages` <- function( m, col=NULL, yMax=NULL, label="", barSplines=FALSE, 
				legend.cex=1.0, max.labels=20, bar.lwd=5) {

	N <- nrow(m)
	NC <- ncol(m)

	if ( all( is.na(m))) {
		cat( "\nWarning:  no non-zero gene intensities !!")
		cat( "\nPerhaps expression data does not match current species...")
		return(NULL)
	}

	par( "mai"=c( 1.25, 0.95, 0.85, 0.2))
	las <- 3
	border <- par( "fg")
	if ( is.null( col)) {
		col=gray( seq( 0.2, 1.0, length.out=N))
	} else {
		if ( N >= 20) border <- NA
	}

	barSpace <- c( 0, N/4)

	if ( is.null( yMax)) yMax <- max(m) * 1.1
	mainText <- paste( "Life Cycle Stage Plot:\n", label)

	mp <- barplot(m, beside=T, col=col, border=border, main=mainText, 
		ylab="Percent of Total Gene Intensity", 
		space=barSpace, las=las, font.lab=2, font.axis=2, cex.lab=1.05, cex.axis=1.05, 
		ylim=c(0,yMax), xlim=c( N*0.1,N*(ncol(m)+4.2)))

	splineAns <- NULL
	if ( is.logical( barSplines)) {
		if ( barSplines) {
			barSplines <- 1:N
		} else {
			barSplines <- vector()
		}
	}
	if ( length(barSplines)) {

		# use the "middle" X for all stages, so large Ncolumns doesn't distort the pic.
		x <- apply( mp, MARGIN=2, FUN=mean)
		splineAns <- vector()
		for( j in round( barSplines)) {
			y <- m[ j, ]
			# x <- mp[ j, ]
			ans <- spline( x, y, n=8*length(x), method="natural")
			lines( ans$x, ans$y, col=par("fg"), lwd=bar.lwd*1.4, lty=1)
			lines( ans$x, ans$y, col=col[j], lwd=bar.lwd, lty=1)
			ans2 <- spline( 1:NC, y, n=20*NC, method="natural")
			splineAns[j] <- ans2$x[ which.max(ans2$y)]
		}
		names(splineAns) <- rownames(m)
		splineAns <- splineAns[ ! is.na(splineAns)]
	}
	
	# limit the legend to a reasonable number
	who <- 1:N
	isSubset <- FALSE
	if ( N > 12) legend.cex <- legend.cex * 0.95
	if ( N > 18) legend.cex <- legend.cex * 0.95
	if ( N > max.labels) {
		who <- seq.int( 1, N, length.out=max.labels)
		isSubset <- TRUE
		legend.cex <- legend.cex * 0.95
	}
	ans <- legend( "topright", rownames(m)[who], fill=col[who], cex=legend.cex, bg="white")
	if ( isSubset) mtext( "not all labeled", side=3, adj=1, cex=legend.cex*0.95)
	dev.flush()

	return( splineAns)
}


`bestLifeCycleStageTimeCall` <- function( splineAns, timePoints=c(0,8,16,24,32,40,48)) {

	# time is currently hard coded for "Blood.Only" mode
	IE_MAX <- 48

	lowLimits <- floor( splineAns)
	lowTime <- timePoints[lowLimits]
	hiLimits <- ceiling( splineAns)
	hiTime <- timePoints[hiLimits]
	if ( any( is.na( c( lowTime, hiTime)))) {
		cat( "\nExpected non-NA time hour values for IE cycle stages")
		cat( "\nRequires a 'spline.maximums' object from a LifeCycle plot result.")
		return( NA)
	}

	fraction <- splineAns - lowLimits
	bestTime <- lowTime + fraction * (hiTime - lowTime)
	isNEG <- which( bestTime < 0)
	isBIG <- which( bestTime > IE_MAX)
	if ( length(isNEG)) bestTime[ isNEG] <- bestTime[isNEG] + IE_MAX
	if ( length(isBIG)) bestTime[ isBIG] <- bestTime[isBIG] - IE_MAX

	return( round( bestTime, digits=2))
}


`plotLifeCycleStageUnitVectors` <- function( gSet, col=1, lwd=1, legend=NA, plot=TRUE, yMax=1,
				legend.cex=1, label="", annotate=FALSE) {

	verifyLifeCycleSetup()

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- LifeCycleEnv[[ "VectorSpace"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	if (plot) {

		par( "mai"=c( 1.25, 0.95, 0.85, 0.4))

		plot( 1,1, type="n", main=paste( "Life Cycle Stage:   Gene Unit Vectors\n",label),
			xlim=c(0.5,N_STAGES+0.5), ylim=c(0,yMax), xaxt="n", xlab=NA,
			ylab="Gene 'Projection' per Stage", las=3, font.axis=2, font.lab=2)
		axis( 1, at=1:N_STAGES, label=STAGE_NAMES, las=3, font=2)
	}

	# draw the lines a bit nicer...as step steps...
	colUse <- rep( col, length.out=length(gSet))
	lwdUse <- rep( lwd, length.out=length(gSet))

	where <- base::match( gSet, vectorSpace$GENE_ID, nomatch=NA)
	who <- which( !is.na(where))

	for ( k in who) {
		y <- as.numeric( vectorSpace[ where[k], 3:ncol(vectorSpace)])
		ans <- drawStageDensityLine( y, col=colUse[k], lwd=lwdUse[k])
		xLocation <- ans$x
		if ( annotate) {
			useY <- which.max( y)
			text( xLocation[useY], y[useY], gSet[k], col=colUse[k], pos=3, font=2, cex=1)
		}
	}
	if ( ! is.na( legend)) {
		legend( x=legend, legend=gSet[who], col=colUse[who], lwd=3, cex=legend.cex)
	}
	
	# send back useful info
	units <- vectorSpace[ where, ]
	out <- list( "unitVectors"=units, "x"=xLocation)
	return( out)
}


`plotLifeCycleStageIntensity` <- function( gSet, col=1, lwd=1, pt.cex=1.25, legend=NA, plot=TRUE, 
				legend.cex=1, label="", minYmin=1, minYmax=NULL, threshold=NULL, annotate=FALSE) {

	verifyLifeCycleSetup()

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	intenSpace <- LifeCycleEnv[[ "IntensitySpace"]]
	intenVectors <- as.matrix( intenSpace[ , 3:ncol(intenSpace)])
	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	where <- base::match( gSet, intenSpace$GENE_ID, nomatch=0)
	if (sum(where > 0) < 1) {
		cat( "\nNo Matching Genes found..")
		return()
	}

	intenSpace <- intenSpace[ where, , drop=F]
	intenVectors <- intenVectors[ where, , drop=F]
	gSet <- gSet[ where > 0]
	NG <- length(gSet)

	bigValue <- max( intenVectors, na.rm=T)
	if ( ! is.null( minYmax)) bigValue <-max( bigValue, minYmax, na.rm=T)
	if ( ! is.null( threshold)) bigValue <- max( bigValue, threshold, na.rm=T)
	yMax <- bigValue * 2
	smallValue <- min( intenVectors, na.rm=T)
	if ( ! is.null( threshold)) minYmin <- min( minYmin, threshold)
	yMin <- min( smallValue, minYmin)
	if ( yMin < 0.01) yMin <- 0.01

	if (plot) {

		par( "mai"=c( 1.25, 0.95, 0.85, 0.4))
		mainText <- paste( "Life Cycle Stage:    Gene Intensity Vectors\n",label)
		if ( NG == 1) {
			mainText <- paste( "Life Cycle Stage:    Gene Intensity:    ", gSet)
			if ( (orig <- gene2OrigID(gSet)) != gSet) mainText <- paste( mainText, "    (",orig,")", sep="")
			mainText <- paste( mainText, "\n", gene2Product( gSet))
		}

		plot( 1,1, type="n", main=mainText, xlim=c(0.5,N_STAGES+0.5), ylim=c(yMin,yMax), 
			log="y", xaxt="n", xlab=NA, ylab="Gene Intensity per Stage  (RPKM)", las=3, font.axis=2, font.lab=2)
		axis( 1, at=1:N_STAGES, label=STAGE_NAMES, las=3, font=2)
	}

	# draw the lines a bit nicer...as step steps...
	colUse <- rep( col, length.out=length(gSet))
	lwdUse <- rep( lwd, length.out=length(gSet))

	if ( ! is.null( threshold)) {
		lines( c(-2,20), rep.int(threshold,2), lwd=1, lty=3, col="brown")
		text( c(0.3,0.3), c( threshold,threshold), c( "Expression", "Threshold"), col="brown", pos=c(3,1))
	}

	# when only one gene, color by stage instead
	if ( NG == 1) {
		colUse <- rainbow( ncol(intenVectors), end=0.78)
		y <- intenVectors[ 1, ]
		y <- pmax( y, yMin)
		ans <- drawStageIntensityLine( y, col=colUse, lwd=lwdUse[1], pt.cex=pt.cex, col2=1, useLog=T)
		xLocation <- ans$x
		if ( annotate) {
			useY <- which.max( y)
			text( xLocation[useY], y[useY], gSet[1], col=colUse[useY], pos=3, font=2, cex=1)
		}
	} else {
	    for ( k in 1:length(gSet)) {
		y <- intenVectors[ k, ]
		y <- pmax( y, yMin)
		ans <- drawStageIntensityLine( y, col=colUse[k], lwd=lwdUse[k], pt.cex=pt.cex, useLog=T)
		xLocation <- ans$x
		if ( annotate) {
			useY <- which.max( y)
			text( xLocation[useY], y[useY], gSet[k], col=colUse[k], pos=3, font=2, cex=1)
		}
	    }
	}
	if ( ! is.na( legend)) {
		if ( NG == 1) colUse <- 1
		legend( x=legend, legend=gSet, col=colUse, lwd=3, cex=legend.cex)
	}
	dev.flush()
	
	# send back useful info
	out <- list( "intensityVectors"=intenVectors, "x"=xLocation)
	return( out)
}


`drawStageDensityLine` <- function(y, col, lwd) {

	yUse <- c( 0, rep( y, each=2), 0)
	xmids <- seq( 0, length(y), by=1) + 0.5
	xUse <- base::sort( c( xmids-0.05, xmids+0.05))
	lines( xUse, yUse, type="l", col=col, lwd=lwd) 
	return( invisible( list( "x"=(xmids+0.5))))
}


`drawStageIntensityLine` <- function(y, at=1:length(y), col, lwd, pt.cex, col2=col, useLog=FALSE) {

	N <- length(y)
	if (useLog) {
		splineAns <- spline( at, log2(y+1), n=N*5)
		xSpline <- splineAns$x
		ySpline <- 2 ^ (splineAns$y) - 1
		# spline in log space can look extreme...
		minY <- min(y, na.rm=T) *.1
		ySpline[ ySpline < minY] <- minY
	} else {
		splineAns <- spline( at, y, n=N*5)
		xSpline <- splineAns$x
		ySpline <- splineAns$y
	}
	lines( xSpline, ySpline, col=col2, lty=2, lwd=lwd)
	points( at, y, pch=21, bg=col, cex=pt.cex)
	return( invisible( list( "x"=at)))
}


`adjustLifeCycleStageMatrix` <- function( geneSet, intenMatrix, fids=colnames(intenMatrix), 
		fcolors=NULL, tolerance=1.0, rate=1.0, max.iterations=100, 
		label="your label here...", legend.cex=1, max.labels=20, pause=FALSE) {

	verifyLifeCycleSetup()

	# turn the matrix into some temp transcript files, then call the file set tool
	nFiles <- ncol( intenMatrix)
	fNames <- vector()

	for ( i in 1:nFiles) {
		f <- paste( "matrix", fids[i], "Pf.Transcript.txt", sep=".")
		fNames[i] <- f
		tmp <- data.frame( "GENE_ID"=geneSet, "INTENSITY"=intenMatrix[ ,i], stringsAsFactors=F)
		write.table( tmp, file=f, sep="\t", quote=F, row.names=F)
	}

	ans <- adjustLifeCycleStageFileSet( fNames, fids, fcolors=fcolors, tolerance=tolerance, 
			rate=rate, max.iterations=max.iterations, label=label, legend.cex=legend.cex,
			max.labels=max.labels, pause=pause)

	# with this final adjustment scaling, modify the original matrix
	intenOut <- intenMatrix
	gPtrs <- match( rownames(intenMatrix), rownames(ans$scaleFactors), nomatch=0)
	for ( j in 1:nrow(intenMatrix)) {
		where <- gPtrs[j]
		if ( where == 0) next
		scal <- ans$scaleFactors[ where, ]
		if (any( scal == 0)) next
		intenOut[ j, ] <- intenMatrix[ j, ] * scal
	}

	# remove those temp files
	for ( i in 1:nFiles) {
		f <- paste( "matrix", fids[i], "Pf.Transcript.txt", sep=".")
		file.delete(f)
	}

	return( intenOut)
}


`adjustLifeCycleStageFileSet` <- function( fnames, fids, fcolors=NULL, geneColumn="GENE_ID", 
		intensityColumn="INTENSITY", tolerance=1.0, rate=1.0, max.iterations=100, 
		label="your label here...", legend.cex=1, max.labels=20, pause=FALSE, yMax=NULL) {

	verifyLifeCycleSetup()

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nLifeCycleStage:  Some transcript files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}

	# build local re-writable files and storage for each
	myFiles <- paste( "tmp.stageAdjusted", basename( fnames), "rda", sep=".")
	nFiles <- length( fnames)
	fileList <- geneIntenList <- vector( mode="list", length=nFiles)

	neededColumns <- c( geneColumn, intensityColumn)
	for( i in 1:nFiles) {
		tbl <- read.delim( fnames[i], as.is=T)
		tbl <- tbl[ , neededColumns]
		fileList[[i]] <- tbl
		save( tbl, file=myFiles[i])
	}

	# make one copy of all gene intensities before and after to create final multipliers
	startingInten <- endingInten <- matrix( 1.0, nrow=nrow(tbl), ncol=nFiles)
	rownames(startingInten) <- rownames(endingInten) <- geneNames <- base::sort( tbl[ , geneColumn])
	colnames(startingInten) <- colnames(endingInten) <- fids
	for( i in 1:nFiles) {
		tbl <- fileList[[i]]
		who <- base::match( geneNames, tbl[ , geneColumn])
		startingInten[ ,i] <- endingInten[ , i] <- tbl[ who, intensityColumn]
	}

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- LifeCycleEnv[[ "VectorSpace"]]
	N_STAGES <- LifeCycleEnv[[ "N_STAGES"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	STAGE_NAMES <- LifeCycleEnv[[ "STAGE_NAMES"]]

	# storage for the resulting stage histograms
	m <- matrix( nrow=nFiles, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# repeatedly measure and adjust til convergence
	nIters <- 0
	curDev <- N_STAGES
	repeat {

		# load each file in turn
		for( i in 1:nFiles) {
			tbl <- fileList[[ i]]
			ans <- calcLifeCycleStage( g=tbl[[ geneColumn]],
				inten=tbl[[ intensityColumn]])
			m[ i, ] <- ans$Stage
			geneIntenList[[ i]] <- ans$IntensityVectors
		}

		# plot it
		plotText <- paste( label, "\nIter: ", nIters, "     Max Intra-Stage Dif:  ", 
				formatC(curDev, format="f", digits=3))
		plotLifeCycleStages(m, col=fcolors, label=plotText, yMax=yMax, legend.cex=legend.cex,
				max.labels=max.labels)

		# calculate the average per stage and the factor set to apply to each file
		stageAvg <- apply( m, MARGIN=2, mean)
		deviation <- m
		for (i in 1:nFiles) deviation[ i, ] <- stageAvg - m[ i, ]
		netDifference <- range( as.vector( deviation))
		curDev <- netDifference[2] - netDifference[1]

		nIters <- nIters + 1
		cat( "\rIter: ", nIters, "\tMax Intra-Stage Dif: ", curDev)
		if ( curDev <= tolerance) break
		if ( nIters > max.iterations) break
		if (pause) locator(1)

		# now use these file to file stage deviations to tweak each transcript
		multiplier <- (100 + (deviation*rate)) / 100
		for( i in 1:nFiles) {
			tbl <- fileList[[ i]]
			newInten <- recombineLifeCycleStageIntensity( geneIntenList[[i]], multiplier[ i, ])

			# put this new intensities back for next pass
			tbl[ , intensityColumn] <- newInten
			fileList[[ i]] <- tbl
			save( tbl, file=myFiles[i])

			who <- base::match( geneNames, names( newInten))
			endingInten[ , i] <- newInten[ who]
		}
	}
	cat( "\nDone.")

	# plot it again to show latest numbers
	plotText <- paste( label, "\nDone.   Iter: ", nIters, "     Max Intra-Stage Dif:  ", 
			formatC(curDev, format="f", digits=3))
	plotLifeCycleStages(m, col=fcolors, label=plotText, yMax=yMax, legend.cex=legend.cex,
			max.labels=max.labels)

	finalScaling <- endingInten / startingInten
	finalScaling[ is.nan(finalScaling)] <- 1.0

	out <- list( "deviation"=curDev, "iterations"=nIters, "Stages"=m, 
			"scaleFactors"=finalScaling)

	# clean up the temp files...
	file.delete( myFiles)

	return( out)
}


`recombineLifeCycleStageIntensity` <- function( geneStageInten, multiplier=rep( 1.0, times=ncol(geneStageInten))) {

	# recombine the stage intensity portions into one gene intensity, allowing for
	# a multiplier to be applied to each stage
	myInten <- geneStageInten
	sumBefore <- sum( myInten, na.rm=T)

	for ( i in 1:ncol( geneStageInten)) myInten[ ,i] <- myInten[ ,i] * multiplier[i]

	sumAfter <- sum( myInten, na.rm=T)
	adjFactor <- sumBefore / sumAfter
	myInten <- myInten * adjFactor

	geneInten <- apply( myInten, MARGIN=1, sum, na.rm=TRUE)

	names(geneInten) <- rownames( geneStageInten)
	return( geneInten)
}

