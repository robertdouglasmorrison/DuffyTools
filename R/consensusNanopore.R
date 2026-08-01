# consensusNanopore.R - create consensus sequences from Nanopore style FASTQ datasets

# try #2:  start every sequence as its own cluster, then merge by declining score threshold

require( Biostrings)
require( pwalign)


consensusNanopore <- function( file, max.reads=500, down.sample=TRUE, final.pct.match=95) {
	
	# save the file name globally
	FASTQ.NAME <<- sub( "\\.f[aq](stq)?(\\.g)?$", "", basename(file))
	
	# read in the raw FASTQ reads
	seqs <- loadNanoporeFastq( file, max.reads=max.reads, down.sample=down.sample)
	
	# try to put them all into coding strand orientation
	seqs <- asCodingStrandSequence(seqs)
	
	# start with every sequence in its own cluster
	makeClusters( seqs)
	
	# show the initial distribution
	plotClusters()
	
	# use the DNA scoring matrix
	DNA.DM <<- nucleotideSubstitutionMatrix()
	
	# repeatedly try to merge clusters, lowering the score threshold as we go
	matchPct <- 100
	N.MERGES <<- 0
	repeat {
		mergeClusters( pct.match=matchPct)
		matchPct <- round( matchPct - 0.1, digits=2)
		if ( matchPct < final.pct.match) break
	}
	
	# when done turn each cluster with enough members into an MSA object
	MSAbyCluster( min.seq=10)
}


loadNanoporeFastq <- function( file, max.reads=500, down.sample=FALSE, verbose=T) {

	# read in the raw sequence data
	if ( ! file.exists( file)) stop( paste( "FASTQ data file not found: ", file))
	fq <- readFastq( file, verbose=verbose)
	ids <- fq$READ_ID
	seqs <- fq$READ_SEQ
	
	# # allow discarding very short & very long reads, to keep the total N reasonable
	N <- length( seqs)
	if ( N > max.reads) {
		if ( down.sample) {
			cat( "\nReducing raw read count to", max.reads, "by random sampling")
			keep <- sample( N, size=max.reads)
		} else {
			lens <- nchar(seqs)
			ord <- order( lens)
			nDropFront <- floor( (N-max.reads)/2)
			nDropRear <- N - max.reads - nDropFront
			keepPtrs <- ord[ (nDropFront+1) : (N-nDropRear)]
			cat( "\nReducing raw read count to", max.reads, "by discarding", nDropFront, "shortest &", nDropRear, "longest sequneces")
			keep <- which( 1:N %in% keepPtrs)
		}
		seqs <- seqs[keep]
		ids <- ids[keep]
	}
	
	# final data will be DNAStringSet object
	out <- DNAStringSet( seqs)
	names(out) <- ids
	return(out)
}


asCodingStrandSequence <- function(seqs) {

	# try to guess the coding strand of each, and reverse complement as needed
	out <- seqs
	nRC <- 0
	for ( i in 1:length(seqs)) {
		dna <- seqs[[i]]
		aa <- DNAtoBestPeptide( as.character(dna))
		if ( names(aa) %in% c( "FR4", "FR5", "FR6")) {
			out[[i]] <- reverseComplement( dna)
			nRC <- nRC + 1
		}
	}
	if ( nRC > 0) cat( "\nSome sequences got RevComp to predicted coding strand. N=", nRC)
	return(out)
}


makeClusters <- function( seqs) {

	# set up global storage for the sequence clusters
	N.Seq <<- N.Clusters <<- N <- length(seqs)
	Sequences <<- seqs
	SequenceLens <<- nchar(seqs)
	
	# start every cluster with one sequence
	ClustMembers <<- vector( mode="list", length=N)
	ClustSize <<- rep.int( 1, N)
	for ( i in 1:N) ClustMembers[[i]] <<- i
	
	# give every cluster a plot location
	Clust.X <<- Clust.Y <<- rep.int(1,N)
	NR <- NC <- ceiling( sqrt( N))
	x <- y <- 1
	for ( i in 1:N) {
		Clust.X[i] <<- x
		Clust.Y[i] <<- y
		x <- x + 1
		if (x > NR) {
			x <- 1
			y <- y + 1
		}
	}
	
	# set up the Score Matrix
	ScoreMatrix <<- matrix( NA, nrow=N, ncol=N)
	
	return(N)
}


plotClusters <- function( scaleFactor=0.25, label="") {
	
	require( plotrix)
	maxColors <- round(N.Seq/2)
	colRamp <- rev( rainbow( maxColors, end=0.7))
	
	# make dots for every cluster, proportional to cluster size
	NR <- NC <- ceiling( sqrt( N.Seq))
	plot( 1, 1, type="n", main=paste( FASTQ.NAME, label, sep="\n"), xlim=c(0,NR+1), ylim=c(0,NC+1), 
			xlab="Clusters", ylab="Clusters")

	# show all that are any non-zero size irst			
	toShow <- which( ClustSize > 0)
	points( Clust.X[toShow], Clust.Y[toShow], pch='.', col=colRamp[1])
	
	for ( i in toShow) {
		mySize <- ClustSize[i]
		if ( mySize > 1) {
			myRadius <- sqrt(mySize) / pi
			myCol <- colRamp[ min( mySize, maxColors)]
			draw.circle( Clust.X[i], Clust.Y[i], radius=(myRadius*scaleFactor), border=1, col=myCol)
		}
	}
	dev.flush()
	return(NULL)
}


mergeClusters <- function( pct.match=99) {	

	nClust <- sum( ClustSize > 0)
	cat( "\nStarting Merges of ", nClust, "Sequence Clusters at ", pct.match,"% match threshold..\n")
	# set up to visit all the clusters in random order, to find other clusters of highly similar sequences
	visitOrder <- sample( N.Clusters)
		
	# visit each one, and count how many merges occur
	nMerges <- 0
	for ( i in 1:N.Clusters) {
		myClust <- visitOrder[i]
		mySize <- ClustSize[ myClust]
		myMembers <- ClustMembers[[ myClust]]
		if ( mySize < 1) next
		myPtr <- myMembers[1]
		myLen <- SequenceLens[myPtr]
		
		# get the set of all other clusters to compare against
		otherClusters <- setdiff( which( ClustSize > 0), myClust)
		N.Others <- length( otherClusters)
		# get the first member of each cluster
		otherPtrs <- sapply( ClustMembers[ otherClusters], `[`, 1)
		otherLens <- SequenceLens[otherPtrs]
		
		# get the pairwise overlap scores between all the others and me
		paScore <- getPairwiseAlignScores( otherPtrs, myPtr)
		
		# the Scoring Matrix is +1 for each match, so the best possible overlap score is the length of the shorter sequence
		maxPossibleScore <- pmin( otherLens, rep.int( myLen, N.Others))
		scoreThreshold <- maxPossibleScore * (pct.match / 100)
		
		# decide which other clusters are good enough to join this one
		goodHits <- which( paScore >= scoreThreshold)
		if ( ! length( goodHits)) next
		
		# when we get 2+ good matches, decide how/if/which to keep
		finalGoodHits <- goodHits
		if ( length( goodHits) > 1) {
			goodOtherLens <- otherLens[ goodHits]
			# if all the other clusters are shorter than my seq, then they all stay.
			# if all are longer than my seq, then only keep the shortest one
			# if some are shorter and others are longer, only keep the shorter ones
			if ( all( goodOtherLens >= myLen)) {
				bestOne <- which.min( abs( goodOtherLens - myLen))
				finalGoodHits <- goodHits[ bestOne]
			} else if ( any( goodOtherLens <= myLen)) {
				bestOnes <- which( goodOtherLens <= myLen)
				finalGoodHits <- goodHits[ bestOnes]
			}
		}
		
		# gather all the members of all the clusters to merge.  Start from mine and add those from each hit.
		clusterPtrSet <- myClust
		seqPtrSet <- myMembers
		for ( k in finalGoodHits) {
			otherClust <- otherClusters[k]
			clusterPtrSet <- c( clusterPtrSet, otherClust)
			otherMembers <- ClustMembers[[ otherClust]]
			seqPtrSet <- c( seqPtrSet, otherMembers)
		}
		# rearrange the members to keep the longest one first within the new cluster
		seqSize <- SequenceLens[ seqPtrSet]
		seqOrd <- order( seqSize, decreasing=T)
		seqPtrSet <- seqPtrSet[seqOrd]
		
		# update this cluster to contain its new members
		ClustSize[ myClust] <<- length( seqPtrSet)
		ClustMembers[[ myClust]] <<- seqPtrSet
		# and then remove/erase all the others
		for ( k in clusterPtrSet[ 2:length(clusterPtrSet)]) {
			ClustSize[ k] <<- 0
			ClustMembers[[ k]] <<- integer(0)
		}
		nMerges <- nMerges + length(goodHits)
		N.MERGES <<- N.MERGES + length(goodHits)
		# done with this cluster
		cat( "\r", i, myClust, nMerges)
		# redraw every time we merge
		plotClusters( label=paste( "N.Merges:", N.MERGES, "  Pct.Match:", pct.match))
	}
	return(nMerges)
}


getPairwiseAlignScores <- function( otherSeqPtrs, mySeqPtr) {

	# given a vector of pointers to sequences that are the 'pattern' seqs and one pointer to the 'subject' seq
	
	# look up those scores in the Score Matrix, to see which need to be calculated vs already exist
	scores <- ScoreMatrix[ otherSeqPtrs, mySeqPtr]
	missing <- which( is.na( scores))
	if ( ! length(missing)) return( scores)
	
	# we need to calculate some.  Get the sequences and do the PA
	otherSeqs <- Sequences[ otherSeqPtrs[missing]]
	mySeq <- Sequences[[ mySeqPtr]]
	
	# measure the pairwise overlap of these
	cat( "\rPairwise Align Scores to calculate=", length(missing), "for seq #", mySeqPtr)
	paScore <- pairwiseAlignment( otherSeqs, mySeq, type="overlap", scoreOnly=T, substitutionMatrix=DNA.DM, gapOpening=4, gapExtension=4)
	
	# stuff these new scores back into the score matrix
	ScoreMatrix[ otherSeqPtrs[missing], mySeqPtr] <<- paScore
	
	# now re-extract and send them back
	scores <- ScoreMatrix[ otherSeqPtrs, mySeqPtr]
	return( scores)
}


MSAbyCluster <- function( min.seq=5) {

	# make DNA Multiple Sequence alignments of each cluster, that has enough members
	hasEnoughMembers <- which( ClustSize >= min.seq)
	
	# visit them by decreasing cluster size
	cSize <- ClustSize[ hasEnoughMembers]
	visitOrd <- order( cSize, decreasing=T)
	hasEnoughMembers <- hasEnoughMembers[visitOrd]
	for ( i in 1:length(hasEnoughMembers)) {
		k <- hasEnoughMembers[i]
		myPtrs <- ClustMembers[[ k]]
		mySeqs <- Sequences[ myPtrs]
		myIDs <- names(Sequences)[myPtrs]
		fa <- as.Fasta( myIDs, mySeqs)
		faname <- paste( "Cluster", i, "fasta", sep=".")
		writeFasta( fa, faname, line=100)
		alnname <- paste( "Cluster", i, "aln", sep=".")
		aln <- mafft( faname, alnname, mode="genaffine")
		aln$alignment <- toupper( aln$alignment)
		writeALN( aln, alnname, line=100)
		consensusDNA <- consensusAlignment( aln$alignment)
		consensusAA <- DNAtoBestPeptide( consensusDNA)
		pctOfSeqs <- round( length(myPtrs) * 100 / N.Seq, digits=1)
		faDNA <- as.Fasta( paste( "Cluster.", i, ".DNA ", "Pct of Reads=", pctOfSeqs, sep=""), consensusDNA)
		writeFasta( faDNA, paste( "Cluster", i, "Consensus.DNA.fasta", sep="."), line=100)
		faAA <- as.Fasta( paste( "Cluster.", i, ".AA ", "Pct of Reads=", pctOfSeqs, sep=""), consensusAA)
		writeFasta( faAA, paste( "Cluster", i, "Consensus.AA.fasta", sep="."), line=100)
	}
}


consensusAlignment <- function( m) {

	# given a matrix of a MSA alignment, select the best single call at every location
	
	# step 1: remove gap symbols from the ends of each sequence, but keep all internal gaps
	for ( i in 1:nrow(m)) {
		j <- 1
		while( m[i,j] == "-") { m[i,j] <- ""; j <- j + 1}
		j <- ncol(m)
		while( m[i,j] == "-") { m[i,j] <- ""; j <- j - 1}
	}
	
	# step 2: get the most common base call at each column
	chAns <- apply( m, 2, function(x) {
		x <- x[ x != ""]
		xt <- sort( table( x), decreasing=T)
		return( names(xt[1]))
	})
	# step 2: remove the ends that are dominated by gaps
	startFrom <- 1
	while (chAns[startFrom] %in% c( "", " ", "-")) startFrom <- startFrom + 1
	endAt <- length(chAns)
	while (chAns[endAt] %in% c( "", " ", "-")) endAt <- endAt - 1
	out <- paste( chAns[ startFrom : endAt], collapse="")
	# remove any embedded gaps
	out <- gsub( "-", "", out, fixed=T)
	return( out)
}
