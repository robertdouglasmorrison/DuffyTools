# nanoporeClusters.R - create cluster-based consensus sequences from a Nanopore style FASTQ dataset

# Overview:  start every sequence as its own cluster, then repeatedly merge by 
#			declining pairwise similarity score threshold, to find the most common 
#			sequences in a raw dataset.

# set up the tool as a closure, to allow encapulation and access to the various parts


nanoporeClusters <- function( file, max.reads=1000, expected.size=NULL, expected.aa.seq=NULL, final.pct.match=97, 
							trim=0.01, prefix=sub("\\.f[aq](stq)?(\\.gz)?$","",basename(file)), 
							results.path=paste( prefix, "NanoporeClusters",sep="."),
							min.seq.per.cluster=5, min.pct.per.cluster=5, verbose=TRUE, plot=TRUE) {

	require( Biostrings)
	require( pwalign)
	require( plotrix)
	checkX11( width=9, height=8)
	

	# use the Biostrings DNA scoring matrix
	DNA.DM <- nucleotideSubstitutionMatrix()

	# define local variables that are globally available inside the closure
	FILE <- file
	MAX.READS <- max.reads
	EXPECT.SIZE <- expected.size
	EXPECT.AA.SEQ <- expected.aa.seq
	FINAL.PCT.MATCH <- final.pct.match
	PREFIX <- gsub( " +", ".", prefix)
	RESULTS.PATH <- results.path
	TRIM <- trim
	MIN.SEQ.PER.CLUSTER <- min.seq.per.cluster
	MIN.PCT.PER.CLUSTER <- min.pct.per.cluster
	VERBOSE <- verbose
	PLOT <- plot
	
	# storage for the Sequences that are members of clusters
	N.Seq <- integer(1)
	Sequences <- character(1)
	SequenceLens <- integer(1)
	# storage for the clusters
	N.Clusters <- integer(1)
	ClustMembers <- vector( mode="list")
	ClustSize <- integer(1)
	Clust.X <- Clust.Y <- numeric(1)
	
	# storage for final MSA objects
	MSAs <- vector( mode="list")
	ConsensusCalls <- vector( mode="list")

	# storage for the pairwise alignment scoring matrix
	ScoreMatrix <- matrix( NA, nrow=1, ncol=1)
	# other storage for convenience 
	seqsRaw <- seqsNet <- character(1)

	# a global counter for showing progress
	N.MERGES <- 0


	# top level function to do everything
	do.all <- function() {
	
		# step 1: load the reads, down-select, and build initial clusters
		setup()
		# step 2: merge exact duplicates
		nMergesD <- mergeDuplicates()	
		N.MERGES <<- N.MERGES + nMergesD
		# step 3: merge perfect substrings
		nMergesS <- mergePerfectSubstrings()	
		N.MERGES <<- N.MERGES + nMergesS
		# step 4:  the main merging step, to iteratively find merges as the score threshold comes down
		N.MERGES <<-  mergeClusters( final.pct.match=FINAL.PCT.MATCH)
			
		# when done, turn each cluster with enough members into an MSA object
		MSAbyCluster( min.seq.per.cluster=MIN.SEQ.PER.CLUSTER, min.pct.per.cluster=MIN.SEQ.PER.CLUSTER, prefix=PREFIX)
		
		# final step: if 2+ clusters, check their similarity
		MultiClusterMSA()
	}


	setup <- function() {
	
		# read in the raw FASTQ reads
		seqsRaw <<- loadNanoporeFastq( file=FILE, max.reads=NULL)
		
		# create the folder where results will land. and if needed, remove older previous results
		if ( ! exists( RESULTS.PATH)) dir.create( RESULTS.PATH, recursive=T, showWarnings=F)
		oldFiles <- dir( RESULTS.PATH, pattern="Cluster\\.[0-9]+\\.(Consensus.AA.fasta|Consensus.DNA.fasta|aln|fasta|Diversity.csv)$", full=T)
		if ( length( oldFiles)) file.delete( oldFiles)

		# down select, giving preference to sequences seen 2+ times
		seqsNet <<- downSampleSequences( seqs=seqsRaw, max.reads=MAX.READS, expected.size=EXPECT.SIZE, trim=TRIM)
		
		# show what we kept/discarded
		plotSeqLenHistogram( seqsRaw, seqsNet)	
		
		# try to put them all into coding strand orientation
		seqsNet <<- asCodingStrandSequence( seqsNet)
	
		# start with every sequence in its own cluster
		makeClusters( seqsNet)
		
		# perhaps show it
		if (PLOT) { 
			plotClusters( label=paste( "Starting with:", N.Seq));  
			Sys.sleep(1)
			printPlot( file.path( RESULTS.PATH, paste(prefix,"Starting.Clusters",sep=".")), width=9, height=8)
		}
	}


	loadNanoporeFastq <- function( file, max.reads=NULL, verbose=T) {

		# read in the raw sequence data
		if ( ! file.exists( file)) stop( paste( "FASTQ data file not found: ", file))
		fq <- readFastq( file, verbose=verbose)
		ids <- fq$READ_ID
		seqs <- fq$READ_SEQ
		N <- length( seqs)
		# allow down sampling to reduce N
		if ( ! is.null( max.reads) && N > max.reads) {
			if (VERBOSE) cat( "\nReducing raw read count to", max.reads, "by random sampling")
			keep <- sample( N, size=max.reads)
			seqs <- seqs[keep]
			ids <- ids[keep]
		}
		# final data will be DNAStringSet object
		out <- DNAStringSet( seqs)
		names(out) <- ids
		return(out)
	}


	asCodingStrandSequence <- function(seqs=seqsNet) {
	
		# try to guess the coding strand of each, and reverse complement as needed
		out <- seqs
		nRC <- 0
		if (VERBOSE) cat( "\nGuessing at coding strand..")
		for ( i in 1:length(seqs)) {
			dna <- seqs[[i]]
			aa <- DNAtoBestPeptide( as.character(dna))
			if ( names(aa) %in% c( "FR4", "FR5", "FR6")) {
				out[[i]] <- reverseComplement( dna)
				nRC <- nRC + 1
			}
		}
		if ( VERBOSE && nRC > 0) cat( "  Some sequences got RevComp to coding strand. N=", nRC)
		return(out)
	}


	downSampleSequences <- function( seqs=seqsRaw, max.reads=1000, expected.size=NULL, trim=0.05) {
	
		# as we down sample to the max count of starting clusters, find sequences that deserve to be discarded or kept
		N <- N.Raw <- length(seqs)	
		slen <- nchar( seqs)
		if ( is.null( expected.size)) {
			expected.size <- round( mean( slen))
			if (VERBOSE) cat( "\nEstimating expected size of raw reads.  Expected =", expected.size,"bp")
		}

		if (VERBOSE) cat( "\nReducing sequences down to", max.reads, "by: ")
		
		# test1:  allow discarding very short & very long reads
		if ( ! is.null(trim)) {
			if ( trim > 0.25) stop( "'trim' value must be below 0.25")
			if ( max.reads > N.Raw) {
				trim <- trim * 0.1
			} else if ( N.Raw < max.reads * (1+trim*2)) {
				trim <- (N.Raw-max.reads) / (N.Raw*2)
			}
			if (VERBOSE) cat( "  Trim % =", trim)
			ord <- order( nchar(seqs))
			from <- round( N * trim)
			to <- round( N * (1-trim))
			keep <- ord[ from : to]
			if (VERBOSE) cat( "  Trimming", (from-1), "reads from each end of length distribution")
			seqs <- seqs[ keep]
			slen <- slen[ keep]
			N <- length( seqs)
			if (VERBOSE) cat("  N_Now: ", N)
		}
		
		# we used to favor perfect duplicates and substrings, but this idea seemed to favor
		# short unwanted artifacts.  Just do random sampling.
		if ( N < max.reads) return( seqs)
		
		# use the expected size to create a normal sampling pattern.  Generate a list of random sizes around the expected size,
		# using a growing wider Std Deviation, to favor the right size but allow outliers
		randLens <- vector()
		for ( mySD in seq( 0.1, 3, by=0.10)) {
			smlRandLens <- round( (rnorm( N, mean=0, sd=mySD) * expected.size / 3) + expected.size)  # scale by 100% of values inside +/-3*SD
			randLens <- c( randLens, smlRandLens)
		}
		drop1 <- which( randLens < min(slen))
		drop2 <- which( randLens > max(slen))
		drops <- sort( unique( c( drop1, drop2)))
		if ( length(drops)) randLens <- randLens[ -drops]
		avail <- 1:N
		chosen <- vector()
		nnow <- 0
		for ( i in 1:length(randLens)) {
			who <- match( randLens[i], slen[avail])
			if ( is.na(who)) next
			nnow <- nnow + 1
			chosen[nnow] <- avail[who]
			avail <- avail[ -who]
			if ( nnow >= max.reads) break
		}
		# we have our final set
		finalSet <- chosen
		if (VERBOSE) cat( "  Normal Random sampling: ", length(finalSet))
		return( seqs[finalSet])
	}


	plotSeqLenHistogram <- function( seqsRaw=seqsRaw, seqsNet=seqsNet) {

		# make a histogram of the sequence lengths, to show effect of trimming and selection
		ncRaw <- nchar( seqsRaw)
		ncNet <- nchar( seqsNet)
		Nraw <- length(seqsRaw)
		Nnet <- length(seqsNet)	
		# set the breaks value to create ~50bp sized bins
		breaksRaw <- max( round( diff( range(ncRaw)) / 100), 30)
		breaksNet <- max( round( diff( range(ncNet)) / 100), 5)
		# draw the original distribution
		hist( ncRaw, breaks=breaksRaw, main=paste( "Nanopore sample: ", PREFIX, "\nHistogram of raw sequence read lengths"),
				xlab="Sequence length (bp)", ylab="Number of sequences", col='grey70')
		# then lay on top the retained distribution
		hist( ncNet, breaks=breaksNet, add=T, col='green3')
		legend( 'topleft', paste( c("Original","Selected"), "  (N=", c(Nraw,Nnet), ")", sep=""), fill=c('grey70','green3'), cex=1.05)
		dev.flush(); Sys.sleep(1)
		printPlot( file.path( RESULTS.PATH, paste(prefix,"Sequence.Length.Histograms",sep=".")), width=9, height=7)
		return(NULL)
	}


	makeClusters <- function( seqs=seqsNet) {

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


	mergeDuplicates <- function() {
		
		# first pass to merge exact duplicates
		nMergesD <- combineDuplicateClusters()
		if (VERBOSE) cat( "\nN.Exact.Duplicates: ", nMergesD)
		if (PLOT) { 
			plotClusters( label=paste( "N.Duplicates:", nMergesD));  
			Sys.sleep(1)
			printPlot( file.path( RESULTS.PATH, paste(prefix,"After.Exact.Duplicates",sep=".")), width=9, height=8)
		}
		return( nMergesD)
	}
	
	
	mergePerfectSubstrings <- function() {
	
		# second pass to merge perfect substrings
		nMergesS <- perfectSubstringClusters()
		if (VERBOSE) cat( "\nN.Perfect.Substrings: ", nMergesS)
		if (PLOT) { 
			plotClusters( label=paste( "N.Perfect.Substrings:", nMergesS));  
			Sys.sleep(1) 
			printPlot( file.path( RESULTS.PATH, paste(prefix,"After.Perfect.Substrings",sep=".")), width=9, height=8)
		}
		return( nMergesS)
	}
		

	mergeClusters <- function( final.pct.match=FINAL.PCT.MATCH, start.pct.match=99.8) {
		
		# repeatedly try to merge clusters, lowering the score threshold as we go
		# because we catch all the perfect duplicates and proper substrings already, start the
		# merging process 1 step below perfect.  A 1bp mismatch in a 500bp sequence would be a 99.8% match
		cur.pct.match <- start.pct.match
		repeat {
			ans <- findAndMergeClusters( pct.match=cur.pct.match, random.order=F)
			if ( is.null(ans)) break
			cur.pct.match <- round( cur.pct.match - 0.1, digits=2)
			if ( cur.pct.match < final.pct.match) break
		}
		if (PLOT) { 
			finalMergeCnt <- N.MERGES
			nOrphans <- sum( ClustSize == 1)
			plotClusters( label=paste( "Final Clustering:   N.Seq:", N.Seq, "   N.Merge:", finalMergeCnt, "   N.Orphan:", nOrphans));  
			Sys.sleep(1) 
			printPlot( file.path( RESULTS.PATH, paste(prefix,"Final.Clusters",sep=".")), width=9, height=8)
		}
		return(N.MERGES)
	}
	
	
	combineDuplicateClusters <- function() {

		# do a pre-pass, to find exact duplicate sequences, so they can become clusters right away
		isDUP <- which( duplicated( Sequences))
		if ( ! length(isDUP)) return(0)	
		# visit each that has 2+ copies, and merge them into clusters immediately
		nMerged <- 0
		alreadyMerged <- vector()
		for ( i in isDUP) {
			mySeq <- Sequences[i]
			allHits <- which( Sequences == mySeq)
			newHits <- setdiff( allHits, alreadyMerged)
			nNew <- length( newHits)
			if ( nNew < 2) next		
			# we have a set of 2+ identical sequences. Put them all into the first
			myClust <- newHits[1]
			ClustSize[ myClust] <<- nNew
			ClustMembers[[ myClust]] <<- newHits		
			# then remove them from the other locations
			for ( k in newHits[2:nNew]) {
				ClustSize[k] <<- 0
				ClustMembers[[k]] <<- integer(0)
			}	
			# do not let these get seen again
			alreadyMerged <- c( alreadyMerged, newHits)
			nMerged <- nMerged + (nNew - 1)
		}
		return( nMerged)
	}


	perfectSubstringClusters <- function() {

		# find exact substring sequences, so they can become clusters right away
		useClust <- which( ClustSize > 0)
		clustSeqPtrs <- sapply( ClustMembers[useClust], FUN=`[`, 1)
		allClustSeqs <- Sequences[ clustSeqPtrs]
		
		nMerged <- 0	
		for ( i in 1:N.Clusters) {
			mySize <- ClustSize[i]
			if ( ! mySize) next
			myMembers <- ClustMembers[[i]]
			mySeqPtr <- myMembers[1]
			mySeq <- Sequences[ mySeqPtr]
			# see which clusters are a superstring of my seq
			myHits <- grep( mySeq, allClustSeqs)
			myClustHits <- useClust[ myHits]
			# remove myself from the set of hits
			myClustHit <- setdiff( myClustHits, mySeqPtr)
			if ( ! length( myClustHit)) next
			# we found a case where this sequence is an exact substring of some other cluster.
			# merge me into the shortest other cluster
			if ( length( myClustHit) > 1) {
				hitSeqLens <- nchar( Sequences[ myClustHit])
				myClustHit <- myClustHit[ which.min( hitSeqLens)]
			}
			# we have a cluster that is the best (shortest) superstring of me.  Add me to that one
			clustPtr <- myClustHit
			clustMembers <- c( ClustMembers[[ clustPtr]], myMembers)
			# no need to reorder by size, as that is already truly in size order
			nNew <- length( clustMembers)
			ClustSize[ clustPtr] <<- nNew
			ClustMembers[[ clustPtr]] <<- clustMembers	
			# then remove them from the current cluster
			ClustSize[i] <<- 0
			ClustMembers[[i]] <<- integer(0)	
			# since the clusters changed a bit, remake the set of all
			useClust <- which( ClustSize > 0)
			clustSeqPtrs <- sapply( ClustMembers[useClust], FUN=`[`, 1)
			allClustSeqs <- Sequences[ clustSeqPtrs]	
			nMerged <- nMerged + 1
		}
		return( nMerged)
	}


	plotClusters <- function( scaleFactor=0.45, label="") {
	
		# show the progress of clustering, by showing size of each cluster changing over time
		maxColors <- 100
		colRamp <- rev( rainbow( maxColors, end=0.6))
		# make dots for every cluster, proportional to cluster size
		NR <- NC <- ceiling( sqrt( N.Seq))
		plot( 1, 1, type="n", main=paste("Nanopore sample: ", PREFIX, "\n", label, sep=" "), xlim=c(0,NR+1), ylim=c(0,NC+1), 
				xlab="Cluster (X axis location)", ylab="Cluster (Y axis location)")
		# show all that are any non-zero size first			
		toShow <- which( ClustSize > 0)
		points( Clust.X[toShow], Clust.Y[toShow], pch='.', col=colRamp[1], cex=3)	
		# now make circles for the larger cluster, use color to show the percentage of the clustered sequences
		# draw from small to big, to prevent visual overlap
		ord <- order( ClustSize[toShow])
		nClust2 <- sum( ClustSize)
		for ( i in toShow[ord]) {
			mySize <- ClustSize[i]
			if ( mySize > 1) {
				myRadius <- sqrt(mySize) / pi
				myPct <- mySize / nClust2
				myColIndex <- round( sqrt(myPct) * 100)
				myCol <- colRamp[ myColIndex]
				draw.circle( Clust.X[i], Clust.Y[i], radius=(myRadius*scaleFactor), border=1, col=myCol)
				if ( mySize >= 10) text( Clust.X[i], Clust.Y[i], paste( round(mySize*100/N.Seq), "%", sep=""), cex=0.75)
			}
		}
		dev.flush(); 
		return(NULL)
	}


	findAndMergeClusters <- function( pct.match=99, random.order=TRUE) {	

		# main workhorse function: given a similarity threshold, combine clusters that are good enough matches
		nClust <- sum( ClustSize > 0)
		if ( nClust < 2) return(NULL)
		if (VERBOSE) cat( "\nStarting Merges of ", nClust, "Sequence Clusters at ", pct.match,"% match threshold..")
	
		# set up to visit all the clusters in random order, to find other clusters of highly similar sequences
		if ( random.order) {
			visitOrder <- sample( N.Clusters)
		} else {
			# since we start with some clusters having 2+ members due to duplicates & perfect substrings, 
			# so allow using size as the ordering preference
			visitOrder <- order( ClustSize, sample(N.Clusters), decreasing=T)
		}
		
		# visit each one, and count how many merges occur
		nMerges <- 0
		nRecentTests <- 0
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

			# the default score threshold was based on the shorter full length sequence, 
			# like when one is a true subset of the other.  But we can also have cases of a
			# true overlap, where both sequences extend past each other.  In this case, the length of the 
			# actual overlap is needed, to decide if the alignment is good enough to keep
			if ( ! length( goodHits)) {
				# see which one other had the very best match, perhaps redo the full PA to know the length of the overlap
				# if the score was somewhat close
				bestBad <- which.max( paScore)
				if ( paScore[bestBad] < (scoreThreshold[bestBad]*0.75)) next
				bestOtherPtr <- otherPtrs[bestBad]
				pa <- pairwiseAlignment( Sequences[bestOtherPtr], Sequences[myPtr], type="overlap", scoreOnly=F, 
									substitutionMatrix=DNA.DM, gapOpening=4, gapExtension=4)
				lenOverlap <- max( width(pattern(pa)), width(subject(pa)))
				scorePA <- score(pa)
				scorePAthreshold <- lenOverlap * (pct.match / 100)
				if( scorePA >= scorePAthreshold) goodHits <- bestBad[1]
			}
			nRecentTests <- nRecentTests + 1
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
			# done with this cluster merge
			#if (VERBOSE) cat( "\r", i, myClust, nMerges)
			# redraw every time we merge, unless we just redrew recently
			if ( nRecentTests < 10) next
			nRecentTests <- 0
			if (PLOT) {
				plotClusters( label=paste( "N.Cluster.Merge:", N.MERGES, "    Pct.Match.Threshold:", pct.match, "%"))
				Sys.sleep(0.01) 
				printPlot( file.path( RESULTS.PATH, paste(prefix,"Clustering.Progress",sep=".")), width=9, height=8)
			}

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
		if (VERBOSE) cat( "\rPairwise Align Scores to calculate=", length(missing), "for seq #", mySeqPtr)
		paScore <- pairwiseAlignment( otherSeqs, mySeq, type="overlap", scoreOnly=T, substitutionMatrix=DNA.DM, gapOpening=4, gapExtension=4)
		# stuff these new scores back into the score matrix
		ScoreMatrix[ otherSeqPtrs[missing], mySeqPtr] <<- paScore
		# now re-extract and send them back
		scores <- ScoreMatrix[ otherSeqPtrs, mySeqPtr]
		return( scores)
	}


	MSAbyCluster <- function( min.seq.per.cluster=5, min.pct.per.cluster=5, prefix=PREFIX) {

		# make DNA Multiple Sequence alignments of every cluster, that has enough members
		hasEnoughMembers <- which( ClustSize >= min.seq.per.cluster)
		# visit them by decreasing cluster size
		cSize <- ClustSize[ hasEnoughMembers]
		visitOrd <- order( cSize, decreasing=T)
		hasEnoughMembers <- hasEnoughMembers[visitOrd]
		# for tallying percentage of used sequence, only count those with 2+ seqs
		totalClustSeqs <- sum( ClustSize)
		nout <- 0	
		for ( i in 1:length(hasEnoughMembers)) {
			k <- hasEnoughMembers[i]
			mySize <- ClustSize[k]
			myPct <- round( mySize * 100 / totalClustSeqs, digits=1)
			if ( myPct < min.pct.per.cluster) next
			if (VERBOSE) cat( "\rMaking MSA for cluster:", k, "  Size:", mySize)
			# part 1: calculate the MSA 
			makeOneClusterMSA( clustPtr=k, clustID=i, prefix=prefix)
			# then extract the consensus info
			ans <- extractMSAconsensus( clustID=i)
			consensusDNA <- ans$seq[1]
			consensusAA <- DNAtoBestPeptide( consensusDNA)	
			myDepth <- ans$depth
			desc <- paste( prefix, "_Cluster", i, "_Reads=", mySize, "_Depth=", myDepth, "_Pct=", myPct, sep="")
			faDNA <- as.Fasta( desc, consensusDNA)
			writeFasta( faDNA, file.path( RESULTS.PATH, paste( prefix, "Cluster", i, "Consensus.DNA.fasta", sep=".")), line=100)
			truncWarn <- if ( nchar(consensusAA) < nchar(consensusDNA)/5) "  Translation Warning:  AA sequence is excessively truncated" else ""
			faAA <- as.Fasta( paste( desc, truncWarn, sep=" "), consensusAA)
			writeFasta( faAA, file.path( RESULTS.PATH, paste( prefix, "Cluster", i, "Consensus.AA.fasta", sep=".")), line=100)
			diversity <- ans$diversity
			write.csv( diversity, file.path( RESULTS.PATH, paste( prefix, "Cluster", i, "Diversity.csv", sep=".")), row.names=F)
			nout <- nout + 1
		}
		return( nout)
	}


	makeOneClusterMSA <- function( clustPtr, clustID=clustPtr, prefix) {

		# do the MSA for one cluster
		myPtrs <- ClustMembers[[ clustPtr]]
		mySeqs <- Sequences[ myPtrs]
		myIDs <- names(Sequences)[myPtrs]
		# make the FASTA of every cluster element sequence
		fa <- as.Fasta( myIDs, mySeqs)
		faname <- file.path( RESULTS.PATH, paste( prefix, "Cluster", clustID, "fasta", sep="."))
		writeFasta( fa, faname, line=100)
		# do the MSA step
		alnname <- file.path( RESULTS.PATH, paste( prefix, "Cluster", clustID, "aln", sep="."))
		aln <- mafft( faname, alnname, mode="genaffine")
		aln$alignment <- toupper( aln$alignment)
		writeALN( aln, alnname, line=100)
		#stuff this MSA object into global storage
		MSAs[[ clustID]] <<- aln
		names(MSAs)[clustID] <<- basename( alnname)
		return( length( mySeqs))
	}


	extractMSAconsensus <- function( clustID) {

		# get the MSA for one cluster, to extract the consensus DNA string
		aln <- MSAs[[ clustID]] 
		ans <- consensusAlignment( aln$alignment)
		#stuff a copy of this Consensus object into global storage
		ConsensusCalls[[ clustID]] <<- ans
		return( ans)
	}


	consensusAlignment <- function( m) {

		# given a matrix of a MSA alignment, select the best single call at every location
		# step 1: remove gap symbols from the ends of each sequence, but keep all internal gaps
		NR <- nrow(m)
		NC <- ncol(m)
		for ( i in 1:NR) {
			j <- 1
			while( m[i,j] == "-") { m[i,j] <- ""; j <- j + 1}
			j <- ncol(m)
			while( m[i,j] == "-") { m[i,j] <- ""; j <- j - 1}
		}
		# step 2: get the most common base call at each column.  Let's expand this ask to get some info about
		# variants seen at each location
		baseCallList <- vector( mode="list")
		baseCntList <- vector( mode="list")
		basePtr <- 0
		chBase <- apply( m, 2, function(x) {
			x <- x[ x != ""]
			if ( ! length(x)) x <- ""
			xt <- sort( table( x), decreasing=T)
			basePtr <<- basePtr + 1
			baseCallList[[basePtr]] <<- names(xt)
			baseCntList[[basePtr]] <<- as.integer(xt)
			return( names(xt[1]))
		})
		chDepth <- apply( m, 2, function(x) {
			x <- x[ x != ""]
			return( length(x))
		})
		# step 3: remove the ends that are dominated by gaps
		startFrom <- 1
		while (chBase[startFrom] %in% c( "", " ", "-")) startFrom <- startFrom + 1
		endAt <- NC
		while (chBase[endAt] %in% c( "", " ", "-")) endAt <- endAt - 1
		finalBase <- chBase[ startFrom : endAt]
		finalDepth <- chDepth[ startFrom : endAt]
		baseCallList <- baseCallList[ startFrom : endAt]
		baseCntList <- baseCntList[ startFrom : endAt]
		# lastly, remove any embedded gaps
		drops <- which( finalBase == "-")
		if ( length(drops)) {
			finalBase <- finalBase[ -drops]
			finalDepth <- finalDepth[ -drops]
			baseCallList <- baseCallList[ -drops]
			baseCntList <- baseCntList[ -drops]
		}
		outSeq <- paste( finalBase, collapse="")
		outDepth <- round( mean( finalDepth, na.rm=T), digits=2)
		# turn the diversity data into something we can see/manipulate.
		# put a floor on how much read depth is needed to trust a minor population
		diversity <- vector()
		for ( k in 1:length(baseCallList)) {
			bset <- baseCallList[[k]]
			bdepth <- baseCntList[[k]]
			myDepth <- finalDepth[k]
			# put a floor on how much read depth is needed to trust a minor population. 
			useDepth <- rep( max(myDepth,100), length(bdepth))
			useDepth[1] <- myDepth
			bpct <- round( bdepth * 100 / useDepth)
			keep <- which( bpct >= 2)
			diversity[k] <- paste( bset[keep], bpct[keep], sep=":", collapse="; ")
		}
		diverseM <- data.frame( "Position"=1:nchar(outSeq), "Base"=finalBase, "Depth"= finalDepth, "Diversity"=diversity, stringsAsFactors=F)
		return( list( "seq"=outSeq, "depth"=outDepth, "diversity"=diverseM))
	}


	MultiClusterMSA <- function() {
	
		# gather all the final AA sequences, and see which can/should be joined
		fset <- dir( RESULTS.PATH, pattern="Cluster\\.[0-9]+\\.Consensus.AA.fasta", full=T)
		if ( length(fset) < 1) {
			cat( "\nToo few AA sequences to do Multi-Cluster MSA.")
			return()
		}
		if ( length(fset) < 2 && is.null(EXPECT.AA.SEQ)) {
			cat( "\nToo few AA sequences to do Multi-Cluster MSA.")
			return()
		}
		faIN <- loadFasta( fset, verbose=F)
		seqs <- faIN$seq
		# simplify the descriptor line
		names(seqs) <- faIN$desc <- sub( paste(PREFIX,"_Cluster[._]?",sep=""), "Clust", faIN$desc)
		names(seqs) <- faIN$desc <- sub( "_Depth=[.0-9]+_", "_", names(seqs))
		# if we were given an expected AA sequence, then only keep those clusters that look like
		# what we are looking for
		if ( ! is.null(EXPECT.AA.SEQ)) {
			# this may be a AA seq or a FASTA filename
			if ( grepl( ".fasta$", EXPECT.AA.SEQ)) {
				tmpFA <- loadFasta( EXPECT.AA.SEQ, verbose=F)
				 EXPECT.AA.SEQ <<- tmpFA$seq[1]
				 names(EXPECT.AA.SEQ) <<- tmpFA$desc[1]
			} else {
				names(EXPECT.AA.SEQ) <<- "Expected.Seq"
			}
			data(BLOSUM62)
			pa <- pairwiseAlignment( seqs, EXPECT.AA.SEQ, type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
			# perfect match is about 5 per AA, so keep anything the is at least half good?
			keep <- which( (pa/nchar(EXPECT.AA.SEQ)) >= 2.5)
			seqs <- seqs[ keep]
			if ( length(seqs) < 1) {
				cat( "\nNo clusters resemble the expected AA sequence.")
				return()
			}
		}
		# when given an expected, add it to the mix
		if ( ! is.null(EXPECT.AA.SEQ)) {
			seqs <- c( EXPECT.AA.SEQ, seqs)
			names(seqs)[1] <- names(EXPECT.AA.SEQ)
		}
		faOUT <- as.Fasta( names(seqs), seqs)
		fastaFile <- file.path( RESULTS.PATH, paste( PREFIX, "All.Consensus.AA.fasta", sep="."))
		writeFasta( faOUT, fastaFile, line=100)
		if ( length(seqs) > 1) {
			alnFile <- file.path( RESULTS.PATH, paste( PREFIX, "All.Consensus.AA.aln", sep="."))
			aln <- mafft( fastaFile, alnFile, mode="genaffine")
			rownames(aln$alignment) <- names(seqs)  # put the fuller names back in place
			writeALN( aln, alnFile, line=100, max.id.width=45)
		} else {
			file.delete( alnFile)
		}
		return( length(seqs))
	}
	

	# all functions local to the closure are defined.
	# now call the top level function to do the clustering and return the environment
	ans <- do.all() 
	
	return( environment())
}
