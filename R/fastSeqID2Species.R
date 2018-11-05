# fatSeqID2Species.R

# turn seqID into species, via the mapSets

`fastSeqID2Species` <-
function( seqid) {

	my_seqin <- seqid

	# set up storage to get filled
	my_speciesOut <- rep( NA, length=length( my_seqin))
	my_sptrOut <- rep( 0, length=length( my_seqin))

	# there are LOTS of duplicates that make this too slow!!
	# use factors to find all the unique <seqid>, do each once, and scatter the answer
	# back to all!!


	# local function for one SeqID at a time...
	`myFast_OneSeqFunc` <- function(x) {

		# given all the pointers to <seqid> tuples shared by one seqid
		thisSeq <- my_seqin[x[1]]
		thisSpecies <- setCurrentSpecies( speciesID=NULL, seqID=thisSeq)

		# if we can't find it...
		if ( is.na( thisSpecies)) {
			my_speciesOut[x] <<- NA
			my_sptrOut[x] <<- 0
			cat( "  No species mapSet has a SEQ_ID matching: ", thisSeq)
			return()
		}

		thisSeqPtr <- base::match( thisSeq, getCurrentSeqMap()$SEQ_ID, nomatch=0)
	
		my_speciesOut[ x] <<- thisSpecies
		my_sptrOut[ x] <<- thisSeqPtr
		return()
	}


	# factor and do it, one seqID at a time
	fseq <- factor( my_seqin)
	tapply( 1:length(my_seqin), INDEX=fseq, FUN=myFast_OneSeqFunc)

	# send back the 'filled in' answer
	out <- list( "SPECIES"=my_speciesOut, "SMAP_PTR"=my_sptrOut)
	return( out)
}
