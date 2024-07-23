# pdbTools.R -- read up a PDB file of atomic XYZ coordinates


`readPDB` <- function( file, verbose=TRUE) {

	txt <- readLines( file)

	# all we want at the ATOM and HETATM rows
	isATOM <- grep( "^ATOM", txt)
	isHET <- grep( "^HETATM", txt)
	keep <- sort( union( isATOM, isHET))
	txt <- txt[ keep]
	if ( ! length(txt)) return(NULL)

	# per the PDB specs, extract the fields by their character locations
	key <- gsub( " +", "", substr( txt, 1, 6))
	atNum <- as.integer( gsub( " +", "", substr( txt, 7, 11)))
	atNam <- gsub( " +", "", substr( txt, 13, 16))
	resNam <- gsub( " +", "", substr( txt, 18, 21))
	chain <- substr( txt, 22, 22)
	resNum <- as.integer( gsub( " +", "", substr( txt, 23, 26)))
	xCoord <- as.numeric( gsub( " +", "", substr( txt, 31, 38)))
	yCoord <- as.numeric( gsub( " +", "", substr( txt, 39, 46)))
	zCoord <- as.numeric( gsub( " +", "", substr( txt, 47, 54)))
	occ <- as.numeric( gsub( " +", "", substr( txt, 55, 60)))
	tempFac <- as.numeric( gsub( " +", "", substr( txt, 61, 66)))
	elNam <- gsub( " +", "", substr( txt, 77, 78))

	out <- data.frame( "Type"=key, "AtomNumber"=atNum, "AtomName"=atNam, "ResidueName"=resNam, "Chain"=chain,
			"ResidueNumber"=resNum, "X"=xCoord, "Y"=yCoord, "Z"=zCoord, "Occupany"=occ,
			"TempFactor"=tempFac, "Element"=elNam, stringsAsFactors=F)
	return( out)
}


`pdbAtomicDistance` <- function( pdb1, pdb2) {

	# given two sets of PDB atom records, calculate the inter-atomic distance in angstroms
	N1 <- nrow( pdb1)
	id1 <- paste( pdb1$AtomName, "_", pdb1$ResidueName, pdb1$ResidueNumber, sep="")
	N2 <- nrow( pdb2)
	id2 <- paste( pdb2$AtomName, "_", pdb2$ResidueName, pdb2$ResidueNumber, sep="")
	dm <- matrix( NA, N1, N2)
	rownames(dm) <- id1
	colnames(dm) <- id2
	for ( j in 1:N2) {
		dm[ ,j] <- xyzDist( pdb1$X, pdb1$Y, pdb1$Z, pdb2$X[j], pdb2$Y[j], pdb2$Z[j])
	}
	return(dm)
}


`pdbResidueDistance` <- function( pdb1, pdb2, mode=c("average","minimum")) {

	# given two sets of PDB atom records, calculate the inter-residue distance in angstroms
	# set up the storage
	resid1 <- paste( pdb1$ResidueName, pdb1$ResidueNumber, sep="")
	uniqRes1 <- unique( resid1)
	NR1 <- length( uniqRes1)
	resid2 <- paste( pdb2$ResidueName, pdb2$ResidueNumber, sep="")
	uniqRes2 <- unique( resid2)
	NR2 <- length( uniqRes2)
	resdm <- matrix( NA, NR1, NR2)
	rownames(resdm) <- uniqRes1
	colnames(resdm) <- uniqRes2
	
	# get the atom distances
	atomdm <- pdbAtomicDistance( pdb1, pdb2)

	mode <- match.arg( mode)
	for ( j in 1:NR2) {
		myCols <- which( resid2 == uniqRes2[j])
		for ( i in 1:NR1) {
			myRows <- which( resid1 == uniqRes1[i])
			dm <- atomdm[ myRows, myCols]
			if ( mode == "average") {
				avg <- mean( as.vector(dm), na.rm=T)
			} else if ( mode == "minimum") {
				avg <- min( as.vector(dm), na.rm=T)
			}
			resdm[i,j] <- avg
		}
	}
	return( resdm)
}


`xyzDist` <- function( x1, y1, z1, x2, y2, z2) {

	dx <- x2 - x1
	dy <- y2 - y1
	dz <- z2 - z1
	d <- sqrt( dx*dx + dy*dy + dz*dz)
	d
}
