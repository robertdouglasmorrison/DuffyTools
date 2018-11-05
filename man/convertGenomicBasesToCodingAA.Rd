\name{convertGenomicBasesToCodingAA}
\alias{convertGenomicBasesToCodingAA}
\alias{convertGenomicDNAtoCodingDNA}
\alias{convertAApositionToGenomicDNAposition}
\title{
Conversions Between DNA and AA in Chromosome or Gene Local Context
}
\description{
Functions to try to convert DNA to/from AA given a local, possibly indel modified, context.
}
\usage{
convertGenomicBasesToCodingAminoAcids(seqID, position, end, strand = "+", 
		dnaQuery, genomeDNA, geneMap = NULL, exonMap = NULL)

convertGenomicDNAtoCodingDNA(geneID, genomeDNA = NULL)

convertAApositionToGenomicDNAposition(geneID, AAposition, AAlength)
}
\arguments{
  \item{seqID, geneID}{
a single character string giving the SeqID or GeneID that is the reference location in the genome
}
  \item{position, end}{
location in reference chromosomal units of the start and stop of the DNA sequence of interest
}
  \item{strand}{
the coding strand to use, for interpreting that region of DNA
}
  \item{dnaQuery}{
the DNA sequence for that region, as called by any non-reference means, such as by a SNP/Indel caller
}
  \item{genomeDNA}{
the DNA sequence for that region, from the genomic reference
}
  \item{geneMap, exonMap}{
the current species GeneMap and ExonMap, or the relavant subset thereof, needed to interpret the DNA
}
  \item{AAposition, AAlength}{
the location and extent of the coding amino acid sequence of interest
}
}
\details{
These functions try to map between DNA and AA given possibly modified genomic information, 
such as INDEL and SNP calling data.  They are not gauranteed to work in all possible contexts!  
Used mostly for adding DNA to AA notation to SNP plots, etc.
}
\value{
For \code{convertGenomicBasesToCodingAA}, a list of 3 elements.  Two are character vectors of 
the same length as the query DNA, having empty strings and amino acid calls that represent the
peptide/protein fragmant in the same physical layout as the given DNA.  The third is an integer
vector of amino acid ordinals giving the relative position in the final protein sequence. Insertions 
and/or deletions may cause elements of the query to be other than a single character.
  \item{genomic }{ the AA sequence of the reference genome in that region}
  \item{query }{ the AA sequence of the query sequence in that region}
  \item{aa.ordinal }{ the protein sequence locations of the query AA sequence in that region}

For \code{convertGenomicDNAtoCodingDNA}, a character string that is the cDNA sequence for the given gene, after removing
any introns and accounting for strand by possible reverse complementation.

For \code{convertAApositionToGenomicDNAposition}, a list of 2 elements:
  \item{SEQ_POSITION }{ the location in the genomic DNA that is the start of the first codon}
  \item{SEQ_END }{ the location in the genomic DNA that is the end of the last codon}
}
\author{
Bob Morrison
}
\note{
These functions are only as valid as the underlying genome and annotation.  They are not foolproof.
}
