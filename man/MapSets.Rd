\name{MapSets}
\alias{MapSets}
\alias{SeqMap}
\alias{GeneMap}
\alias{ExonMap}
\alias{RrnaMap}
\title{
Duffy Lab "MapSets" Overview
}
\description{
MapSets are the fundamental unit of genomic annotation for all the Duffy
Lab R packages.  They encapsulate all the information about genes, exons, 
and chromosomes for any organism, in a manner that is self-contained and
user configurable.
}
\details{
At the highest level, one MapSet defines all genomic annotation details for one
organism.  We use the terms 'species' and 'organism' interchangably, to mean
a unique strain of some genome with it's own unique annotation.  Each MapSet 
contains components of detailed information for the following layers of the
genome's annotation.
}

\section{Notation}{
\subsection{SpeciesID}{
Each species is given its own unique species identifier (SpeciesID), to 
distiquish it from any other species being used in any Duffy Lab R packages. 
Predefined SpeciesIDs include:  
\preformatted{
    'Pf3D7'   - for Plasmodium falciparum 3D7
    'Py17X'   - for Plasmodium yoelii 17XNL
    'PyYM'    - for Plasmodium yoelii YM
    'PbANKA'  - for Plasmodium bergheii ANKA
    'PvSAL1'  - for Plasmodium vivax SAL1
    'PchAS'   - for Plasmodium chaubadii AS
    'PkH'     - for Plasmodium knowlesi strain H
    'PyAABL'  - for Plasmodium yoelii (old non-chromosomal assembly)
    'Hs_grc'  - for Homo sapiens GRC37.p5
    'Mmu_grc' - for Mouse musculus 57BL/J6 strain, GRC38.p3
    'Agam'    - for Anopheles gambiae strain PEST
    'MT_H37'  - for M.tuberculosis H37RV_V2
}
Species are always referenced by their SpeciesID.  To see the set of all currently 
defined species identifiers, use \code{\link{getAllSpecies}}.  To get the current
SpeciesID, use \code{\link{getCurrentSpecies}}.  To change to a different species,
use \code{\link{setCurrentSpecies}}.  Even when doing mixed-organism analyses, 
only one species is actually 'current' at any point in time.
}

\subsection{SeqMap}{
The Sequence Map (SeqMap) is a table that identifies each chromosome or other supercontig in
the current species.  It has 2 columns (SEQ_ID, LENGTH).  The length is the length of
that contig in nucleotides, including any 'N's in the genomic FASTA sequence.
Every SequenceID (SeqID) must begin with the organism's SpeciesID, and then 
contain enough additional characters to uniquely identify each chromosome.  By forcing
the SeqIDs to contain the SpeciesID, we can guarantee sequence identifiers will
always be unique, even in mixed-organism analyses. To get the sequence map for
the current species, use \code{\link{getCurrentSeqMap}}.  

NOTE:  When building complementary FASTA genomic files, as for NextGen sequencing, 
be sure that the descriptor line for each FASTA chromosome matches the SeqID in 
the SeqMap.  This will guarantee that alignments will map exactly back to named 
sequences in the MapSet.  For example, the FASTA file for chromosome 1 of 
speciesID 'Pf3D7' begins with ">Pf3D7_01 ".
}

\subsection{GeneMap}{
The Gene Map (GeneMap) is a table that describes each gene level feature to be tracked for
the current species.  It has several required columns:
  \subsection{GENE_ID}{ the unique gene identifier (GeneID), that is unique over all defined
  species}
  \subsection{POSITION}{ the starting nucleotide position on its chromosome.}
  \subsection{END}{ the ending nucleotide position on its chromosome.  These are always
  given as forward strand positions, base 1 encoding, with POSITION < END, 
  regardless of the gene's coding strand.}
  \subsection{SEQ_ID}{ the sequence identifier (SeqID) for this gene.}
  \subsection{REAL_G}{ logical, is this gene real, or just an imaginary placeholder.  By default,
  all DuffyTools gene maps include extra 'non-genes' as placeholders for intergenic regions,
  to capture potential expression information from areas of the chromomsomes that 
  are not currently annotated as genes.}
  \subsection{STRAND}{ the coding strand '+' or '-' for this gene.  Genes without any strand
  information and non-genes use \code{NA}. }
  \subsection{NAME}{ the common name, nickname, or gene symbol.  This need not be unique.}
  \subsection{PRODUCT}{ the gene product descriptor, as supplementary info for text 
  files, plots, etc.}
  \subsection{N_EXON_BASES}{ the number of coding exonic bases, after splicing out introns.
  This is used as a proxy for the genes true size when normalizing read abundance 
  for RNA-seq.}

To get the gene map for the current species, use \code{\link{getCurrentGeneMap}}.  To get
the gene map for any other species, you must first change to that species 
via \code{\link{setCurrentSpecies}}

NOTE:  the requirement of GeneID uniqueness can be difficult to implement.  A good example
is the human genome, where multiple copies of the same gene is quite common.  For that case,
it took a combination of 4 terms to build a truly unique identifier.  Thus the human
GeneIDs are of the form:  '<gene_name>:<EntrezID>:<chromosome>:<nucleotide_location>'
This is essential for guaranteeing that every gene level feature on the genome is uniquely
identifiable.
}

\subsection{ExonMap}{
The Exon Map (ExonMap) is a table that describes each exon level feature to be tracked for
the current species.  It has several required columns:
  \subsection{GENE_ID}{ the gene identifier (GeneID), that this exon belongs to.}
  \subsection{POSITION}{ the starting nucleotide position on its chromosome.}
  \subsection{END}{ the ending nucleotide position on its chromosome.  These are always
  given as forward strand positions, base 1 encoding, with POSITION < END, 
  regardless of the gene's coding strand.}
  \subsection{SEQ_ID}{ the sequence identifier (SeqID) for this exon.}
  \subsection{STRAND}{ the coding strand '+' or '-' of this gene.}

To get the exon map for the current species, use \code{\link{getCurrentExonMap}}.  To get
the exon map for any other species, you must first change to that species 
via \code{\link{setCurrentSpecies}}
}

\subsection{RrnaMap}{
The RibosomalRNA Map (RrnaMap) is poorly name.  What it really represents is the gene
features that are targeted for special treatment, such as the Ribosomal Clearing step
of the RNA-seq pipeline, for removing unwanted ribosomal RNA and human globin from 
NextGen sequencing data.  It is a table that describes each gene level feature that may
warrant special attention in the current species.  It has several required columns:
  \subsection{GENE_ID}{ the gene identifier (GeneID).}
  \subsection{POSITION}{ the starting nucleotide position on its chromosome.}
  \subsection{END}{ the ending nucleotide position on its chromosome.  These are always
  given as forward strand positions, base 1 encoding, with POSITION < END, 
  regardless of the gene's coding strand.}
  \subsection{SEQ_ID}{ the sequence identifier (SeqID) for this exon.}
  \subsection{STRAND}{ the coding strand '+' or '-' for this gene.}
  \subsection{GROUP}{ a named grouping identifier that this gene is a member of.}
  \subsection{CLEAR}{ a logical, to flag this gene from removal and/or special
  treatment.  Currently used for building the pre-alignment index for the Ribosomal
  Clearing step of the RNA-seq pipeline.}

To get the ribosomal RNA map for the current species, use \code{\link{getCurrentRrnaMap}}.  To get
the RrnaMap for any other species, you must first change to that species 
via \code{\link{setCurrentSpecies}}
}

}
\section{Modifying MapSets}{
... to do...
}
\section{Defining New MapSets}{
GFFtoTargetMapSet.R
}
