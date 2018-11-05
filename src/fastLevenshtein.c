/* fastLevenshtein.c

	do the Levenshtein edit distance 
*/

#include <stdlib.h>
#include <stdio.h>


void fastLevenshtein( int *id1, int *n1, int *id2, int *n2, int *score) {

	//  Inputs --
      	//  id1, id2: 	vectors of integer characters,  1=A,2=B,etc. These represent the 2 sequences
	//  n1,n2:	length of each read

	//  Outputs --
	//  score:	number of edits ( substitutions, insertions, deletions)

	int i, j, l1, l2, size1, size2;
	int *F;
	int MATCH, MISMATCH, GAP_COST;
	int im1, im1offset, i1offset, ij;
	int thisGap1, thisGap2, thisF, bestF;

	l1 = *n1;
	l2 = *n2;
	size1 = l1 + 1;
	size2 = l2 + 1;

	F = (int *) malloc( (size_t) size1*size2*sizeof(int));
	for( i=0; i < (size1*size2); i++) F[i] = 0;

	// lets use string 1 going down the rows, string 2 across the columns
	#define true (1)
	#define false (0)

	// start the scores and trace at zero
	for(i = 0; i < size2; i++) {
		F[i] = 0;		// across the top all Js of string 2
	}
	for(i = 0; i < size1; i++) {
		F[i*size2] = 0; 	// down the left all Is of string 1
	}
	
	// fill the scores in
	MATCH = 0;
	MISMATCH = 1;
	GAP_COST = 1;

	for(i = 1; i < size1; i++) {
		im1 = i - 1;
		im1offset = im1 * size2;
		i1offset = i * size2;
		for(j = 1; j < size2; j++) {

			// get the 3 choices
			thisGap2 = F[ i1offset + j-1] + GAP_COST;
			thisGap1 = F[ im1offset + j] + GAP_COST;
			thisF = F[im1offset + (j-1)];
			if( id1[i-1] == id2[j-1]) {
				thisF += MATCH;
			} else {
				thisF += MISMATCH;
			}

			// whos best?
			ij = i1offset + j;
			if ( thisF < thisGap1 && thisF < thisGap2) {
				F[ij] = thisF;
			} else if ( thisGap1 < thisGap2) {
				F[ij] = thisGap1;
			} else {
				F[ij] = thisGap2;
			}
		}
	}

	// scores are all known, final edit distance is at the far corner
	ij = (size1 * size2) - 1;
	bestF = F[ ij];

	// ok, we now know where it started... so clean up...
	if (F) free( F);

	score[0] = bestF;
}
