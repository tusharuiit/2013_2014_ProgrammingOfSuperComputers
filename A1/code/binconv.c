#include <stdio.h>
#include <stdlib.h>

#include "xread.h"

static int write_binary( char *fileName, int NINTCI, int NINTCF, int NEXTCI,
                         int NEXTCF, int **LCC, double *BS, double *BE, double *BN,
                         double *BW, double *BL, double *BH, double *BP, double *SU,
                         int *NBOARD );

int main( int argc, char *argv[] ) {
    if( argc < 3 ) {
        printf( "Usage: %s input_file output_file\n", argv[0] );
        return EXIT_FAILURE;
    }

    char *file_in = argv[1];
    char *file_out = argv[2];

    int status;

    /** internal cells start and end index*/
    int nintci, nintcf;
    /** external cells start and end index.
     * The external cells are only ghost cells. They are accessed only through internal cells*/
    int nextci, nextcf;
    /** link cell-to-cell array. Stores topology information*/
    int **lcc;
    /** red-black colouring of the cells*/
    int *nboard;

    /** boundary coefficients for each volume cell */
    double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;

    // read-in the input file
    status = read_formatted( file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                             &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard );

    if( status != 0 ) {
        printf( "failed to read data!\n" );
        return EXIT_FAILURE;
    }

    // write-out the output file
    status = write_binary( file_out, nintci, nintcf, nextci, nextcf, lcc,
                           bs, be, bn, bw, bl, bh, bp, su, nboard );

    if( status != 0 ) {
        printf( "failed to write data!\n" );
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static int write_binary( char *fileName, int NINTCI, int NINTCF, int NEXTCI,
                         int NEXTCF, int **LCC, double *BS, double *BE, double *BN,
                         double *BW, double *BL, double *BH, double *BP, double *SU,
                         int *NBOARD ) {
    int i;
    FILE *fp = fopen( fileName, "wb" );
    if( fp == NULL ) {
        printf( "Error opening file %s\n", fileName );
        return -1;
    }

    // 4 variables in total!!!
    fwrite( &NINTCI, sizeof( int ), 1, fp );
    fwrite( &NINTCF, sizeof( int ), 1, fp );
    fwrite( &NEXTCI, sizeof( int ), 1, fp );
    fwrite( &NEXTCF, sizeof( int ), 1, fp );

    for( i = 0; i < 6; i++ ) {
        fwrite( LCC[i] + NINTCI, sizeof( int ), NINTCF - NINTCI + 1, fp );
    }

    fwrite( BS + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BE + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BN + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BW + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BL + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BH + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( BP + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );
    fwrite( SU + NINTCI, sizeof( double ), NINTCF - NINTCI + 1, fp );

    fwrite( NBOARD + NINTCI, sizeof( int ), NINTCF - NINTCI + 1, fp );

    fclose( fp );
    return 0;
}
