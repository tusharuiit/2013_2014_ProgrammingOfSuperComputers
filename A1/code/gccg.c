#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <papi.h>

#include "xread.h"
#include "xwrite.h"
#include "vol2mesh.h"

#define SIZE( arr ) ( sizeof( arr ) / sizeof( arr[0] ) )

static void handle_error( int exit_code );

static void log_counters( FILE *out, const char *phase, long long *pTic, long long values[] );

static void write_vtk( char *prefix, char *outFileName, int startIdx, int endIdx,
                       int nodeCnt, int **points, int **elems, double *vector );


int main( int argc, char *argv[] ) {
    int Events[] = {
#ifdef CACHE_PROFILE
        PAPI_L2_TCM,
        PAPI_L3_TCM,
        PAPI_L2_TCA,
        PAPI_L3_TCA
#else
        PAPI_FP_OPS
#endif
    };
    long long values[SIZE( Events )];
    long long tic;

    if( argc != 4 ) {
        printf( "Usage: %s input_format input_file output_prefix\n", argv[0] );
        return EXIT_FAILURE;
    }

    char *input_format = argv[1];
    char *input_file = argv[2];
    char *output_prefix = argv[3];

    int status = 0;

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

    char pstats_filename[strlen( output_prefix ) + strlen( "pstats.dat" ) + 1];
    strcpy( pstats_filename, output_prefix );
    strcat( pstats_filename, "pstats.dat" );

    FILE *pstats = fopen( pstats_filename, "w" );
    if( pstats == NULL ) {
        printf( "Cannot open file for writing: %s\n", pstats_filename );
        return EXIT_FAILURE;
    }

    /* Start counting events */
    if( PAPI_start_counters( Events, SIZE( Events ) ) != PAPI_OK ) {
        handle_error( 1 );
    }

    /** start measuring wall clock time */
    tic = PAPI_get_real_usec();

    /* initialization  */
    // read-in the input file
    if( !strcmp( "bin", input_format ) ) {
        status = read_binary( input_file, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                              &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard );
    } else if( !strcmp( "text", input_format ) ) {
        status = read_formatted( input_file, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                 &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard );
    } else {
        printf( "valid input_format values: text, bin\n" );
        return EXIT_FAILURE;
    }

    if( status != 0 ) {
        printf( "failed to initialize data!\n" );
        return EXIT_FAILURE;
    }

    /* Print profile data for phase INPUT */
    log_counters( pstats, "INPUT", &tic, values );

    // allocate arrays used in gccg
    int nomax = 3;
    /** the reference residual*/
    double resref = 0.0;
    /** the ratio between the reference and the current residual*/
    double ratio;

    /** array storing residuals */
    double *resvec = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    /** the variation vector -> keeps the result in the end */
    double *var = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );

    /** the computation vectors */
    double *direc1 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double *direc2 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );

    /** additional vectors */
    double *cgup = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double *oc = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *cnorm = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *adxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *adxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *dxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *dxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );

    // initialize the reference residual
    for( int nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    resref = sqrt( resref );
    if( resref < 1.0e-15 ) {
        printf( "i/o - error: residue sum less than 1.e-15 - %lf\n", resref );
        return EXIT_FAILURE;
    }

    // initialize the arrays
    for( int nc = 0; nc <= 10; nc++ ) {
        oc[nc] = 0.0;
        cnorm[nc] = 1.0;
    }

    for( int nc = nintci; nc <= nintcf; nc++ ) {
        cgup[nc] = 0.0;
        var[nc] = 0.0;
    }

    for( int nc = nextci; nc <= nextcf; nc++ ) {
        var[nc] = 0.0;
        cgup[nc] = 0.0;
        direc1[nc] = 0.0;
        bs[nc] = 0.0;
        be[nc] = 0.0;
        bn[nc] = 0.0;
        bw[nc] = 0.0;
        bl[nc] = 0.0;
        bh[nc] = 0.0;
    }

    for( int nc = nintci; nc <= nintcf; nc++ ) {
        cgup[nc] = 1.0 / bp[nc];
    }

    int if1 = 0;
    int if2 = 0;
    int iter = 1;
    int nor = 1;
    int nor1 = nor - 1;
    /* finished initalization */

    /* start computation loop */
    while( iter < 10000 ) {
        /* start phase 1 */

        // update the old values of direc
        for( int nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess (approximation) for direc
        for( int nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
                         - bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
                         - bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
                         - bh[nc] * direc1[lcc[5][nc]];
        } /* end phase 1 */

        /*  start phase 2 */
        // execute normalization steps
        double oc1, oc2, occ;
        if( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;
            for( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }
            oc1 = occ / cnorm[1];
            for( int nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }
            if1++;

        } else if( nor1 == 2 ) {
            oc1 = 0;
            occ = 0;
            for( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }

            oc1 = occ / cnorm[1];
            oc2 = 0;
            occ = 0;
            for( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor2[nc] * direc2[nc];
            }

            oc2 = occ / cnorm[2];
            for( int nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
            }

            if2++;
        }

        cnorm[nor] = 0;
        double omega = 0;

        // compute the new residual
        for( int nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        omega = omega / cnorm[nor];

        double resnew = 0.0;
        for( int nc = nintci; nc <= nintcf; nc++ ) {
            var[nc] = var[nc] + omega * direc1[nc];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            resnew = resnew + resvec[nc] * resvec[nc];
        }
        resnew = sqrt( resnew );
        ratio = resnew / resref;

        // exit on no improvements of residual
        if( ratio <= 1.0e-10 ) {
            break;
        }

        iter++;

        // prepare additional arrays for the next iteration step
        if( nor == nomax ) {
            nor = 1;
        } else {
            if( nor == 1 ) {
                for( int nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }

            } else if( nor == 2 ) {
                for( int nc = nintci; nc <= nintcf; nc++ ) {
                    dxor2[nc] = direc1[nc];
                    adxor2[nc] = direc2[nc];
                }
            }
            nor++;
        }
        nor1 = nor - 1;
    }/* end phase 2 */

    /* finished computation loop */

    /* Print profile data for phase CALC */
    log_counters( pstats, "CALC", &tic, values );

    /* write output file  */
    int nodeCnt;
    int **points, **elems;

    if( vol2mesh( nintci, nintcf, lcc, &nodeCnt, &points, &elems ) != 0 ) {
        printf( "error during conversion from volume to mesh\n" );
    }

    write_vtk( output_prefix, "VAR.vtk", nintci, nintcf, nodeCnt, points, elems, var );
    write_vtk( output_prefix, "CGUP.vtk", nintci, nintcf, nodeCnt, points, elems, cgup );
    write_vtk( output_prefix, "SU.vtk", nintci, nintcf, nodeCnt, points, elems, su );

    /* Print profile data for phase OUTPUT */
    log_counters( pstats, "OUTPUT", &tic, values );

    /* Stop counting events */
    if( PAPI_stop_counters( values, SIZE( values ) ) != PAPI_OK ) {
        handle_error( 1 );
    }

    fclose( pstats );

#if 0
    /* Free all the dynamically allocated memory */
    free( direc2 );
    free( direc1 );
    free( dxor2 );
    free( dxor1 );
    free( adxor2 );
    free( adxor1 );
    free( cnorm );
    free( oc );
    free( var );
    free( cgup );
    free( resvec );
    free( su );
    free( bp );
    free( bh );
    free( bl );
    free( bw );
    free( bn );
    free( be );
    free( bs );
#endif

    printf( "Simulation completed successfully!\n" );
    return EXIT_SUCCESS;
}

static void handle_error( int exit_code ) {
    puts( "PAPI fails miserably!" );
    exit( exit_code );
}

/** Log performance counters into "pstats.dat" */
static void log_counters( FILE *out, const char *phase, long long *pTic, long long values[] ) {
    long long toc = PAPI_get_real_usec();

#ifdef CACHE_PROFILE
    if( PAPI_read_counters( values, 4 ) != PAPI_OK ) {
        handle_error( 1 );
    }

    fprintf( out, "%s PAPI_L2_TCM %lld\n", phase, values[0] );
    fprintf( out, "%s PAPI_L2_TCA %lld\n", phase, values[2] );
    fprintf( out, "%s L2MissRate %.4lf%\n", phase, ( double )values[0] / ( double )values[2] );
    fprintf( out, "%s PAPI_L3_TCM %lld\n", phase, values[1] );
    fprintf( out, "%s PAPI_L3_TCA %lld\n", phase, values[3] );
    fprintf( out, "%s L3MissRate %.4lf%\n", phase, ( double )values[1] / ( double )values[3] );
#else
    if( PAPI_read_counters( values, 1 ) != PAPI_OK ) {
        handle_error( 1 );
    }

    fprintf( out, "%s PAPI_FP_OPS %lld\n", phase, values[0] );
#endif

    fprintf( out, "%s RealTime %.4lfs\n", phase, ( toc - *pTic ) * 1e-6 );
    *pTic = PAPI_get_real_usec();
}

/** Utility function for writing .vtk files with prefix */
static void write_vtk( char *prefix, char *name, int start_idx, int end_idx,
                       int node_cnt, int **points, int **elems, double *vector ) {
    char file_name[strlen( prefix ) + strlen( name ) + 1];
    strcpy( file_name, prefix );
    strcat( file_name, name );

    if( write_result_vtk( file_name, start_idx, end_idx, node_cnt, points, elems, vector ) != 0 ) {
        printf( "error when trying to write to file %s\n", file_name );
    }
}
