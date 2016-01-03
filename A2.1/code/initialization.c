/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL
 * computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <metis.h>
#include <mpi.h>

#include "util_read_files.h"
#include "util_write_files.h"
#include "initialization.h"

int localize_points( int nintcf, int *local_global_index, int *g_elems, int **elems,
                     int g_points_count, int **g_points, int *points_count, int ***points );

int initialization( char *file_in, char *part_type, int *nintci, int *nintcf, int *nextci,
                    int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw,
                    double **bl, double **bh, double **bp, double **su, int *points_count,
                    int ***points, int **elems, double **var, double **cgup, double **oc,
                    double **cnorm, int **local_global_index, int **global_local_index,
                    int *neighbors_count, int **send_count, int ***send_list, int **recv_count,
                    int ***recv_list, idx_t **epart, idx_t **npart, idx_t *objval ) {
    /********** START INITIALIZATION **********/
    int my_rank, num_procs;
    int i = 0, j;

    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &num_procs );

    // Global variables for reading from file
    int g_nintci, g_nintcf, g_nextci, g_nextcf;
    int **g_lcc;
    double *g_bs, *g_be, *g_bn, *g_bw, *g_bh, *g_bl, *g_bp, *g_su;

    int g_points_count, **g_points, *g_elems;

    // read-in the input file
    int f_status = read_binary_geo( file_in, &g_nintci, &g_nintcf, &g_nextci, &g_nextcf, &g_lcc,
                                    &g_bs, &g_be, &g_bn, &g_bw, &g_bl, &g_bh, &g_bp, &g_su,
                                    &g_points_count, &g_points, &g_elems );

    if( f_status != 0 ) {
        return f_status;
    }

    // Partition computational domain
    idx_t ne = g_nintcf - g_nintci + 1;
    idx_t nn = g_points_count;
    idx_t nparts = num_procs;
    idx_t *eptr = ( idx_t * ) malloc( ( ne + 1 ) * sizeof( *eptr ) );
    for( i = 0; i <= ne; i++ ) {
        eptr[i] = 8 * i;
    }
    idx_t *eind = ( idx_t * ) malloc( ne * 8 * sizeof( *eind ) );
    for( i = 0; i < 8 * ne; i++ ) {
        eind[i] = g_elems[i];
    }
    *epart = ( idx_t * ) malloc( ne * sizeof( **epart ) );
    *npart = ( idx_t * ) malloc( g_points_count * sizeof( **npart ) );

    int status;
    if( nparts == 1 || !strcmp( part_type, "classical" ) ) {
        int length = ( ne - 1 ) / nparts + 1;
        for( i = 0; i < ne; i++ ) {
            ( *epart )[i] = i / length;
        }
        *objval = 0;
        status = METIS_OK;
    } else if( !strcmp( part_type, "dual" ) ) {
        idx_t ncommon = 3;
        status = METIS_PartMeshDual( &ne, &nn, eptr, eind, NULL, NULL, &ncommon,
                                     &nparts, NULL, NULL, objval, *epart, *npart );
    } else if( !strcmp( part_type, "nodal" ) ) {
        status = METIS_PartMeshNodal( &ne, &nn, eptr, eind, NULL, NULL,
                                      &nparts, NULL, NULL, objval, *epart, *npart );
    } else {
        fprintf( stderr, "Invalid partition type! Valid options: classical, dual, nodal.\n" );
        MPI_Abort( MPI_COMM_WORLD, my_rank );
    }

    if( status != METIS_OK ) {
        fprintf( stderr, "METIS failed with return value %d.\n", status );
        MPI_Abort( MPI_COMM_WORLD, my_rank );
    }
    free( eptr );
    free( eind );

#if 0
    // debugging
    if( my_rank == 0 ) {
        vtk_write_unstr_grid_header( "debug", "partition.vtk", g_nintci, g_nintcf, g_points_count,
                                     g_points, g_elems );
        int *part = ( int * ) malloc( ne * sizeof( int ) );
        for( i = 0; i < ne; i++ ) {
            part[i] = ( *epart )[i];
        }
        vtk_append_integer( "partition.vtk", "partition", 0, ne - 1, part );
        free( part );
    }
#endif

    // global_local_index
    if( ( *global_local_index = ( int * ) malloc( ( g_nextcf + 1 ) * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate global_local_index\n" );
        return -1;
    }

    for( i = 0; i <= g_nextcf; i++ ) {
        ( *global_local_index )[i] = -1;
    }

    *nintci = 0;
    *nintcf = *nintci - 1;
    for( i = 0; i < ne; i++ ) {
        if( ( *epart )[i] == my_rank ) {
            ( *global_local_index )[i] = ++*nintcf;
        }
    }

    *nextci = *nintcf + 1;
    *nextcf = *nextci - 1;
    for( i = 0; i < ne; i++ ) {
        if( ( *epart )[i] == my_rank ) {
            for( j = 0; j < 6; j++ ) {
                int idx = g_lcc[i][j];
                if( idx < g_nintci || idx > g_nintcf || ( *epart )[idx] != my_rank ) {
                    if( ( *global_local_index )[idx] == -1 ) {
                        ( *global_local_index )[idx] = ++*nextcf;
                    }
                }
            }
        }
    }

    // local_global_index
    if( ( *local_global_index = ( int * ) malloc( ( *nextcf + 1 ) * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate local_global_index\n" );
        return -1;
    }

    for( i = 0; i <= g_nextcf; i++ ) {
        int j = ( *global_local_index )[i];
        if( j != -1 ) {
            ( *local_global_index )[j] = i;
        }
    }

    // LCC
    if( ( *lcc = ( int ** ) malloc( ( *nintcf + 1 ) * sizeof( int * ) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate first dimension of LCC\n" );
        return -1;
    }

    for( i = 0; i < *nintcf + 1; i++ ) {
        if( ( ( *lcc )[i] = ( int * ) malloc( 6 * sizeof( int ) ) ) == NULL ) {
            fprintf( stderr, "malloc failed to allocate second dimension of LCC\n" );
            return -1;
        }
    }

    for( i = *nintci; i <= *nintcf; i++ ) {
        ( *lcc )[i][0] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][0]];
        ( *lcc )[i][1] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][1]];
        ( *lcc )[i][2] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][2]];
        ( *lcc )[i][3] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][3]];
        ( *lcc )[i][4] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][4]];
        ( *lcc )[i][5] = ( *global_local_index )[g_lcc[( *local_global_index )[i]][5]];
    }

    // allocate other arrays
    if( ( *bs = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BS) failed\n" );
        return -1;
    }

    if( ( *be = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BE) failed\n" );
        return -1;
    }

    if( ( *bn = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BN) failed\n" );
        return -1;
    }

    if( ( *bw = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BW) failed\n" );
        return -1;
    }

    if( ( *bl = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BL) failed\n" );
        return -1;
    }

    if( ( *bh = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BH) failed\n" );
        return -1;
    }

    if( ( *bp = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(BP) failed\n" );
        return -1;
    }

    if( ( *su = ( double * ) malloc( ( *nextcf + 1 ) * sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(SU) failed\n" );
        return -1;
    }

    for( i = *nintci; i <= *nintcf; i++ ) {
        ( *bs )[i] = g_bs[( *local_global_index )[i]];
        ( *be )[i] = g_be[( *local_global_index )[i]];
        ( *bn )[i] = g_bn[( *local_global_index )[i]];
        ( *bw )[i] = g_bw[( *local_global_index )[i]];
        ( *bl )[i] = g_bl[( *local_global_index )[i]];
        ( *bh )[i] = g_bh[( *local_global_index )[i]];
        ( *bp )[i] = g_bp[( *local_global_index )[i]];
        ( *su )[i] = g_su[( *local_global_index )[i]];
    }

    int l_status = localize_points( *nintcf, *local_global_index, g_elems, elems,
                                    g_points_count, g_points, points_count, points );
    if( l_status != 0 ) {
        return l_status;
    }

    // Free global variables
    for( i = g_nintci; i <= g_nintcf; i++ ) {
        free( g_lcc[i] );
    }
    free( g_lcc );
    free( g_bs );
    free( g_be );
    free( g_bn );
    free( g_bw );
    free( g_bh );
    free( g_bl );
    free( g_bp );
    free( g_su );
    free( g_elems );

    *var = ( double * ) calloc( ( *nextcf + 1 ), sizeof( double ) );
    *cgup = ( double * ) calloc( ( *nextcf + 1 ), sizeof( double ) );
    *cnorm = ( double * ) calloc( ( *nintcf + 1 ), sizeof( double ) );

    // initialize the arrays
    for( i = 0; i <= 10; i++ ) {
        ( *cnorm )[i] = 1.0;
    }

    for( i = ( *nintci ); i <= ( *nintcf ); i++ ) {
        ( *var )[i] = 0.0;
    }

    for( i = ( *nintci ); i <= ( *nintcf ); i++ ) {
        ( *cgup )[i] = 1.0 / ( ( *bp )[i] );
    }

    for( i = ( *nextci ); i <= ( *nextcf ); i++ ) {
        ( *var )[i] = 0.0;
        ( *cgup )[i] = 0.0;
        ( *bs )[i] = 0.0;
        ( *be )[i] = 0.0;
        ( *bn )[i] = 0.0;
        ( *bw )[i] = 0.0;
        ( *bh )[i] = 0.0;
        ( *bl )[i] = 0.0;
    }

    return 0;
}

int localize_points( int nintcf, int *local_global_index, int *g_elems, int **elems,
                     int g_points_count, int **g_points, int *points_count, int ***points ) {
    int *global2local;
    int i, j;

    // create global2local mapping for point indices
    if( ( global2local = ( int * ) malloc( g_points_count * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate global2local" );
        return -1;
    }

    for( i = 0; i < g_points_count; i++ ) {
        global2local[i] = -1;
    }

    *points_count = 0;
    for( i = 0; i <= nintcf; i++ ) {
        for( j = 0; j < 8; j++ ) {
            int idx = g_elems[8 * local_global_index[i] + j];
            if( global2local[idx] == -1 ) {
                global2local[idx] = ( *points_count )++;
            }
        }
    }

    // allocate and fill elems
    if( ( *elems = ( int * ) malloc( ( nintcf + 1 ) * 8 * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate elems" );
        return -1;
    }

    for( i = 0; i <= nintcf; i++ ) {
        for( j = 0; j < 8; j++ ) {
            int idx = g_elems[8 * local_global_index[i] + j];
            ( *elems )[8 * i + j] = global2local[idx];
        }
    }

    // allocate and fill points
    if( ( *points = ( int ** ) calloc( *points_count, sizeof( int * ) ) ) == NULL ) {
        fprintf( stderr, "malloc() POINTS 1st dim. failed\n" );
        return -1;
    }

    for( i = 0; i < g_points_count; i++ ) {
        int li = global2local[i];
        if( li == -1 ) {
            free( g_points[i] );
        } else {
            ( *points )[li] = g_points[i];
        }
    }

    free( g_points );
    free( global2local );

    return 0;
}
