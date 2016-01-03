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

static int localize_points( int nintcf, int *local_global_index, int *g_elems, int **elems,
                            int g_points_count, int **g_points, int *points_count, int ***points );

static int compute_neighbors_list( int nintci, int nintcf, int nextci, int nextcf, int **lcc,
                                   int *local_global_index, int *global_local_index, idx_t *epart,
                                   int *neighbors_count, int **send_count, int ***send_list,
                                   int **recv_count, int ***recv_list, idx_t nparts );

static void add_neighbor( int arr[], int *n, int x );
static void shrink( void **p, size_t size );

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

    status = compute_neighbors_list( *nintci, *nintcf, g_nextci, g_nextcf, *lcc,
                                     *local_global_index, *global_local_index, *epart,
                                     neighbors_count, send_count, send_list, recv_count, recv_list,
                                     nparts );
    if( status != 0 ) {
        return status;
    }

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

static int localize_points( int nintcf, int *local_global_index, int *g_elems, int **elems,
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

static int compute_neighbors_list( int nintci, int nintcf, int g_nextci, int g_nextcf, int **lcc,
                                   int *local_global_index, int *global_local_index, idx_t *epart,
                                   int *neighbors_count, int **send_count, int ***send_list,
                                   int **recv_count, int ***recv_list, idx_t nparts ) {
    int *rank;
    int i, j;

    if( ( *send_count = ( int * ) calloc( nparts, sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc(send_count) failed\n" );
        return -1;
    }

    // Count neighboring elements inside to send
    for( i = nintci; i <= nintcf; i++ ) {
        int neighbors[6], nn = 0;

        for( j = 0; j < 6; j++ ) {
            int idx = lcc[i][j];
            if( nintci <= idx && idx <= nintcf ) {
                // neighbor element also belongs to this node
                continue;
            }

            int g_idx = local_global_index[idx];
            if( g_nextci <= g_idx && g_idx <= g_nextcf ) {
                // neighbor element is external
                continue;
            }

            int n_idx = epart[g_idx];
            add_neighbor( neighbors, &nn, n_idx );
        }

        for( j = 0; j < nn; j++ ) {
            ( *send_count )[neighbors[j]]++;
        }
    }

    // Creating send_list
    if( ( *send_list = ( int ** ) calloc( nparts, sizeof( int * ) ) ) == NULL ) {
        fprintf( stderr, "malloc(send_list) failed\n" );
        return -1;
    }

    for( i = 0; i < nparts; i++ ) {
        if( ( *send_count )[i] != 0 ) {
            ( *send_list )[i] = ( int * ) malloc( ( *send_count )[i] * sizeof( int ) );
            if( ( *send_list )[i] == NULL ) {
                fprintf( stderr, "malloc(send_list[%d]) failed\n", i );
                return -1;
            }
        }
    }

    memset( *send_count, 0, nparts * sizeof( int ) );

    for( i = nintci; i <= nintcf; i++ ) {
        int neighbors[6], nn = 0;

        for( j = 0; j < 6; j++ ) {
            int idx = lcc[i][j];
            if( nintci <= idx && idx <= nintcf ) {
                // neighbor element also belongs to this node
                continue;
            }

            int g_idx = local_global_index[idx];
            if( g_nextci <= g_idx && g_idx <= g_nextcf ) {
                // neighbor element is external
                continue;
            }

            int n_idx = epart[g_idx];
            add_neighbor( neighbors, &nn, n_idx );
        }

        for( j = 0; j < nn; j++ ) {
            int n_idx = neighbors[j];
            ( *send_list )[n_idx][( *send_count )[n_idx]++] = i;
        }
    }

    // Count neighbors_count, create rank, compress send_count and send_list
    if( ( rank = ( int * ) malloc( nparts * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc(rank) failed\n" );
        return -1;
    }

    *neighbors_count = 0;
    for( i = 0; i < nparts; i++ ) {
        if( ( *send_count )[i] != 0 ) {
            int j = *neighbors_count;
            rank[j] = i;
            ( *send_list )[j] = ( *send_list )[i];
            ( *send_count )[j] = ( *send_count )[i];
            ( *neighbors_count )++;
        }
    }
    shrink( ( void ** )send_count, *neighbors_count * sizeof( int ) );
    shrink( ( void ** )send_list, *neighbors_count * sizeof( int * ) );

#if 0
    int my_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    printf( "neighbors_count(%d): %d\n", my_rank, *neighbors_count );
    for( i = 0; i < *neighbors_count; i++ ) {
        printf( "send_count(%d)[%d]: %d\n", my_rank, i, ( *send_count )[i] );
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        printf( "rank(%d)[%d]: %d\n", my_rank, i, rank[i] );
    }
#endif

    // Send and receive communication information
    MPI_Request *request;
    int **send_buf;

    if( ( request = malloc( *neighbors_count * sizeof( MPI_Request ) ) ) == NULL ) {
        fprintf( stderr, "malloc(request) failed\n" );
        return -1;
    }

    if( ( *recv_count = ( int * ) malloc( *neighbors_count * sizeof( int ) ) ) == NULL ) {
        fprintf( stderr, "malloc(recv_count) failed\n" );
        return -1;
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        MPI_Isend( &( *send_count )[i], 1, MPI_INT, rank[i], 0, MPI_COMM_WORLD, &request[i] );
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        MPI_Recv( &( *recv_count )[i], 1, MPI_INT, rank[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        MPI_Wait( &request[i], MPI_STATUS_IGNORE );
    }

    if( ( send_buf = malloc( *neighbors_count * sizeof( int * ) ) ) == NULL ) {
        fprintf( stderr, "malloc(send_buf) failed\n" );
        return -1;
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        send_buf[i] = ( int * ) malloc( ( *send_count )[i] * sizeof( int ) );
        if( send_buf[i] == NULL ) {
            fprintf( stderr, "malloc(send_buf[%d]) failed\n", i );
            return -1;
        }
        for( j = 0; j < ( *send_count )[i]; j++ ) {
            send_buf[i][j] = local_global_index[( *send_list )[i][j]];
        }
        MPI_Isend( send_buf[i], ( *send_count )[i], MPI_INT, rank[i], 0, MPI_COMM_WORLD,
                   &request[i] );
    }

    if( ( *recv_list = ( int ** ) malloc( *neighbors_count * sizeof( int * ) ) ) == NULL ) {
        fprintf( stderr, "malloc(recv_list) failed\n" );
        return -1;
    }
    for( i = 0; i < *neighbors_count; i++ ) {
        ( *recv_list )[i] = ( int * ) malloc( ( *recv_count )[i] * sizeof( int ) );
        if( ( *recv_list )[i] == NULL ) {
            fprintf( stderr, "malloc(recv_list[%d]) failed\n", i );
            return -1;
        }

        MPI_Recv( ( *recv_list )[i], ( *recv_count )[i], MPI_INT, rank[i], 0, MPI_COMM_WORLD,
                  MPI_STATUS_IGNORE );
        for( j = 0; j < ( *recv_count )[i]; j++ ) {
            ( *recv_list )[i][j] = global_local_index[( *recv_list )[i][j]];
        }
    }
    free( rank );

    for( i = 0; i < *neighbors_count; i++ ) {
        MPI_Wait( &request[i], MPI_STATUS_IGNORE );
        free( send_buf[i] );
    }
    free( request );
    free( send_buf );

    return 0;
}

static void add_neighbor( int arr[], int *n, int x ) {
    int i;
    for( i = 0; i < *n; i++ ) {
        if( arr[i] == x ) {
            break;
        }
    }
    if( i == *n ) {
        arr[*n] = x;
        ( *n )++;
    }
}

static void shrink( void **p, size_t size ) {
    void *r = realloc( *p, size );
    if( size == 0 || r != NULL ) {
        *p = r;
    }
}
