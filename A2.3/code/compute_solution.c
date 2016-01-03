/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <metis.h>
#include <mpi.h>

static int exchange_direc1( double *direc1, int *local_global_index, idx_t *epart,
                            int neighbors_count, int *send_count, int **send_list,
                            int *recv_count, int **recv_list );

int compute_solution( const int max_iters, int nintci, int nintcf, int nextcf, int **lcc,
                      double *bp, double *bs, double *bw, double *bl, double *bn, double *be,
                      double *bh, double *cnorm, double *var, double *su, double *cgup,
                      double *residual_ratio, int *local_global_index, idx_t *epart,
                      int neighbors_count, int *send_count, int **send_list,
                      int *recv_count, int **recv_list ) {
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    /** the reference residual */
    double l_resref = 0.0, resref;

    /** array storing residuals */
    double *resvec = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );

    // initialize the reference residual
    for( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        l_resref = l_resref + resvec[nc] * resvec[nc];
    }
    MPI_Allreduce( &l_resref, &resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    resref = sqrt( resref );
    if( resref < 1.0e-15 ) {
        fprintf( stderr, "Residue sum less than 1.e-15 - %lf\n", resref );
        return 0;
    }

    /** the computation vectors */
    double *direc1 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double *direc2 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double *adxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *adxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *dxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double *dxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );

    while( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        exchange_direc1( direc1, local_global_index, epart, neighbors_count, send_count, send_list,
                         recv_count, recv_list );

        // compute new guess (approximation) for direc
        for( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        if( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }

            oc1 = occ / cnorm[1];
            for( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else if( nor1 == 2 ) {
            oc1 = 0;
            occ = 0;

            for( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }

            oc1 = occ / cnorm[1];
            oc2 = 0;
            occ = 0;
            for( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor2[nc];
            }

            oc2 = occ / cnorm[2];
            for( nc = nintci; nc <= nintcf; nc++ ) {
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
            }

            if2++;
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        omega = omega / cnorm[nor];
        double l_res_updated = 0.0, res_updated;
        for( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            l_res_updated = l_res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }
        MPI_Allreduce( &l_res_updated, &res_updated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        res_updated = sqrt( res_updated );
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if( *residual_ratio <= 1.0e-10 ) {
            break;
        }

        iter++;

        // prepare additional arrays for the next iteration step
        if( nor == nomax ) {
            nor = 1;
        } else {
            if( nor == 1 ) {
                for( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if( nor == 2 ) {
                    for( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

    free( direc1 );
    free( direc2 );
    free( adxor1 );
    free( adxor2 );
    free( dxor1 );
    free( dxor2 );
    free( resvec );

    return iter;
}

static int exchange_direc1( double *direc1, int *local_global_index, idx_t *epart,
                            int neighbors_count, int *send_count, int **send_list,
                            int *recv_count, int **recv_list ) {
    MPI_Request *send_request, *recv_request;
    double **send_buf, **recv_buf;
    int n, i;

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Start sending
    if( ( send_request = malloc( neighbors_count * sizeof( MPI_Request ) ) ) == NULL ) {
        return -1;
    }
    if( ( send_buf = malloc( neighbors_count * sizeof( double * ) ) ) == NULL ) {
        return -1;
    }
    for( n = 0; n < neighbors_count; n++ ) {
        if( ( send_buf[n] = malloc( send_count[n] * sizeof( double ) ) ) == NULL ) {
            return -1;
        }
        for( i = 0; i < send_count[n]; i++ ) {
            send_buf[n][i] = direc1[send_list[n][i]];
        }

        int dest = epart[local_global_index[recv_list[n][0]]];
        MPI_Isend( send_buf[n], send_count[n], MPI_DOUBLE, dest, 0, MPI_COMM_WORLD,
                   &send_request[n] );
    }

    // Start receiving
    if( ( recv_request = malloc( neighbors_count * sizeof( MPI_Request ) ) ) == NULL ) {
        return -1;
    }
    if( ( recv_buf = malloc( neighbors_count * sizeof( double * ) ) ) == NULL ) {
        return -1;
    }
    for( n = 0; n < neighbors_count; n++ ) {
        if( ( recv_buf[n] = malloc( recv_count[n] * sizeof( double ) ) ) == NULL ) {
            return -1;
        }

        int source = epart[local_global_index[recv_list[n][0]]];
        MPI_Irecv( recv_buf[n], recv_count[n], MPI_DOUBLE, source, 0, MPI_COMM_WORLD,
                   &recv_request[n] );
    }

    // Wait for data to be received
    for( n = 0; n < neighbors_count; n++ ) {
        MPI_Wait( &recv_request[n], MPI_STATUS_IGNORE );
        for( i = 0; i < recv_count[n]; i++ ) {
            direc1[recv_list[n][i]] = recv_buf[n][i];
        }
        free( recv_buf[n] );
    }
    free( recv_buf );
    free( recv_request );

    // Wait for data to be sent
    for( n = 0; n < neighbors_count; n++ ) {
        MPI_Wait( &send_request[n], MPI_STATUS_IGNORE );
        free( send_buf[n] );
    }
    free( send_buf );
    free( send_request );

    return 0;
}
