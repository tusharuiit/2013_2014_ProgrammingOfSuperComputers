/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>

#include "util_read_files.h"
#include "util_write_files.h"

int test_distribution( char *file_in, char *file_vtk_out, int *local_global_index,
                       int local_num_elems, double *cgup_local ) {
    // Global variables for reading from file
    int nintci, nintcf, nextci, nextcf;
    int **lcc;
    double *bs, *be, *bn, *bw, *bh, *bl, *bp, *su;

    int points_count, **points, *elems;
    double *distr;

    int i;

    // read-in the input file
    int f_status = read_binary_geo( file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                    &bs, &be, &bn, &bw, &bl, &bh, &bp, &su,
                                    &points_count, &points, &elems );

    if( f_status != 0 ) {
        return f_status;
    }

    // Free unnecessary global data
    for( i = nintci; i <= nintcf; i++ ) {
        free( lcc[i] );
    }
    free( lcc );
    free( bs );
    free( be );
    free( bn );
    free( bw );
    free( bh );
    free( bl );
    free( bp );
    free( su );

    // Allocate and fill distr
    if( ( distr = ( double * ) calloc( nintcf - nintci + 1, sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(distr) failed\n" );
        return -1;
    }

    for( i = 0; i < local_num_elems; i++ ) {
        distr[local_global_index[i]] = cgup_local[i];
    }

    // Write vtk file
    vtk_write_unstr_grid_header( "test_distribution", file_vtk_out, nintci, nintcf, points_count,
                                 points, elems );
    vtk_append_double( file_vtk_out, "CGUP", nintci, nintcf, distr );

    // Free the rest of global data
    for( i = 0; i < points_count; i++ ) {
        free( points[i] );
    }
    free( points );
    free( elems );
    free( distr );

    return 0;
}

int test_communication( char *file_in, char *file_vtk_out, int *local_global_index,
                        int local_num_elems, int neighbors_count, int *send_count, int **send_list,
                        int *recv_count, int **recv_list ) {
    // Global variables for reading from file
    int nintci, nintcf, nextci, nextcf;
    int **lcc;
    double *bs, *be, *bn, *bw, *bh, *bl, *bp, *su;

    int points_count, **points, *elems;
    double *commlist;

    int i, j;

    // read-in the input file
    int f_status = read_binary_geo( file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                    &bs, &be, &bn, &bw, &bl, &bh, &bp, &su,
                                    &points_count, &points, &elems );

    if( f_status != 0 ) {
        return f_status;
    }

    // Free unnecessary global data
    for( i = nintci; i <= nintcf; i++ ) {
        free( lcc[i] );
    }
    free( lcc );
    free( bs );
    free( be );
    free( bn );
    free( bw );
    free( bh );
    free( bl );
    free( bp );
    free( su );

    // Allocate and fill commlist
    if( ( commlist = ( double * ) calloc( nintcf - nintci + 1, sizeof( double ) ) ) == NULL ) {
        fprintf( stderr, "malloc(commlist) failed\n" );
        return -1;
    }

    for( i = 0; i < local_num_elems; i++ ) {
        // internal cells
        commlist[local_global_index[i]] = 15;
    }

    for( i = 0; i < neighbors_count; i++ ) {
        for( j = 0; j < send_count[i]; j++ ) {
            // ghost cell to be sent
            commlist[local_global_index[send_list[i][j]]] = 10;
        }
    }

    for( i = 0; i < neighbors_count; i++ ) {
        for( j = 0; j < recv_count[i]; j++ ) {
            // ghost cell to be received
            commlist[local_global_index[recv_list[i][j]]] = 5;
        }
    }

    // Write vtk file
    vtk_write_unstr_grid_header( "test_communication", file_vtk_out, nintci, nintcf, points_count,
                                 points, elems );
    vtk_append_double( file_vtk_out, "commlist", nintci, nintcf, commlist );

    // Free the rest of global data
    for( i = 0; i < points_count; i++ ) {
        free( points[i] );
    }
    free( points );
    free( elems );
    free( commlist );

    return 0;
}

