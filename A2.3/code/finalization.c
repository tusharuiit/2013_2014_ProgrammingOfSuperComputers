/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <mpi.h>
#include "util_write_files.h"

void finalization( char *file_in, char *out_prefix, int total_iters, double residual_ratio,
                   int nintci, int nintcf, int points_count, int **points, int *elems, double *var,
                   double *cgup, double *su ) {
    int rank, status = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if( rank == 0 ) {
        char file_out[100];
        sprintf( file_out, "%s_summary.out", out_prefix );

        // TODO: Merge number of active cells
        // TODO: How to merge maxima and minima?
        status = store_simulation_stats( file_in, file_out, nintci, nintcf, var, total_iters,
                                         residual_ratio );
    }

    char file_out[100];
    sprintf( file_out, "%s_data_proc%d.vtk", out_prefix, rank + 1 );
    vtk_write_unstr_grid_header( file_in, file_out, nintci, nintcf, points_count, points, elems );
    vtk_append_double( file_out, "CGUP", nintci, nintcf, cgup );

    if( status != 0 ) {
        fprintf( stderr, "Error when trying to write to file %s\n", file_out );
    }
}

