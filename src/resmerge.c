/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This file is part of FORCE - Framework for Operational Radiometric 
Correction for Environmental monitoring.

Copyright (C) 2013-2022 David Frantz

FORCE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FORCE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FORCE.  If not, see <http://www.gnu.org/licenses/>.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions to enhance spatial resolution
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "resmerge.h"

/** GNU Scientific Library (GSL) **/
#include <gsl/gsl_multifit.h>          // multi-parameter fitting



int resolution_merge(img_t *images, args_t *args){
int b = 0;
int i, j, p, ii, jj, ni, nj, np;
int w, nw, k, nv = images[PCA].meta.dim.band + 1;
bool nodata;
gsl_matrix *X, **cov;
gsl_vector *x, **y, **c;
gsl_multifit_linear_workspace **work;
double chisq, est, err;
time_t TIME;

  
  time(&TIME);

  printf("Starting Resolution Merge\n")  ;



  // kernel size
  w = 2 * args->radius + 1;
  nw = w * w;

//gsl_set_error_handler_off();
  #pragma omp parallel private(k,b,j,p,ii,jj,ni,nj,np,X,x,y,c,cov,work,chisq,est,err,nodata) shared(w,nw,nv,images,args) default(none)
  {

    /** initialize and allocate
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    // nw-by-nv predictor variables; kernel + central pixel
    X = gsl_matrix_calloc(nw, nv);
    x = gsl_vector_calloc(nv);

    // set first column of X to 1 -> intercept c0
    for (k=0; k<nw; k++) gsl_matrix_set(X, k, 0, 1.0);
    gsl_vector_set(x, 0, 1.0);

    // vector of nw observations
    alloc((void**)&y, images[LOWRES].meta.dim.band, sizeof(gsl_vector*));
    for (b=0; b<images[LOWRES].meta.dim.band; b++) y[b] = gsl_vector_calloc(nw);

    // nv regression coefficients
    alloc((void**)&c, images[LOWRES].meta.dim.band, sizeof(gsl_vector*));
    for (b=0; b<images[LOWRES].meta.dim.band; b++) c[b] = gsl_vector_calloc(nv);

    // nv-by-nv covariance matrix
    alloc((void**)&cov, images[LOWRES].meta.dim.band, sizeof(gsl_matrix*));
    for (b=0; b<images[LOWRES].meta.dim.band; b++) cov[b] = gsl_matrix_calloc(nv, nv);

    // workspace
    alloc((void**)&work, images[LOWRES].meta.dim.band, sizeof(gsl_multifit_linear_workspace*));
    for (b=0; b<images[LOWRES].meta.dim.band; b++) work[b] = gsl_multifit_linear_alloc(nw, nv);

    // sharpened dataset
    alloc_2D((void***)&images[SHARPENED].data, images[LOWRES].meta.dim.band, images[LOWRES].meta.dim.cell, sizeof(float));


    /** do regression for every valid pixel, and for each 20m band
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    #pragma omp for schedule(guided)  
    for (i=0; i<images[PCA].meta.dim.row; i++){
    for (j=0; j<images[PCA].meta.dim.col; j++){

      p = i*images[PCA].meta.dim.col+j;

      if (images[NODATA].data[0][p] < 0){

        for (b=0; b<images[LOWRES].meta.dim.band; b++){
          images[SHARPENED].data[b][p] = images[LOWRES].meta.nodata;
        }

        continue;

      }

      // add central pixel
      for (b=0; b<images[PCA].meta.dim.band; b++) gsl_vector_set(x, b+1, images[PCA].data[b][p]);
      
      k = 0;

      // add neighboring pixels
      for (ii=-args->radius; ii<=args->radius; ii++){
      for (jj=-args->radius; jj<=args->radius; jj++){

        
        if (ii < 0) ni = i-ii*ii; else ni = i+ii*ii;
        if (jj < 0) nj = j-jj*jj; else nj = j+jj*jj;

        if (ni < 0 || ni >= images[PCA].meta.dim.row || nj < 0 || nj >= images[PCA].meta.dim.col) continue;
        np = ni*images[PCA].meta.dim.col+nj;

        if (images[NODATA].data[0][np] < 0) continue;

        for (b=0, nodata=0; b<images[LOWRES].meta.dim.band; b++){

          if (fequal(images[LOWRES].data[0][np], images[LOWRES].meta.nodata)){
            nodata = true;
            break;
          }

          gsl_vector_set(y[b], k, images[LOWRES].data[b][np]);

        }

        if (!nodata){
          for (b=0; b<images[PCA].meta.dim.band; b++) gsl_matrix_set(X, k, b+1, images[PCA].data[b][np]);
          k++;
        }

      }
      }

      if (k < nw/2){

        for (b=0; b<images[LOWRES].meta.dim.band; b++){
          images[SHARPENED].data[b][p] = images[LOWRES].meta.nodata;
        }

        images[NODATA].data[0][p] = -10000.0;

        continue;
        
      }


      // append zeros, if less than nw neighboring pixels were added
      while (k < nw){
        for (b=0; b<images[PCA].meta.dim.band; b++) gsl_matrix_set(X, k, b+1, 0.0);
        for (b=0; b<images[LOWRES].meta.dim.band; b++) gsl_vector_set(y[b], k, 0.0);
        k++;
      }

      // solve model, and predict central pixel
      for (b=0; b<images[LOWRES].meta.dim.band; b++){

        gsl_multifit_linear(X, y[b], c[b], cov[b], &chisq, work[b]);
        gsl_multifit_linear_est(x, c[b], cov[b], &est, &err);
        images[SHARPENED].data[b][p] = est;

      }

    }
    }


    /** clean
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
    gsl_matrix_free (X); gsl_vector_free (x);
    for (b=0; b<images[LOWRES].meta.dim.band; b++) gsl_vector_free(y[b]); 
    for (b=0; b<images[LOWRES].meta.dim.band; b++) gsl_vector_free (c[b]); 
    for (b=0; b<images[LOWRES].meta.dim.band; b++) gsl_matrix_free (cov[b]); 
    for (b=0; b<images[LOWRES].meta.dim.band; b++) gsl_multifit_linear_free(work[b]); 
    free((void*)y);      free((void*)c);
    free((void*)cov);    free((void*)work);

  }

//  gsl_set_error_handler(NULL);


  memcpy(&images[SHARPENED].meta, &images[LOWRES].meta, sizeof(meta_t));

  proctime_print("Resolution merge", TIME);

  
  return SUCCESS;
}
