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


#include "spectralfit.h"

/** GNU Scientific Library (GSL) **/
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>


int spectral_fit(img_t *images, table_t *bandlist, args_t *args){
int b, b_vector, b_highres, b_sharpened, b_spectralfit;
int p;
int nb = images[HIGHRES].meta.dim.band + images[SHARPENED].meta.dim.band;
gsl_vector *x, *y, *c;
gsl_bspline_workspace *work;
double chisq, est;
time_t TIME;

double min_wavelength;
double max_wavelength;

size_t control_points;

gsl_rng *rng;

int col_use, col_band, col_wavelength;

  time(&TIME);

  printf("Starting Spectral Fit\n")  ;


  if ((col_use  = find_table_col(bandlist, "use")) < 0){
    printf("there is no column 'use' in csv-file\n");
    exit(FAILURE);
  }
  if ((col_band = find_table_col(bandlist, "band")) < 0){
    printf("there is no column 'use' in csv-file\n");
    exit(FAILURE);
  }
  if ((col_wavelength = find_table_col(bandlist, "wavelength")) < 0){
    printf("there is no column 'wavelength' in csv-file\n");
    exit(FAILURE);
  }

  memcpy(&images[SPECTRALFIT].meta, &images[HIGHRES].meta, sizeof(meta_t));

  for (b=0, images[SPECTRALFIT].meta.dim.band=0; b<bandlist->nrow; b++){
    if ((int)bandlist->data[b][col_use] == 0) images[SPECTRALFIT].meta.dim.band++;
  }

  // nothing to do here
  if (images[SPECTRALFIT].meta.dim.band == 0) return(SUCCESS);

  for (b=0, min_wavelength=DBL_MAX, max_wavelength=DBL_MIN; b<bandlist->nrow; b++){
    if ((int)bandlist->data[b][col_use] == 1 || 
        (int)bandlist->data[b][col_use] == 2 ||
        (int)bandlist->data[b][col_use] == 0){
      if (bandlist->data[b][col_wavelength] < min_wavelength) min_wavelength = bandlist->data[b][col_wavelength];
      if (bandlist->data[b][col_wavelength] > max_wavelength) max_wavelength = bandlist->data[b][col_wavelength];
    }
  }

  alloc_2D((void***)&images[SPECTRALFIT].data, images[SPECTRALFIT].meta.dim.band, images[SPECTRALFIT].meta.dim.cell, sizeof(float));


  #pragma omp parallel private(b,b_highres,b_sharpened,b_spectralfit,b_vector,p,x,y,c,work,chisq,est,rng,control_points) shared(nb,images,bandlist,col_use, col_wavelength,args,min_wavelength,max_wavelength,gsl_rng_default) default(none)
  {

    /** initialize and allocate
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);

    // workspace
    work = gsl_bspline_alloc(args->order, args->nbreak);
    gsl_bspline_init_uniform(min_wavelength, max_wavelength, work);

    // number of control points
    control_points = gsl_bspline_ncontrol(work);


    // vector of nb observations
    x = gsl_vector_alloc(nb);
    y = gsl_vector_alloc(nb);

    // nv regression coefficients
    c = gsl_vector_calloc(control_points);


    #pragma omp for schedule(guided)  
    for (p=0; p<images[HIGHRES].meta.dim.cell; p++){

      if (images[NODATA].data[0][p] < 0){

        for (b=0, b_highres=0, b_sharpened=0, b_spectralfit=0; b<bandlist->nrow; b++){

          if ((int)bandlist->data[b][col_use] == 0){
          }
          if ((int)bandlist->data[b][col_use] == 1){
            images[HIGHRES].data[b_highres++][p] = images[HIGHRES].meta.nodata;
          } else if ((int)bandlist->data[b][col_use] == 2){
            images[SHARPENED].data[b_sharpened++][p] = images[SHARPENED].meta.nodata;
          } else if ((int)bandlist->data[b][col_use] == 0){
            images[SPECTRALFIT].data[b_spectralfit++][p] = images[SPECTRALFIT].meta.nodata;
          }
          
        }

        continue;
      } 

      for (b=0, b_vector=0, b_highres=0, b_sharpened=0; b<bandlist->nrow; b++){

        if ((int)bandlist->data[b][col_use] == 1){
          gsl_vector_set(y, b_vector, images[HIGHRES].data[b_highres++][p]);
        } else if ((int)bandlist->data[b][col_use] == 2){
          gsl_vector_set(y, b_vector, images[SHARPENED].data[b_sharpened++][p]);
        } else {
          continue;
        }

        gsl_vector_set(x, b_vector++, bandlist->data[b][col_wavelength]); 

      }

      gsl_bspline_lssolve(x, y, c, &chisq, work);

      for (b=0, b_highres=0, b_sharpened=0, b_spectralfit=0; b<bandlist->nrow; b++){
        gsl_bspline_calc(bandlist->data[b][col_wavelength], c, &est, work);
        if ((int)bandlist->data[b][col_use] == 1){
          images[HIGHRES].data[b_highres++][p] = (float)est;
        } else if ((int)bandlist->data[b][col_use] == 2){
          images[SHARPENED].data[b_sharpened++][p] = (float)est;
        } else if ((int)bandlist->data[b][col_use] == 0){
          images[SPECTRALFIT].data[b_spectralfit++][p] = (float)est;
        }
      }

    }


    /** clean
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
    gsl_rng_free(rng);
    gsl_vector_free(x);
    gsl_vector_free(y); 
    gsl_vector_free (c); 
    gsl_bspline_free(work);

  }


  proctime_print("Spectral fit", TIME);

  
  return SUCCESS;
}
