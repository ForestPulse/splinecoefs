#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include <omp.h>

// include stuff
#include "dtype.h"
#include "alloc.h"
#include "utils.h"
#include "dir.h"
#include "image_io.h"
#include "date.h"
#include "quality.h"
#include "table.h"
#include "usage.h"


/** GNU Scientific Library (GSL) **/
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>


enum { BOA, QAI, THT, N_PROD };


// example parameters
// a B-spline of order k is a polynomial of degree k-1, 4 = cubic
// k = 4
// number of control points / coefficients to estimate
// n = 22
// smoothing parameter (adjust as needed)
// l = 1000

// -i input-table = csv file with input images
// 3 columns: BOA, QAI, THT
// each row corresponds to one date, and the 3 columns must have the same date
// dates should be ordered from earliest to latest
// file basenames must start with date in format YYYYMMDD, e.g., 20200731_LEVEL2_SEN2A_BOA.tif




int determine_annual_weights(int order, int n_control, double max_weight, int target_year, date_t **dates, int n_dates, int n_years, image_t **input, int band_nir, int band_red, image_t *mask, image_t *weight){


  #pragma omp parallel shared(order, n_control, max_weight, target_year, dates, n_dates, n_years, input, band_nir, band_red, mask, weight, stderr) default(none)
  {

    double min_x = 0.0;
    double max_x = 365.0;
    int n_doys = 365;
    
    enum { START, END, RANGE };


    // workspace
    gsl_bspline_workspace *work;
    work = gsl_bspline_alloc_ncontrol(order, n_control);
    gsl_bspline_init_uniform(min_x, max_x, work);
    
    // coefficients
    gsl_vector *c= gsl_vector_calloc(n_control);
    
    // data vectors
    int vector_size = n_dates + n_control; // max possible number of data points (all dates plus possible padded points for each year)

    
    gsl_vector **x = NULL;
    gsl_vector **y = NULL;
    gsl_vector **w = NULL;
    int         *n = NULL;
    alloc((void**)&x, n_years, sizeof(gsl_vector*));
    alloc((void**)&y, n_years, sizeof(gsl_vector*));
    alloc((void**)&w, n_years, sizeof(gsl_vector*));
    alloc((void**)&n, n_years, sizeof(int));
    for (int year=0; year<n_years; year++){
      x[year] = gsl_vector_alloc(vector_size);
      y[year] = gsl_vector_alloc(vector_size);
      w[year] = gsl_vector_alloc(vector_size);
    }

    double **predicted = NULL;
    alloc_2D((void***)&predicted, n_years, n_doys, sizeof(double));

    double *mean = NULL;
    alloc((void**)&mean, n_years, sizeof(double));

    #pragma omp for
    for (int p=0; p<mask->nc; p++){

      if (mask->data[0][p] == mask->nodata || mask->data[0][p] == 0) continue;

      // reset variables and counter for each year
      for (int year=0; year<n_years; year++){
        weight->data[year][p] = 1;
        for (int d=0; d<vector_size; d++){
          x[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
          y[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
          w[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
        }
        n[year] = 0;
      }

      for (int d=0; d<n_dates; d++){

        if (!use_this_pixel(input[d][QAI].data[0][p])) continue;

        int year = target_year - dates[d][BOA].year;

        x[year]->data[n[year]] = input[d][THT].data[0][p];
        double denom = input[d][BOA].data[band_nir][p] + input[d][BOA].data[band_red][p];
        if (denom == 0) denom = 1;
        double ndvi = (input[d][BOA].data[band_nir][p] - input[d][BOA].data[band_red][p]) / denom;
        y[year]->data[n[year]] = ndvi;
        if (ndvi < 0){
          w[year]->data[n[year]] = 0.01;
        } else {
          w[year]->data[n[year]] = sqrt(ndvi);
        }

        n[year]++;

      }

            
      for (int year=0; year<n_years; year++){
      
        if (n[year] < 1){
          fprintf(stderr, "At least one datapoiunt is required (pixel %d, year %d): %d points.\n", p, year, n[year]);
          if (year == 0) break; else continue; // need at least as many data points as control points for fitting
        } 

        int added_points = 0;

        if (x[year]->data[0] > 365/n_control){
          
          int seg_first = (int)(x[year]->data[0] / (365/n_control));

          for (int s=0; s<seg_first; s++){
            x[year]->data[n[year]+added_points] = s * 365 / n_control + (365 / n_control / 2); // add a point at regular interval for gap
            y[year]->data[n[year]+added_points] = y[year]->data[0]; // add a point at regular interval for gap
            w[year]->data[n[year]+added_points] = w[year]->data[0] * 0.1; // add a point at regular interval for gap
            added_points++;
          }

        }

        
        if (365 - x[year]->data[n[year]-1] > 365/n_control){

          int seg_last = (int)((365 - x[year]->data[n[year]-1]) / (365/n_control));
        
          for (int s=0; s<seg_last; s++){
            x[year]->data[n[year]+added_points] = 365 - s * 365 / n_control - 365 / n_control / 2; // add a point at regular interval for gap
            y[year]->data[n[year]+added_points] = y[year]->data[n[year]-1]; // add a point at regular interval for gap
            w[year]->data[n[year]+added_points] = w[year]->data[n[year]-1] * 0.1; // add a point at regular interval for gap
            added_points++;
          }
        }

        for (int d=1; d<n[year]; d++){

          if (x[year]->data[d] - x[year]->data[d-1] > 365/n_control){

            int seg_prev = (int)(x[year]->data[d-1] / (365/n_control));
            int seg_now  = (int)(x[year]->data[d]   / (365/n_control));

            for (int s=seg_prev+1; s<seg_now; s++){
              x[year]->data[n[year]+added_points] = s * 365 / n_control + (365 / n_control / 2); // add a point at regular interval for gap
              y[year]->data[n[year]+added_points] = y[year]->data[d-1]; // add a point at regular interval for gap
              w[year]->data[n[year]+added_points] = w[year]->data[d-1] * 0.1; // add a point at regular interval for gap
              added_points++;
            }

          }
        }
        
        n[year] += added_points;

        double chisq, est;
        gsl_bspline_wlssolve(x[year], y[year], w[year], c, &chisq, work);

        mean[year] = 0.0;
        for (int doy=0; doy<n_doys; doy++){
          gsl_bspline_calc(doy, c, &est, work);
          predicted[year][doy] = est;
          mean[year] += est;
        }

        mean[year] /= n_doys;

      }

      if (n[0] < n_control){
        fprintf(stderr, "No data for target year %d for pixel %d.\n", target_year, p);
        continue; // need data for target year
      } 

      double dist_max = 0.0;
      for (int doy=0; doy<n_doys; doy++){
        dist_max += (predicted[0][doy] - mean[0]) * 
                    (predicted[0][doy] - mean[0]);
      }
      dist_max = sqrt(dist_max/n_doys);

      if (dist_max == 0.0){
        fprintf(stderr, "No variation in predicted values for target year %d for pixel %d.\n", target_year, p);
        continue; // avoid division by zero, means all values for target year are the same, so no point in weighting
      }

      weight->data[0][p] = 10000; // maximum weight for target year, scaled to 0-10000 for storage as short

      for (int year=1; year<n_years; year++){

        if (n[year] < n_control){
          fprintf(stderr, "No data for year %d for pixel %d.\n", target_year - year, p);
          continue; // skip years with no data
        } 
        
        double dist = 0.0;
        for (int doy=0; doy<n_doys; doy++){
          dist += (predicted[year][doy] - predicted[0][doy]) * 
                  (predicted[year][doy] - predicted[0][doy]);
        }
        dist = sqrt(dist/n_doys);

        //double w = max_weight * (1.0 - dist / dist_max);
        //double w = 0.5;
        //double w = 1.0/((year+1.0)*(year+1.0));
        double w = 0.0;
        if (year == 1){
          w = 0.25;
        }
        if (year == 2){
          w = 0.0;
        }

        if (w > 0){
          weight->data[year][p] = w * 10000; // scale weight to 0-10000 for storage as short
        }

      }

    } // end pixel-loop
      

    for (int year=0; year<n_years; year++){
      gsl_vector_free(x[year]);
      gsl_vector_free(y[year]);
      gsl_vector_free(w[year]);
    }
    free((void*)x);
    free((void*)y);
    free((void*)w);
    free((void*)n);
    gsl_vector_free(c);
    gsl_bspline_free(work);
    free_2D((void**)predicted, n_years);
    free((void*)mean);

  } // end parallel region

  return SUCCESS;
}

//####################################
// Hilffunktion zum setzten der Basisgewichte pro Jahr
static void set_base_weights(
    int p,
    int target_year,
    date_t **dates,
    int n_dates,
    image_t **input,
    image_t *weight)
{
    for(int d=0; d<n_dates; d++){

        if(!use_this_pixel(input[d][QAI].data[0][p]))
            continue;

        int year_diff = target_year - dates[d][BOA].year;

        double base_w = 0.0;

        if(year_diff==0) base_w = 1.0;
        else if(year_diff==1) base_w = 0.3;
        else if(year_diff==2) base_w = 0.15;
        else continue;

        weight->data[d][p] = base_w * 10000;
    }
}

int determine_point_weights(
    int order,
    int n_control,
    int target_year,
    date_t **dates,
    int n_dates,
    image_t **input,
    int band_nir,
    int band_red,
    image_t *mask,
    image_t *weight){

  #pragma omp parallel shared(order,n_control,target_year,dates,n_dates,input,band_nir,band_red,mask,weight) default(none)
  {
    int n_pad = 3;
    //Defines the range of the spline
    double min_x = -n_pad;
    double max_x = 365 + n_pad;

    
    gsl_vector *c = gsl_vector_calloc(n_control);
    gsl_vector *tmpB = gsl_vector_calloc(n_control);

    // defines Workspace and initializes the B-spline basis functions
    gsl_bspline_workspace *work;
    work = gsl_bspline_alloc_ncontrol(order,n_control);
    gsl_bspline_init_uniform(min_x,max_x,work);


    int n_x_max = n_dates + 2*n_pad;

    double *x=NULL;
    double *y=NULL;
    double *w=NULL;

    alloc((void**)&x,n_x_max,sizeof(double));
    alloc((void**)&y,n_x_max,sizeof(double));
    alloc((void**)&w,n_x_max,sizeof(double));

    int *obs_index=NULL;
    alloc((void**)&obs_index,n_x_max,sizeof(int));

    #pragma omp for
    for(int p=0;p<mask->nc;p++){
        if(mask->data[0][p]==mask->nodata || mask->data[0][p]==0) continue;

        for(int d=0; d<n_dates; d++)      // initialize weights
          weight->data[d][p] = 0;

        int n_x=0;
        for(int d=0; d<n_dates; d++){

            if(!use_this_pixel(input[d][QAI].data[0][p])) continue;

            int year_diff = target_year - dates[d][BOA].year;
            double base_w = 0.0;
            if(year_diff==0) base_w = 1.0;
            else if(year_diff==1) base_w = 0.3;
            else if(year_diff==2) base_w = 0.15;
            else continue;

            x[n_x] = input[d][THT].data[0][p];

            double denom = input[d][BOA].data[band_nir][p] +
                          input[d][BOA].data[band_red][p];

            if(denom==0) denom=1;

            double ndvi = (input[d][BOA].data[band_nir][p] -
                          input[d][BOA].data[band_red][p]) / denom;

            y[n_x] = ndvi;
            w[n_x] = base_w;

            obs_index[n_x] = d; 

            n_x++;
        }

        if(n_x < n_control){
          set_base_weights(p,target_year,dates,n_dates,input,weight);
          continue;
        }
        /* check number of unique x values */
        int n_unique = 1;
        for(int i=0;i<n_x;i++){
          int unique = 1;
          for(int j=0;j<i;j++){
              if(fabs(x[i]-x[j]) < 1e-6){
                  unique = 0;
                  break;
              }
          }
          if(unique) n_unique++;
      }

        if(n_unique < order){
            set_base_weights(p,target_year,dates,n_dates,input,weight);
            continue;
        }
        // first b-spline fit with fixed year's weights
        double chisq;
        gsl_vector_view xv = gsl_vector_view_array(x,n_x);
        gsl_vector_view yv = gsl_vector_view_array(y,n_x);
        gsl_vector_view wv = gsl_vector_view_array(w,n_x);
        
        // if any error happens, use the year weights
        gsl_set_error_handler_off();
        int status = gsl_bspline_wlssolve(&xv.vector,&yv.vector,&wv.vector,c,&chisq,work);

        if(status != GSL_SUCCESS){
          set_base_weights(p,target_year,dates,n_dates,input,weight);
          continue;
        }
        //gsl_bspline_wlssolve(&xv.vector,&yv.vector,&wv.vector,c,&chisq,work);

        // calculate residuals and adjust weights with penalty
        for(int i=0;i<n_x;i++){

            double est;
            gsl_bspline_calc(x[i],c,&est,work);
            double res = fabs(y[i]-est);

            if(fabs(est)<1e-6) est = 1e-6;
            double rel = res/fabs(est);
            double penalty;

            // bestimmung des panalties, z.B. linear zwischen 0.1 und 0.5, darüber 0.5, darunter 1.0
            double max_rel = 0.4;
            if(rel<=0.1) penalty = 1.0;
            else if(rel>=max_rel) penalty = 0.1; 
            //Wenn die abweichung größer als 40% ist, wird das Gewicht auf 10% reduziert
            else{
                double t = (rel-0.1)/(max_rel-0.1);
                penalty = 1.0 - 0.9*t;
            }

            w[i] *= penalty;

            int d = obs_index[i];
            if(d >= 0 && d < weight->nb){   
              weight->data[d][p] = w[i] * 10000;
            }

        }
    }

    gsl_vector_free(c);
    gsl_vector_free(tmpB);
    gsl_bspline_free(work);

    free(x);
    free(y);
    free(w);
    free(obs_index);              // NEW

  }
  return SUCCESS;
}
//####################################

int main( int argc, char *argv[] ){
time_t TIME;


time(&TIME);

  int n_pad = 3; // number of padded points before first and after last date for stabilization of spline at edges
  double w_pad = 0.5; // weight of padded points (between 0 and 1, smaller means less influence on fit)

  args_t args;
  parse_args(argc, argv, &args);

  GDALAllRegister();

  omp_set_num_threads(args.n_cpus);


  image_t mask;
  read_image(args.path_mask, NULL, &mask, args.pixel_size);

  table_t input_table = read_str_table(args.path_input, false, false);

  if (input_table.ncol != N_PROD){
    fprintf(stderr, "Error: number of columns in input table must be 3 (BOA, QAI, THT)\n");
    return FAILURE;
  }

  int n_dates = input_table.nrow;
  if (n_dates < args.n_control_points){
    fprintf(stderr, "Error: at least %d rows (dates) in input table are required.\n", args.n_control_points);
    return FAILURE;
  }

  int n_years = -1;
  date_t **dates = NULL;
  image_t **input = NULL;
  alloc_2D((void***)&input, n_dates, N_PROD, sizeof(image_t));
  alloc_2D((void***)&dates, n_dates, N_PROD, sizeof(date_t));

  for (int i=0; i<n_dates; i++){

    for (int j=0; j<N_PROD; j++){

      char basename[STRLEN];
      basename_with_ext(input_table.str_data[i][j], basename, STRLEN);
      date_from_string(&dates[i][j], basename);
  
      if (i > 0 && j == 0){
        if (dates[i][j].ce < dates[i-1][j].ce){
          fprintf(stderr, "Input images must be ordered by date (earliest to latest).\n");
          exit(FAILURE);
        }
        if (dates[i][j].year > args.target_year){
          fprintf(stderr, "Input images must not exceed target year.\n");
          exit(FAILURE);
        }
      }

      if (j == 0){
        int diff_years = args.target_year - dates[i][j].year + 1;
        if (diff_years > n_years) n_years = diff_years;
      }

      if (j > 0){
        if (dates[i][j].ce != dates[i][0].ce){
          fprintf(stderr, "All columns in input table must have the same date for each row.\n");
          exit(FAILURE);
        }
      }
      
    }
    
  }

  if (n_years < 1){
    fprintf(stderr, "Error: no valid years found in input data.\n");
    return FAILURE;
  }
  

  #pragma omp parallel shared(input_table, dates, n_dates, input, mask, args, stderr) default(none)
  {

    #pragma omp for
    for (int i=0; i<n_dates; i++){

      for (int j=0; j<N_PROD; j++){

        read_image(input_table.str_data[i][j], NULL, &input[i][j], args.pixel_size);
        compare_images(&mask, &input[i][j]);

        if (i > 0 && j == BOA){
          if (input[i][j].nb < input[i-1][j].nb){
            fprintf(stderr, "Input images must have the same number of bands.\n");
            fprintf(stderr, "Image %s has fewer bands than previous image %s.\n", input_table.str_data[i][j], input_table.str_data[i-1][j]);
            exit(FAILURE);
          }
        }
        
      }
      
    }
    
  }

  int n_bands = input[0][BOA].nb;
  
  if (args.band_nir >= n_bands){
    fprintf(stderr, "Error: NIR band is out of bounds for the input images: %d.\n", n_bands);
    return FAILURE;
  }

  if (args.band_red >= n_bands){
    fprintf(stderr, "Error: RED band is out of bounds for the input images: %d.\n", n_bands);
    return FAILURE;
  }

  proctime_print("Reading time", TIME);
  
  image_t weights;
  //copy_image(&mask, &weights, n_years, SHRT_MIN, "/home/ahsoka/frantz/temp/weights.tif");
  //###############################
  // "alte" geupdatete weights version
  //copy_image(&mask, &weights, n_dates, SHRT_MIN, "/home/ahsoka/frantz/temp/weights.tif");
  //###############################
  // nach git merge eingefügte masken version
  copy_image(&mask, &weights, n_years, SHRT_MIN, "NULL");
  

  //if (determine_annual_weights(4, 6, args.max_weight, args.target_year, 
  //      dates, n_dates, n_years, input, args.band_nir, args.band_red, &mask, &weights) != SUCCESS){
  //  fprintf(stderr, "Error determining annual weights.\n");
  //  return FAILURE;
  //}
  //###############################
  if (determine_point_weights(args.order,args.n_control_points,args.target_year,
        dates,n_dates,input,args.band_nir,args.band_red,&mask,&weights) != SUCCESS){
    fprintf(stderr,"Error determining weights\n");
    return FAILURE;
  }
  //###############################
  
  //write_image(&weights);

  proctime_print("Weighting time", TIME);

  image_t coefficients;
  copy_image(&mask, &coefficients, args.n_control_points * n_bands, SHRT_MIN, args.path_output);


  
  
  
  #pragma omp parallel shared(mask, weights, input, dates, n_dates, n_pad, w_pad, args, n_bands, coefficients) default(none)
  {

    double min_x = 0 - n_pad;
    double max_x = 365 + n_pad;

    gsl_vector *c = gsl_vector_calloc(args.n_control_points);

    gsl_bspline_workspace *work;
    work = gsl_bspline_alloc_ncontrol(args.order, args.n_control_points);
    gsl_bspline_init_uniform(min_x, max_x, work);
    
    gsl_matrix *AtWA = gsl_matrix_calloc(args.n_control_points, args.n_control_points); // normal matrix
    gsl_vector *AtWy = gsl_vector_calloc(args.n_control_points); // right-hand side
    gsl_vector *tmpB = gsl_vector_calloc(args.n_control_points); // temp for basis
    
    int n_x_max = n_dates + n_pad*2;

    double *x = NULL;
    double **y = NULL;
    double *w = NULL;

    alloc((void**)&x, n_x_max, sizeof(double));
    alloc_2D((void***)&y, n_x_max, n_bands, sizeof(double));
    alloc((void**)&w, n_x_max, sizeof(double));
    
    // Allocate B-spline basis matrix A (size n_dates x n_control)
    double **A = NULL;
    alloc_2D((void ***)&A, n_x_max, args.n_control_points, sizeof(double)); // B-spline basis matrix
    
    double **P = NULL;
    alloc_2D((void***)&P, args.n_control_points, args.n_control_points, sizeof(double)); // penalty matrix

    // Build penalty matrix P (integral of product of 2nd derivatives of basis functions)
    // Approximate by finite differences over a fine grid
    int n_grid = 365; // number of points in grid for numerical integration
    double dx = (max_x - min_x) / (n_grid - 1);
    for (int i=0; i<args.n_control_points; i++) {
      for (int j=0; j<args.n_control_points; j++) {
        double integral = 0.0;
        for (int g=0; g<n_grid; g++) {
          double xg = min_x + g * dx;
          double eps = 1e-5;
          // Clamp xg+eps and xg-eps to [min_x, max_x]
          double xg_peps = xg + eps;
          double xg_meps = xg - eps;
          if (xg_peps > max_x) xg_peps = max_x;
          if (xg_peps < min_x) xg_peps = min_x;
          if (xg_meps > max_x) xg_meps = max_x;
          if (xg_meps < min_x) xg_meps = min_x;
          gsl_bspline_eval(xg_peps, tmpB, work);
          double f1 = tmpB->data[i]; // gsl_vector_get(tmpB, i);
          double f2 = tmpB->data[j];
          gsl_bspline_eval(xg, tmpB, work);
          double f0i = tmpB->data[i];
          double f0j = tmpB->data[j];
          gsl_bspline_eval(xg_meps, tmpB, work);
          double fm1i = tmpB->data[i];
          double fm1j = tmpB->data[j];
          double d2Bi = (f1 - 2*f0i + fm1i) / (eps*eps);
          double d2Bj = (f2 - 2*f0j + fm1j) / (eps*eps);
          integral += d2Bi * d2Bj * dx;
        }
        P[i][j] = integral;
      }
    }

    #pragma omp for
    for (int p=0; p<mask.nc; p++){

      if (mask.data[0][p] == mask.nodata || mask.data[0][p] == 0) continue;

      int n_x = 0;
      int target_year_start = -1;
      int target_year_end = -1;

      for (int i=0; i<n_dates; i++){

        if (dates[i][BOA].year > args.target_year) continue;
        if (!use_this_pixel(input[i][QAI].data[0][p])) continue;

        x[n_x] = input[i][THT].data[0][p];
        for (int b=0; b<n_bands; b++){
          y[n_x][b] = input[i][BOA].data[b][p];
        }

        double denom = input[i][BOA].data[args.band_nir][p] + input[i][BOA].data[args.band_red][p];
        if (denom == 0) denom = 1;
        double ndvi = (input[i][BOA].data[args.band_nir][p] - input[i][BOA].data[args.band_red][p]) / denom;
        if (ndvi < 0) ndvi = 0.01;
        //######################
        //w[n_x] = sqrt(ndvi) * weights.data[args.target_year - dates[i][BOA].year][p]/10000.0;
        w[n_x] = sqrt(ndvi) * weights.data[i][p] / 10000.0;
        //######################
        if (dates[i][BOA].year == args.target_year){
          if (target_year_start < 0) target_year_start = n_x;
          target_year_end = n_x;
        }

        n_x++;

      }

      if (target_year_start < 0 || target_year_end < 0){
        for (int c=0; c<args.n_control_points*n_bands; c++) coefficients.data[c][p] = coefficients.nodata; // set coefficients to nodata for this pixel
        continue; // no data for target year
      } 

      // old Padding: just take the first and last point and add n_pad points before and after with same value and weight
      /*for (int i=0; i<n_pad; i++){

        // Pad before first date
        x[n_x] = -1 - i;
        for (int b=0; b<n_bands; b++){
          y[n_x][b] = y[target_year_start][b];
        }
        w[n_x] = w[target_year_start]*w_pad; // same weight as first point
        n_x++;

        // Pad after last date
        x[n_x] = 366 + i;
        for (int b=0; b<n_bands; b++){
          y[n_x][b] = y[target_year_end][b];
        }
        w[n_x] = w[target_year_end]*w_pad; // same weight as last point
        n_x++;

      }*/
      // new padding: take the median of the first/last 3 points for value and weight, to be more robust against outliers at the edges, and add n_pad points before and after with same value and weight
      for (int i=0; i<n_pad; i++){
        // --- Pad before first date: Median der ersten 3 Observations ---
        x[n_x] = -1 - i;
        for (int b=0; b<n_bands; b++){
            double vals[3];
            int count = 0;
            for (int j=0; j<3 && (target_year_start+j)<n_x; j++){
                vals[count++] = y[target_year_start + j][b];
            }
            // sortiere und nimm Median
            if(count==3){
                if(vals[0]>vals[1]) { double t=vals[0]; vals[0]=vals[1]; vals[1]=t; }
                if(vals[1]>vals[2]) { double t=vals[1]; vals[1]=vals[2]; vals[2]=t; }
                if(vals[0]>vals[1]) { double t=vals[0]; vals[0]=vals[1]; vals[1]=t; }
                y[n_x][b] = vals[1]; // Median
            } else if(count>0){
                y[n_x][b] = vals[count/2];
            } else {
                y[n_x][b] = y[target_year_start][b];
            }
        }
        w[n_x] = w[target_year_start]*w_pad;
        n_x++;

        // --- Pad after last date: Median der letzten 3 Observations ---
        x[n_x] = 366 + i;
        for (int b=0; b<n_bands; b++){
            double vals[3];
            int count = 0;
            for (int j=0; j<3 && (target_year_end-j)>=0; j++){
                vals[count++] = y[target_year_end - j][b];
            }
            // sortiere und nimm Median
            if(count==3){
                if(vals[0]>vals[1]) { double t=vals[0]; vals[0]=vals[1]; vals[1]=t; }
                if(vals[1]>vals[2]) { double t=vals[1]; vals[1]=vals[2]; vals[2]=t; }
                if(vals[0]>vals[1]) { double t=vals[0]; vals[0]=vals[1]; vals[1]=t; }
                y[n_x][b] = vals[1]; // Median
            } else if(count>0){
                y[n_x][b] = vals[count/2];
            } else {
                y[n_x][b] = y[target_year_end][b];
            }
        }
        w[n_x] = w[target_year_end]*w_pad;
        n_x++;
      }
      //new padding end


      if (n_x < args.n_control_points){
        for (int c=0; c<args.n_control_points*n_bands; c++) coefficients.data[c][p] = coefficients.nodata; // set coefficients to nodata for this pixel
        continue; // need at least n_control_points data points for fitting
      }

      // Build B-spline basis matrix A
      for (int i=0; i<n_x; i++) {
        gsl_bspline_eval(x[i], tmpB, work);
        for (int j=0; j<args.n_control_points; j++) {
          A[i][j] = tmpB->data[j];
        }
      }

      // Compute AtWA
      for (int i=0; i<args.n_control_points; i++) {
        for (int j=0; j<args.n_control_points; j++) {
          double sum = 0.0;
          for (int k=0; k<n_x; k++) {
            double wij = w[k] * A[k][i] * A[k][j];
            sum += wij;
          }
          AtWA->data[i*AtWA->tda+j] = sum;
        }
      }


      // Add lambda*P to AtWA
      for (int i=0; i<args.n_control_points; i++) {
        for (int j=0; j<args.n_control_points; j++) {
          AtWA->data[i*AtWA->tda+j] += args.lambda * P[i][j];
        }
      }

      gsl_linalg_cholesky_decomp(AtWA);

      
      for (int b=0; b<n_bands; b++){
        
        // Compute AtWy
        for (int i=0; i<args.n_control_points; i++) {
          double sum2 = 0.0;
          for (int k=0; k<n_x; k++) {
            sum2 += w[k] * A[k][i] * y[k][b];
          }
          AtWy->data[i] = sum2;
        }

        // Solve (AtWA + lambda*P) c = AtWy
        //gsl_linalg_LU_solve(AtWA, perm, AtWy, c);
        gsl_linalg_cholesky_solve(AtWA, AtWy, c);

        for (int i=0; i<args.n_control_points; i++) coefficients.data[i*n_bands + b][p] = (short)c->data[i];

      } // end band-loop
      
    } // end pixel-loop
    //##########################
    //free_2D((void**)A, n_dates);
    free_2D((void**)A, n_x_max);
    //#########################
    gsl_vector_free(c);
    gsl_bspline_free(work);
    gsl_matrix_free(AtWA);
    free_2D((void**)P, args.n_control_points);
    gsl_vector_free(AtWy);
    gsl_vector_free(tmpB);

  } // end parallel region


  proctime_print("Splines calculated", TIME);

  write_image(&coefficients);

  
  free_table(&input_table);
  for (int i=0; i<n_dates; i++){
    for (int j=0; j<N_PROD; j++){
       free_image(&input[i][j]);
    }
  }
  free_2D((void**)input, n_dates);
  free_image(&mask);
  free_image(&coefficients);
  free_image(&weights);
  free((void*)dates);

  proctime_print("Total time", TIME);

  return SUCCESS;
}

