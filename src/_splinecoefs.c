#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>

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



typedef struct {
  int n_cpus;
  char path_input[STRLEN];
  char path_mask[STRLEN];
  char path_output[STRLEN];
  int order;
  int n_control_points;
  double lambda;
  int target_year;
  double max_weight;
  int band_nir;
  int band_red;
} args_t;


// -i input-table = csv file with input images
// 3 columns: BOA, QAI, THT
// each row corresponds to one date, and the 3 columns must have the same date
// dates should be ordered from earliest to latest
// file basenames must start with date in format YYYYMMDD, e.g., 20200731_LEVEL2_SEN2A_BOA.tif


void usage(char *exe, int exit_code){
  printf("Usage: %s -j cpus -i input-table -x mask-image -o output-image\n", exe);
  printf("          -k order -c control-points -l lambda -y target-year\n");
  printf("          -r = RED-band -n = NIR-band\n");
  printf("\n");
  printf("  -j = number of CPUs to use\n");
  printf("\n");
  printf("  -i = input-table = csv file with input images\n");
  printf("  -x = mask image\n");
  printf("  -o = output file (.tif)\n");
  printf("\n");  
  printf("  -k = order of the spline\n");
  printf("  -c = number of control points (equiv. to number of output coefficients)\n");
  printf("  -l = lambda (smoothing parameter)\n");
  printf("  -y = target year\n");
  printf("  -w = maximum weight for past data points\n");
  printf("  -r = RED band number\n");
  printf("  -n = NIR band number\n");
  printf("\n");
  exit(exit_code);
  return;
}

void parse_args(int argc, char *argv[], args_t *args){
  int opt, received_n = 0, expected_n = 11;
  opterr = 0;

  while ((opt = getopt(argc, argv, "j:i:x:o:k:c:l:y:w:r:n:")) != -1){
    switch(opt){
      case 'j':
        args->n_cpus = atoi(optarg);
        received_n++;
        break;
      case 'i':
        copy_string(args->path_input, STRLEN, optarg);
        received_n++;
        break;
      case 'x':
        copy_string(args->path_mask, STRLEN, optarg);
        received_n++;
        break;
      case 'o':
        copy_string(args->path_output, STRLEN, optarg);
        received_n++;
        break;
      case 'k':
        args->order = atoi(optarg);
        received_n++;
        break;
      case 'c':
        args->n_control_points = atoi(optarg);
        received_n++;
        break;
      case 'l':
        args->lambda = atof(optarg);
        received_n++;
        break;
      case 'y':
        args->target_year = atoi(optarg);
        received_n++;
        break;
      case 'w':
        args->max_weight = atof(optarg);
        received_n++;
        break;
      case 'r':
        args->band_red = atoi(optarg) - 1;
        received_n++;
        break;
      case 'n':
        args->band_nir = atoi(optarg) - 1;
        received_n++;
        break;
      case '?':
        if (isprint(optopt)){
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        } else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        usage(argv[0], FAILURE);
      default:
        fprintf(stderr, "Error parsing arguments.\n");
        usage(argv[0], FAILURE);
    }
  }

  if (received_n != expected_n){
    fprintf(stderr, "Not all arguments received.\n");
    usage(argv[0], FAILURE);
  }

  if (!fileexist(args->path_input)){
    fprintf(stderr, "Input file %s does not exist.\n", args->path_input);
    usage(argv[0], FAILURE);
  }
  
  if (!fileexist(args->path_mask)){
    fprintf(stderr, "Mask file %s does not exist.\n", args->path_mask);
    usage(argv[0], FAILURE);
  }

  if (fileexist(args->path_output)){
    fprintf(stderr, "Output file %s already exists.\n", args->path_output);
    usage(argv[0], FAILURE);
  }

  if (args->n_cpus < 1){
    fprintf(stderr, "Number of CPUs must be at least 1.\n");
    usage(argv[0], FAILURE);
  }
  
  if (args->order < 1){
    fprintf(stderr, "order must be at least 1.\n");
    usage(argv[0], FAILURE);
  }

  if (args->n_control_points < 1){
    fprintf(stderr, "number of control points must be at least 1.\n");
    usage(argv[0], FAILURE);
  }

  if (args->lambda < 0){
    fprintf(stderr, "lambda must be non-negative.\n");
    usage(argv[0], FAILURE);
  }

  if (args->target_year < 1960){
    fprintf(stderr, "target year must be >= 1960.\n");
    usage(argv[0], FAILURE);
  }

  if (args->target_year > 2099){
    fprintf(stderr, "target year must be <= 2099.\n");
    usage(argv[0], FAILURE);
  }

  if (args->max_weight < 0){
    fprintf(stderr, "maximum weight must be non-negative.\n");
    usage(argv[0], FAILURE);
  }

  if (args->max_weight > 1){
    fprintf(stderr, "maximum weight must be <= 1.\n");
    usage(argv[0], FAILURE);
  }

  if (args->band_red < 0){
    fprintf(stderr, "RED band number must be at least 1.\n");
    usage(argv[0], FAILURE);
  }

  if (args->band_nir < 0){
    fprintf(stderr, "NIR band number must be at least 1.\n");
    usage(argv[0], FAILURE);
  }

  return;
}



int determine_annual_weights(int order, int n_control, double max_weight, int target_year, date_t **dates, int n_dates, int n_years, image_t **input, int band_nir, int band_red, image_t *mask, image_t *weight){

  //image_t predicted_img[3];
  //copy_image(mask, &predicted_img[0], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/predicted_0.tif");
  //copy_image(mask, &predicted_img[1], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/predicted_1.tif");
  //copy_image(mask, &predicted_img[2], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/predicted_2.tif");

  //image_t original_img[3];
  //copy_image(mask, &original_img[0], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/original_0.tif");
  //copy_image(mask, &original_img[1], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/original_1.tif");
  //copy_image(mask, &original_img[2], 52, SHRT_MIN, "/home/ahsoka/frantz/temp/original_2.tif");

  #pragma omp parallel shared(order, n_control, max_weight, target_year, dates, n_dates, n_years, input, band_nir, band_red, mask, weight, stderr) default(none)
  {

    double min_x = 0.0;
    double max_x = 365.0;
    
    enum { START, END, RANGE };

    int n_doys = 365;

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

    //int *valid_years = NULL;
    //alloc((void**)&valid_years, n_years, sizeof(int));

    double **predicted = NULL;
    alloc_2D((void***)&predicted, n_years, n_doys, sizeof(double));

    double *mean = NULL;
    alloc((void**)&mean, n_years, sizeof(double));

    #pragma omp for
    for (int p=0; p<mask->nc; p++){

      if (mask->data[0][p] == mask->nodata || mask->data[0][p] == 0) continue;

      // reset variables and counter for each year
      for (int year=0; year<n_years; year++){
        for (int d=0; d<vector_size; d++){
          x[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
          y[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
          w[year]->data[d] = 0; // reset size of vector to 0 for each year, will be incremented as data points are added
        }
        n[year] = 0;
      }

      for (int d=0; d<n_dates; d++){

        if (!use_this_pixel(input[d][QAI].data[0][p])) continue;
/**
        if (d > 0 && dates[d][BOA].year == dates[d-1][BOA].year &&
          (dates[d][BOA].doy - dates[d-1][BOA].doy) > (365/n_control)){
          x->data[k] = x->data[k-1] + 365 / n_control / 2; // add a point at regular interval for gap
          y->data[k] = y->data[k-1]; // repeat last value for gap
          w->data[k] = w->data[k-1] * 0.1; // reduce weight for gap point
          k++;
        }

        if (d > 0 && dates[d][BOA].year != dates[d-1][BOA].year &&
          dates[d][BOA].doy > (365/n_control)){
          x->data[k] = 365 / n_control / 2; // add a point at regular interval for gap
          y->data[k] = 0.5; // repeat last value for gap
          w->data[k] = w->data[k-1] * 0.01; // reduce weight for gap point
          k++;
        }
**/

        int year = target_year - dates[d][BOA].year;

        x[year]->data[n[year]] = input[d][THT].data[0][p];
        double denom = input[d][BOA].data[band_nir][p] + input[d][BOA].data[band_red][p];
        if (denom == 0) denom = 1;
        double ndvi = (input[d][BOA].data[band_nir][p] - input[d][BOA].data[band_red][p]) / denom;
        y[year]->data[n[year]] = ndvi;
        w[year]->data[n[year]] = sqrt(ndvi);

        n[year]++;
/**
        if (d != (n_dates-1) && dates[d][BOA].year != dates[d+1][BOA].year &&
          (365 - dates[d][BOA].doy) > (365/n_control)){
          x->data[k] = dates[d][BOA].doy + 365 / n_control / 2; // add a point at regular interval for gap
          y->data[k] = y->data[k-1]; // repeat last value for gap
          w->data[k] = w->data[k-1] * 0.1; // reduce weight for gap point
          k++;
        }

        if (d == (n_dates-1) && 
          (365 - dates[d][BOA].doy) > (365/n_control)){
          x->data[k] = dates[d][BOA].doy + 365 / n_control / 2; // add a point at regular interval for gap
          y->data[k] = y->data[k-1]; // repeat last value for gap
          w->data[k] = w->data[k-1] * 0.1; // reduce weight for gap point
          k++;
        }
**/
      }

if (p == 3600000) printf("Pixel %d: year 0 has %d data points, year 1 has %d data points, year 2 has %d data points.\n", p, n[0], n[1], n[2]);
            
      for (int year=0; year<n_years; year++){
      
        if (n[year] < n_control){
          fprintf(stderr, "Not enough data points for pixel %din year %d: %d points, but %d control points required.\n", p, year, n[year], n_control);
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

        if (p == 4797133){
          printf("Pixel %d, year %d: added %d points for gaps, total now %d.\n", p, year, added_points, n[year]);
          for (int d=0; d<n[year]; d++){
            printf("  Point %d: x = %f, y = %f, w = %f\n", d, x[year]->data[d], y[year]->data[d], w[year]->data[d]);
          }
        }

        //gsl_vector_view x_sub;
        //gsl_vector_view y_sub;
        //gsl_vector_view w_sub;
        
        //printf("Year %d for pixel %d: using data points from index %d to %d.\n", target_year - year, p, year_indices[year][START], year_indices[year][END]);

        //x_sub = gsl_vector_subvector(x, year_indices[year][START], year_indices[year][END] - year_indices[year][START] + 1);
        //y_sub = gsl_vector_subvector(y, year_indices[year][START], year_indices[year][END] - year_indices[year][START] + 1);
        //w_sub = gsl_vector_subvector(w, year_indices[year][START], year_indices[year][END] - year_indices[year][START] + 1);
/**
        int n[52] = {0};
        for (int i=0; i<x_sub.vector.size; i++){
          int week = (int)(x_sub.vector.data[i] / 7);
          if (week >= 52) week = 51; // cap at week 51 to avoid out-of-bounds

          if (p == 556) printf("Adding value %f to week %d for year %d, pixel %d.\n", y_sub.vector.data[i], week, target_year - year, p);
          original_img[year].data[week][p] += y_sub.vector.data[i]*10000;
          n[week]++;
        }
        for (int w=0; w<52; w++){
          if (n[w] > 0){
            original_img[year].data[w][p] /= n[w];
          }
        }
**/

        double chisq, est;
        gsl_bspline_wlssolve(x[year], y[year], w[year], c, &chisq, work);
        //gsl_bspline_wlssolve(&x_sub.vector, &y_sub.vector, &w_sub.vector, c, &chisq, work);

        mean[year] = 0.0;
        for (int doy=0; doy<n_doys; doy++){
          gsl_bspline_calc(doy, c, &est, work);
          predicted[year][doy] = est;
          //if (doy % 7 == 3) { // only store every 30th day for debugging
          //  predicted_img[year].data[(int)(doy/7)][p] = (short)(est * 10000); // scale to 0-10000 for storage as short
          //}
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

        double w = max_weight * (1.0 - dist / dist_max);

        if (w < 0){
          weight->data[year][p] = 0;
          //fprintf(stderr, "Negative weight for year %d for pixel %d, setting to 0, d_max = %f, dist = %f.\n", target_year - year, p, dist_max, dist);
        } else {
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


  //write_image(&predicted_img[0]);
  //write_image(&predicted_img[1]);
  //write_image(&predicted_img[2]);

  //write_image(&original_img[0]);
  //write_image(&original_img[1]);
  //write_image(&original_img[2]);

  return SUCCESS;
}




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
  read_image(args.path_mask, NULL, &mask);

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
  

  #pragma omp parallel shared(input_table, dates, n_dates, input, mask, stderr) default(none)
  {

    #pragma omp for
    for (int i=0; i<n_dates; i++){

      for (int j=0; j<N_PROD; j++){

        read_image(input_table.str_data[i][j], NULL, &input[i][j]);
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
  copy_image(&mask, &weights, n_years, SHRT_MIN, "/home/ahsoka/frantz/temp/weights.tif");
  
  if (determine_annual_weights(4, 6, args.max_weight, args.target_year, 
        dates, n_dates, n_years, input, args.band_nir, args.band_red, &mask, &weights) != SUCCESS){
    fprintf(stderr, "Error determining annual weights.\n");
    return FAILURE;
  }
  
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

    //printf("Computing penalty matrix P...\n");
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

        w[n_x] = sqrt(ndvi) * weights.data[args.target_year - dates[i][BOA].year][p]/10000.0;
        n_x++;

        if (dates[i][BOA].year == args.target_year){
           if (target_year_start < 0) target_year_start = i;
           target_year_end = i;
        }

      }

      if (target_year_start < 0 || target_year_end < 0) continue; // no data for target year

      for (int i=0; i<n_pad; i++){

        // Pad before first date
        x[n_x] = -1 - i;
        for (int b=0; b<n_bands; b++){
          y[n_x][b] = input[target_year_start][BOA].data[b][p];
        }
        w[n_x] = w_pad; // same weight as first point
        n_x++;

        // Pad after last date
        x[n_x] = 366 + i;
        for (int b=0; b<n_bands; b++){
          y[n_x][b] = input[target_year_end][BOA].data[b][p];
        }
        w[n_x] = w_pad; // same weight as last point
        n_x++;

      }

      if (n_x < args.n_control_points) continue; // need at least 3 data points for fitting


      //printf("B-spline basis matrix A:\n");
      // Build B-spline basis matrix A
      for (int i=0; i<n_x; i++) {
        gsl_bspline_eval(x[i], tmpB, work);
        for (int j=0; j<args.n_control_points; j++) {
          A[i][j] = tmpB->data[j];
        }
      }

      //printf("Computing AtWA...\n");
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


      //printf("Solving linear system...\n");
      // Add lambda*P to AtWA
      for (int i=0; i<args.n_control_points; i++) {
        for (int j=0; j<args.n_control_points; j++) {
          AtWA->data[i*AtWA->tda+j] += args.lambda * P[i][j];
        }
      }

      //int signum;
      //gsl_linalg_LU_decomp(AtWA, perm, &signum);
      gsl_linalg_cholesky_decomp(AtWA);

      
      for (int b=0; b<n_bands; b++){
        
        //printf("Computing AtWy...\n");
        // Compute AtWy
        for (int i=0; i<args.n_control_points; i++) {
          double sum2 = 0.0;
          for (int k=0; k<n_x; k++) {
            sum2 += w[k] * A[k][i] * y[k][b];
          }
          AtWy->data[i] = sum2;
        }



        //printf("AtWA + lambda*P:\n");
        // Solve (AtWA + lambda*P) c = AtWy
        //gsl_linalg_LU_solve(AtWA, perm, AtWy, c);
        gsl_linalg_cholesky_solve(AtWA, AtWy, c);

        for (int i=0; i<args.n_control_points; i++) coefficients.data[i*n_bands + b][p] = (short)c->data[i];

        if (p==80564 && b == 8){
          // Print coefficients and fitted values
          printf("\nSmoothing spline coefficients (lambda=%g):\n", args.lambda);
          for (int i=0; i<args.n_control_points; i++) {
            printf("c [%02d] = %g\n", i, c->data[i]);
          }
          printf("\nFitted values:\n");
          for (int i=0; i<n_x; i++) {
            gsl_bspline_eval(x[i], tmpB, work);
            double est = 0.0;
            for (int j=0; j<args.n_control_points; j++) {
              est += c->data[j] * tmpB->data[j];
            }
            printf("%g %g %g %g\n", x[i], w[i], y[i][b], est);
          }
        }


      } // end band-loop
      
    } // end pixel-loop
    
    free_2D((void**)A, n_dates);
    gsl_vector_free(c);
    gsl_bspline_free(work);
    gsl_matrix_free(AtWA);
    free_2D((void**)P, args.n_control_points);
    gsl_vector_free(AtWy);
    gsl_vector_free(tmpB);

  } // end parallel region




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

