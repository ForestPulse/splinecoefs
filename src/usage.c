#include "usage.h"


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

