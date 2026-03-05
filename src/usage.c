#include "usage.h"


void usage(char *exe, int exit_code){


  printf("Usage: %s [-h] [-o] [-p] [-f] [-v] [-s] [-r] [-n] [-d] [-j] input-image input-bands\n", exe);
  printf("\n");
  printf("  -h  = show this help\n");
  printf("\n");
  printf("  -o output-file  = output file path with extension,\n");
  printf("     defaults to 'sharpened.tif'\n");
  printf("  -p pca-file = output file path of PCA transformation,\n");
  printf("     when not given, file is not written\n");
  printf("  -f format  = output format (GDAL vector driver short name)\n");
  printf("     defaults to GTiff\n");
  printf("  -v variance = how much percent of the PCA-variance should be retained for the target bands?\n");
  printf("     defaults to 99\n");
  printf("  -s sampling = sampling factor to speed up computation of PCA\n");
  printf("     defaults to 10\n");
  printf("  -r radius = how many neighboring cells to use for sharpening?\n");
  printf("     defaults to 2\n");
  printf("  -n nbreaks = number of breaks for B-Spline\n");
  printf("     defaults to 10\n");
  printf("  -d order = order of the B-Spline\n");
  printf("     defaults to 4\n");
  printf("  -j ncpu = How many CPUs to use?\n");
  printf("     defaults to all\n");
  
  printf("\n");
  printf("  Positional arguments:\n");
  printf("  - input-image: well, the input image...\n");
  printf("  - input-bands: band definition\n");
  printf("     csv table [en], three (or more) named columns\n");
  printf("       band: band number\n");
  printf("       wavelength: band's central wavelength\n");
  printf("       use:  usage code\n");
  printf("         1: target band (highres)\n");
  printf("         2: spatial prediction band (lowres)\n");
  printf("         0: spectral prediction band (anyres)\n");
  printf("        -1: ignore, bad band\n");
  printf("\n");

  exit(exit_code);
  return;
}


void parse_args(int argc, char *argv[], args_t *args){
int opt;
bool o = false, f = false, p = false;


  opterr = 0;

  // default parameters
  args->ncpu = omp_get_max_threads();
  args->radius = 2;
  args->minvar = 0.99;
  args->sample = 10;
  args->nbreak = 10;
  args->order  = 4;
  copy_string(args->f_output, STRLEN, "sharpened.tif");
  copy_string(args->f_pca, STRLEN, "NULL");
  copy_string(args->format, STRLEN, "GTiff");

  // optional parameters
  while ((opt = getopt(argc, argv, "ho:f:j:r:v:p:s:n:d:")) != -1){
    switch(opt){
      case 'h':
        usage(argv[0], SUCCESS);
      case 'o':
        copy_string(args->f_output, STRLEN, optarg);
        o = true;
        break;
      case 'f':
        copy_string(args->format, STRLEN, optarg);
        f = true;
        break;
      case 'j':
        args->ncpu = atoi(optarg);
        break;
      case 'r':
        args->radius = atoi(optarg);
        break;
      case 'v':
        args->minvar = atof(optarg)/100.0;
        break;
      case 'p':
        copy_string(args->f_pca, STRLEN, optarg);
        p = true;
        break;
      case 's':
        args->sample = atoi(optarg);
        break;
      case 'n':
        args->nbreak = atoi(optarg);
        break;
      case 'd':
        args->order = atoi(optarg);
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

  // non-optional parameters
  args->n = 2;

  if (optind < argc){
    if (argc-optind == args->n){
      copy_string(args->f_input, STRLEN, argv[optind++]);
      copy_string(args->f_bands, STRLEN, argv[optind++]);
    } else if (argc-optind < args->n){
      fprintf(stderr, "some non-optional arguments are missing.\n");
      usage(argv[0], FAILURE);
    } else if (argc-optind > args->n){
      fprintf(stderr, "too many non-optional arguments.\n");
      usage(argv[0], FAILURE);
    }
  } else {
    fprintf(stderr, "non-optional arguments are missing.\n");
    usage(argv[0], FAILURE);
  }

  if ((!o && f) || (!f && o)){
    fprintf(stderr, "If -f is given, -o needs to be given, too.\n"); 
    usage(argv[0], FAILURE);
  }

  if (p && !f){
    fprintf(stderr, "If -p is given, -f needs to be given, too. Suggestion: -f GTiff\n"); 
    usage(argv[0], FAILURE);
  }

  return;
}

