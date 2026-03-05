#ifndef USAGE_H
#define USAGE_H

#include <stdio.h>   // core input and output functions
#include <stdlib.h>  // standard general utilities library
#include <unistd.h>
#include <ctype.h>

#include "dtype.h"
#include "dir.h"
#include "string.h"

#ifdef __cplusplus
extern "C" {
#endif

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

void usage(char *exe, int exit_code);
void parse_args(int argc, char *argv[], args_t *args);

#ifdef __cplusplus
}
#endif

#endif
