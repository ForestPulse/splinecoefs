#ifndef USAGE_H
#define USAGE_H

#include <stdio.h>   // core input and output functions
#include <stdlib.h>  // standard general utilities library
#include <stdbool.h> // boolean data type

#include <omp.h>

#include "dtype.h"
#include "string.h"

#ifdef __cplusplus
extern "C" {
#endif

void parse_args(int argc, char *argv[], args_t *args);

#ifdef __cplusplus
}
#endif

#endif
