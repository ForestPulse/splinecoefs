/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Image I/O header
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <stdio.h>    // core input and output functions
#include <stdlib.h>   // standard general utilities library
#include <string.h>   // string handling functions
#include <stdbool.h>  // boolean data type

#include "alloc.h"
#include "const.h"
#include "string.h"

/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h"       // public (C callable) GDAL entry points
#include "cpl_conv.h"   // various convenience functions for CPL
#include "cpl_string.h" // various convenience functions for strings


#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  int *number;    // array of band numbers
  float *wavelengths; // array of band wavelengths
  int n;   // number of bands
} bandlist_t;

typedef struct {
  char path[STRLEN];    // file path
  char proj[STRLEN];    // directory name
  double geotran[6];    // geotransform
  int nx, ny, nc, nb;   // dimensions
  short **data;
  short nodata;
} image_t;

void read_image(char *path, bandlist_t *bands, image_t *image);
void copy_image(image_t *from, image_t *to, int nbands, short nodata, char *path);
void write_image(image_t *image);
void free_image(image_t *image);
void compare_images(image_t *image_1, image_t *image_2);

#ifdef __cplusplus
}
#endif

#endif

