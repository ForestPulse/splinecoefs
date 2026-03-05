#ifndef WRITE_H
#define WRITE_H

/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h"  // public (C callable) GDAL entry points
#include "cpl_string.h"  // Various convenience functions for working with strings and string lists


#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include "dtype.h"
#include "table.h"

int write_pca(img_t *images, args_t *args);
int write_output(img_t *images, table_t *bandlist, args_t *args);

#ifdef __cplusplus
}
#endif

#endif
