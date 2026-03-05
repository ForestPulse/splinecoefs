#ifndef READ_H
#define READ_H

/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h" // public (C callable) GDAL entry points

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "dtype.h"
#include "alloc.h"
#include "string.h"
#include "table.h"

#ifdef __cplusplus
extern "C" {
#endif

int read_dataset(img_t *images, table_t *bandlist, args_t *args);

#ifdef __cplusplus
}
#endif

#endif
