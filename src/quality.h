/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Higher Level quality screening
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#ifndef QUALITY_H
#define QUALITY_H

#include <stdio.h>   // core input and output functions
#include <stdlib.h>  // standard general utilities library
#include <stdbool.h>  // boolean data type


#ifdef __cplusplus
extern "C" {
#endif

bool use_this_pixel(short qai);

#ifdef __cplusplus
}
#endif

#endif

