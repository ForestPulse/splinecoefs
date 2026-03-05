#include "read.h"


int read_dataset(img_t *images, table_t *bandlist, args_t *args){
GDALDatasetH dataset = NULL;
GDALRasterBandH band = NULL;
int b, b_highres, b_lowres;
int has_nodata, n_band;
int col_use, col_band, col_wavelength;
time_t TIME;

  
  time(&TIME);

  printf("Starting Image Read\n")  ;


  if ((dataset = GDALOpen(args->f_input, GA_ReadOnly)) == NULL){
    printf("unable to open %s\n", args->f_input);
    exit(FAILURE);
  }

  if ((n_band = GDALGetRasterCount(dataset)) != bandlist->nrow){
    printf("number of bands in bandlist (%d) and input image (%d) do not match\n",
      bandlist->nrow, n_band); 
    exit(FAILURE);
  }

  images[HIGHRES].meta.dim.col = GDALGetRasterXSize(dataset);
  images[HIGHRES].meta.dim.row = GDALGetRasterYSize(dataset);
  images[HIGHRES].meta.dim.cell = images[HIGHRES].meta.dim.col * images[HIGHRES].meta.dim.row;
  images[HIGHRES].meta.dim.band = 0;

  GDALGetGeoTransform(dataset, images[HIGHRES].meta.transformation);
  copy_string(images[HIGHRES].meta.projection, STRLEN, GDALGetProjectionRef(dataset));
  images[HIGHRES].meta.datatype = GDALGetDataTypeByName(dataset);
  
  memcpy(&images[LOWRES].meta, &images[HIGHRES].meta, sizeof(meta_t));

  //print_table(bandlist, false, false);

  if ((col_use  = find_table_col(bandlist, "use")) < 0){
    printf("there is no column 'use' in csv-file\n");
    exit(FAILURE);
  }
  if ((col_band = find_table_col(bandlist, "band")) < 0){
    printf("there is no column 'use' in csv-file\n");
    exit(FAILURE);
  }
  if ((col_wavelength = find_table_col(bandlist, "wavelength")) < 0){
    printf("there is no column 'wavelength' in csv-file\n");
    exit(FAILURE);
  }

  for (b=0; b<bandlist->nrow; b++){
    if ((int)bandlist->data[b][col_use] == 1) images[HIGHRES].meta.dim.band++;
    if ((int)bandlist->data[b][col_use] == 2) images[LOWRES].meta.dim.band++;
  }

  alloc_2D((void***)&images[HIGHRES].data, images[HIGHRES].meta.dim.band, images[HIGHRES].meta.dim.cell, sizeof(float));
  alloc_2D((void***)&images[LOWRES].data, images[LOWRES].meta.dim.band, images[LOWRES].meta.dim.cell, sizeof(float));

  for (b=0, b_highres=0, b_lowres=0; b<bandlist->nrow; b++){

    if ((int)bandlist->data[b][col_band] > n_band){
      printf("band %d in bandlist is higher than bands (%d) in dataset\n", (int)bandlist->data[b][col_band], n_band);
      exit(FAILURE);
    }

    if ((int)bandlist->data[b][col_band] < 1){
      printf("band %d in bandlist is smaller than 1\n", (int)bandlist->data[b][col_band]);
      exit(FAILURE);
    }

    //printf("read band # %d\n", (int)bandlist->data[b][col_band]);

    band = GDALGetRasterBand(dataset, (int)bandlist->data[b][col_band]);

    if (b == 0){
      images[HIGHRES].meta.nodata = (float) GDALGetRasterNoDataValue(band, &has_nodata);
      images[LOWRES].meta.nodata = images[HIGHRES].meta.nodata;
      if (!has_nodata){
        printf("input image has no nodata value in band %d.\n", (int)bandlist->data[b][col_band]); 
        exit(FAILURE);
      }
    }

    if ((int)bandlist->data[b][col_use] == 1){
      if (GDALRasterIO(band, GF_Read, 0, 0, images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, images[HIGHRES].data[b_highres++], 
        images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("could not read band #%d from %s.\n", b+1, args->f_input); 
        exit(FAILURE);
      }
    } else if ((int)bandlist->data[b][col_use] == 2){
      if (GDALRasterIO(band, GF_Read, 0, 0, images[LOWRES].meta.dim.col, images[LOWRES].meta.dim.row, images[LOWRES].data[b_lowres++], 
        images[LOWRES].meta.dim.col, images[LOWRES].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("could not read band #%d from %s.\n", b+1, args->f_input); 
        exit(FAILURE);
      }
    }

  }

  GDALClose(dataset);


  proctime_print("Reading", TIME);

	return SUCCESS;
}
