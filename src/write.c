#include "write.h"



int write_pca(img_t *images, args_t *args){
GDALDatasetH file = NULL;
GDALRasterBandH band = NULL;
GDALDriverH driver = NULL;
char **options = NULL;
int b;
time_t TIME;

  
  time(&TIME);

  printf("Starting PCA Write\n")  ;

  if (strcmp(args->f_pca, "NULL") == 0) return(SUCCESS);

  if ((driver = GDALGetDriverByName(args->format)) == NULL){
    printf("%s driver not found\n", args->format);
    exit(FAILURE);
  }

  if (strcmp(args->format, "GTiff") == 0){
    options = CSLSetNameValue(options, "TILED", "YES");
    options = CSLSetNameValue(options, "COMPRESS", "LZW");
    options = CSLSetNameValue(options, "PREDICTOR", "2");
    options = CSLSetNameValue(options, "INTERLEAVE", "BAND");
    options = CSLSetNameValue(options, "BIGTIFF", "YES");
  }

  if ((file = GDALCreate(driver, args->f_pca, images[PCA].meta.dim.col, images[PCA].meta.dim.row, images[PCA].meta.dim.band, GDT_Float32, options)) == NULL){
    printf("Error creating file %s. ", args->f_pca);
    exit(FAILURE);
  }

  for (b=0; b<images[PCA].meta.dim.band; b++){

    band = GDALGetRasterBand(file, b+1);

    if (GDALRasterIO(band, GF_Write, 0, 0, 
          images[PCA].meta.dim.col, images[PCA].meta.dim.row, images[PCA].data[b], 
          images[PCA].meta.dim.col, images[PCA].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
      printf("Unable to write a band into %s. ", args->f_pca);
      exit(FAILURE);
    }

    GDALSetDescription(band, "band name here");
    GDALSetRasterNoDataValue(band, images[PCA].meta.nodata);

  }


  #pragma omp critical
  {
    GDALSetGeoTransform(file, images[PCA].meta.transformation);
    GDALSetProjection(file, images[PCA].meta.projection);
  }

  GDALClose(file);

  CSLDestroy(options);

  proctime_print("writing PCA", TIME);

  return SUCCESS;
}


int write_output(img_t *images, table_t *bandlist, args_t *args){
GDALDatasetH file = NULL;
GDALRasterBandH band = NULL;
GDALDriverH driver = NULL;
char **options = NULL;
int b_list, b_highres, b_sharpened, b_spectralfit, b_output;
int col_use;
time_t TIME;

  
  time(&TIME);


  printf("Starting Image Write\n")  ;

  col_use  = find_table_col(bandlist, "use");


  if ((driver = GDALGetDriverByName(args->format)) == NULL){
    printf("%s driver not found\n", args->format);
    exit(FAILURE);
  }

  if (strcmp(args->format, "GTiff") == 0){
    options = CSLSetNameValue(options, "TILED", "YES");
    options = CSLSetNameValue(options, "COMPRESS", "LZW");
    options = CSLSetNameValue(options, "PREDICTOR", "2");
    options = CSLSetNameValue(options, "INTERLEAVE", "BAND");
    options = CSLSetNameValue(options, "BIGTIFF", "YES");
  }

  if ((file = GDALCreate(driver, args->f_output, images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, images[HIGHRES].meta.dim.band+images[SHARPENED].meta.dim.band+images[SPECTRALFIT].meta.dim.band, GDT_Int16, options)) == NULL){
    printf("Error creating file %s. ", args->f_output);
    exit(FAILURE);
  }

  for (b_list=0, b_highres=0, b_sharpened=0, b_spectralfit=0, b_output=1; b_list<bandlist->nrow; b_list++){

    if ((int)bandlist->data[b_list][col_use] == 1){

      band = GDALGetRasterBand(file, b_output++);

      if (GDALRasterIO(band, GF_Write, 0, 0, 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, images[HIGHRES].data[b_highres++], 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("Unable to write a band into %s. ", args->f_output);
        exit(FAILURE);
      }

      GDALSetDescription(band, "band name here");
      GDALSetRasterNoDataValue(band, images[HIGHRES].meta.nodata);

    } else if ((int)bandlist->data[b_list][col_use] == 2){

      band = GDALGetRasterBand(file, b_output++);

      if (GDALRasterIO(band, GF_Write, 0, 0, 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, images[SHARPENED].data[b_sharpened++], 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("Unable to write a band into %s. ", args->f_output);
        exit(FAILURE);
      }

      GDALSetDescription(band, "band name here");
      GDALSetRasterNoDataValue(band, images[HIGHRES].meta.nodata);

    } else if  ((int)bandlist->data[b_list][col_use] == 0){

      band = GDALGetRasterBand(file, b_output++);

      if (GDALRasterIO(band, GF_Write, 0, 0, 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, images[SPECTRALFIT].data[b_spectralfit++], 
            images[HIGHRES].meta.dim.col, images[HIGHRES].meta.dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("Unable to write a band into %s. ", args->f_output);
        exit(FAILURE);
      }

      GDALSetDescription(band, "band name here");
      GDALSetRasterNoDataValue(band, images[HIGHRES].meta.nodata);

    }
    
  }


  #pragma omp critical
  {
    GDALSetGeoTransform(file, images[HIGHRES].meta.transformation);
    GDALSetProjection(file, images[HIGHRES].meta.projection);
  }

  GDALClose(file);

  CSLDestroy(options);

  proctime_print("writing", TIME);

  return SUCCESS;
}
