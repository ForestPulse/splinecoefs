/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions for I/O-ing images
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "image_io.h"


void read_image(char *path, bandlist_t *bands, image_t *image){


  GDALDatasetH  fp_dataset;
  
  copy_string(image->path, STRLEN, path);
  if ((fp_dataset = GDALOpen(path, GA_ReadOnly))== NULL){ 
    fprintf(stderr, "could not open %s\n", path); exit(FAILURE);}

  copy_string(image->proj, STRLEN, GDALGetProjectionRef(fp_dataset));
  GDALGetGeoTransform(fp_dataset, image->geotran);

  image->nx = GDALGetRasterXSize(fp_dataset);
  image->ny = GDALGetRasterYSize(fp_dataset);
  image->nc = image->nx*image->ny;

  image->nb = GDALGetRasterCount(fp_dataset);

  if (bands != NULL){
    if (bands->n < 1){
      fprintf(stderr, "no bands specified for %s\n", path); exit(FAILURE);}
    for (int b=0; b<bands->n; b++){
      if (bands->number[b] < 1 || bands->number[b] > image->nb){
        fprintf(stderr, "band number %d out of range for %s\n", bands->number[b], path); exit(FAILURE);}
    }
    image->nb = bands->n;
  } 

  alloc_2D((void***)&image->data, image->nb, image->nc, sizeof(short));

  for (int b=0; b<image->nb; b++){

    GDALRasterBandH band;
    if (bands != NULL){
      band = GDALGetRasterBand(fp_dataset, bands->number[b]);
    } else {
      band = GDALGetRasterBand(fp_dataset, b+1);
    }

    int has_nodata;
    image->nodata = (short)GDALGetRasterNoDataValue(band, &has_nodata);
    if (!has_nodata){
      fprintf(stderr, "%s has no nodata value.\n", path); 
      exit(FAILURE);
    }

    if (GDALRasterIO(band, GF_Read, 0, 0, 
      image->nx, image->ny, image->data[b], 
      image->nx, image->ny, GDT_Int16, 0, 0) == CE_Failure){
      printf("could not read band %d from %s.\n", b+1, path); exit(FAILURE);}
  }

  GDALClose(fp_dataset);

  return;
}

void copy_image(image_t *from, image_t *to, int nbands, short nodata, char *path){

  copy_string(to->path, STRLEN, path);
  copy_string(to->proj, STRLEN, from->proj);
  for (int i=0; i<6; i++) to->geotran[i] = from->geotran[i];
  to->nx = from->nx;
  to->ny = from->ny;
  to->nc = from->nc;
  to->nb = nbands;
  to->nodata = nodata;

  alloc_2D((void***)&to->data, to->nb, to->nc, sizeof(short));

  return;
}

void write_image(image_t *image){


  GDALDriverH driver = NULL;
  if ((driver = GDALGetDriverByName("GTiff")) == NULL){
    printf("%s driver not found\n", "GTiff"); exit(FAILURE);}
    
  char **options = NULL;
  options = CSLSetNameValue(options, "COMPRESS", "ZSTD");
  options = CSLSetNameValue(options, "PREDICTOR", "2");
  options = CSLSetNameValue(options, "INTERLEAVE", "BAND");
  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "TILED", "YES");
  options = CSLSetNameValue(options, "BLOCKXSIZE", "256");
  options = CSLSetNameValue(options, "BLOCKYSIZE", "256");
  
  GDALDatasetH fp_dataset = NULL;
  if ((fp_dataset = GDALCreate(driver, image->path, image->nx, image->ny, image->nb, GDT_Int16, options)) == NULL){
    printf("Error creating file %s.\n", image->path); exit(FAILURE);}

  for (int b=0; b<image->nb; b++){

    GDALRasterBandH band = GDALGetRasterBand(fp_dataset, b+1);
    if (GDALRasterIO(band, GF_Write, 0, 0, 
      image->nx, image->ny, image->data[b], 
      image->nx, image->ny, GDT_Int16, 0, 0) == CE_Failure){
      printf("Unable to write band %d to %s.\n", b+1, image->path); 
      exit(FAILURE);
    }

    GDALSetRasterNoDataValue(band, image->nodata);

  }

  GDALSetGeoTransform(fp_dataset, image->geotran);
  GDALSetProjection(fp_dataset,   image->proj);

  GDALClose(fp_dataset);

  if (options != NULL) CSLDestroy(options);   

  return;
}

void free_image(image_t *image){
  free_2D((void**)image->data, image->nb);
  image->data = NULL;
  return;
}

void compare_images(image_t *image_1, image_t *image_2){
bool equal = true;

  if (image_1->nx != image_2->nx){
    fprintf(stderr, "Image dimensions nx do not match: %d vs %d\n", image_1->nx, image_2->nx);
    equal = false;
  }

  if (image_1->ny != image_2->ny){
    fprintf(stderr, "Image dimensions ny do not match: %d vs %d\n", image_1->ny, image_2->ny);
    equal = false;
  }

  if (image_1->nc != image_2->nc){
    fprintf(stderr, "Image dimensions nc do not match: %d vs %d\n", image_1->nc, image_2->nc);
    equal = false;
  }

  if (strcmp(image_1->proj, image_2->proj) != 0){
    fprintf(stderr, "Image projections do not match: %s vs %s\n", image_1->proj, image_2->proj);
    equal = false;
  }

  for (int i=0; i<6; i++){
    if (fabs(image_1->geotran[i] - image_2->geotran[i]) > 1e-6){
      fprintf(stderr, "Image geotransform parameters do not match at index %d: %f vs %f\n", i, image_1->geotran[i], image_2->geotran[i]);
      equal = false;
    }
  }

  if (equal == false){
    fprintf(stderr, "Images %s and %s are not compatible.\n", image_1->path, image_2->path);
    exit(FAILURE);
  }

  return;
}


