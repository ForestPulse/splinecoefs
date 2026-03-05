/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions for screening quality bit files
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "quality.h"

enum { _QAI_BIT_OFF_ =  0, _QAI_BIT_CLD_ =  1, _QAI_BIT_SHD_ =  3, 
       _QAI_BIT_SNW_ =  4, _QAI_BIT_WTR_ =  5, _QAI_BIT_AOD_ =  6, 
       _QAI_BIT_SUB_ =  8, _QAI_BIT_SAT_ =  9, _QAI_BIT_SUN_ = 10, 
       _QAI_BIT_ILL_ = 11, _QAI_BIT_SLP_ = 13, _QAI_BIT_WVP_ = 14};

short get_qai_from_value(short value, int index, int bitfields);
bool get_off_from_value(short value);
char get_cloud_from_value(short value);
bool get_shadow_from_value(short value);
bool get_snow_from_value(short value);
bool get_water_from_value(short value);
char get_aerosol_from_value(short value);
bool get_subzero_from_value(short value);
bool get_saturation_from_value(short value);
bool get_lowsun_from_value(short value);
char get_illumination_from_value(short value);
bool get_slope_from_value(short value);
bool get_vaporfill_from_value(short value);

/** Decide whether to use this pixel
+++ This function checks the QAI layer against the user-defined QAI crite-
+++ ria.
--- qai:      Quality Assurance Information
--- p:        pixel
--- qai_rule: ruleset for QAI filtering
+++ Return:   true/false
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
bool use_this_pixel(short qai){

  if (get_off_from_value(qai))               return false; // on/off flag
  if (get_cloud_from_value(qai) == 1)        return false; // cloud uncertain
  if (get_cloud_from_value(qai) == 2)        return false; // cloud opaque
  if (get_cloud_from_value(qai) == 3)        return false;  // cloud cirrus
  if (get_shadow_from_value(qai))            return false; // cloud shadow
  if (get_snow_from_value(qai))              return false; // snow
  //if (get_water_from_value(qai))             return false; // water
  //if (get_aerosol_from_value(qai) == 1)      return false; // aerosol interpolated
  //if (get_aerosol_from_value(qai) == 2)      return false; // aerosol high
  //if (get_aerosol_from_value(qai) == 3)      return false; // aerosol fill
  if (get_subzero_from_value(qai))           return false; // subzero reflectance
  if (get_saturation_from_value(qai))        return false; // saturated reflectance
  //if (get_lowsun_from_value(qai))            return false; // low sun angle
  //if (get_illumination_from_value(qai) == 1) return false; // low illumination
  //if (get_illumination_from_value(qai) == 2) return false; // poor illumination
  if (get_illumination_from_value(qai) == 3) return false; // shadow illumination
  //if (get_slope_from_value(qai))             return false; // high slope
  //if (get_vaporfill_from_value(qai))         return false; // water vapor fill

  return true;
}



/** This function reads any quality bit in the QAI layer
+++ The same as get_qai, but reads the value from a short value directly
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
short get_qai_from_value(short value, int index, int bitfields){
int i;
short val = 0;

  for (i=0; i<bitfields; i++) val |= (short)(1 << i);

  return (short)(value >> index) & val;
}


/** read off/on flag, directly from value **/
bool get_off_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_OFF_, 1);
}

/** read cloud flag, directly from value **/
char get_cloud_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_CLD_, 2);
}

/** read cloud shadow flag, directly from value **/
bool get_shadow_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SHD_, 1);
}

/** read snow flag, directly from value **/
bool get_snow_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SNW_, 1);
}

/** read water flag, directly from value **/
bool get_water_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_WTR_, 1);
}

/** read aerosol flag, directly from value **/
char get_aerosol_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_AOD_, 2);
}

/** read subzero reflectance flag, directly from value **/
bool get_subzero_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SUB_, 1);
}

/** read saturated reflectance flag, directly from value **/
bool get_saturation_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SAT_, 1);
}

/** read low sun angle flag, directly from value **/
bool get_lowsun_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SUN_, 1);
}

/** read illumination flag, directly from value **/
char get_illumination_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_ILL_, 2);
}

/** read slope flag, directly from value **/
bool get_slope_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SLP_, 1);
}

/** read water vapor fill flag, directly from value **/
bool get_vaporfill_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_WVP_, 1);
}
