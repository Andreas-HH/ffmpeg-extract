#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#include <stdio.h>
// #include "mpegvideo.h"

#define FEATURE_DIMENSION 100;

typedef int feature_elem;

typedef struct H264FeatureVector {
  feature_elem min;
  feature_elem max;
  int *mb_t;
  int *qp;
  feature_elem *v;
} H264FeatureVector;

typedef struct H264FeatureContext {
  int N;
  int ldc_range, lac_range, cdc_range, cac_range, qp_range;
  FILE *file; 
  H264FeatureVector *vec;
} H264FeatureContext;

void addCounts(H264FeatureContext *feature_context, int *level, int *run_before, int total_coeff, int totalZeros, int mb_type, int mb_xy, int qp, int n);
void storeFeatures(H264FeatureContext *feature_context);
void refreshFeatures(H264FeatureContext *feature_context);

#endif /* AVCODEC_H264_EXTRACT */