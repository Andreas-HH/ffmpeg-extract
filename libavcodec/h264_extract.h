#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#include <stdio.h>

#define FEATURE_DIMENSION 100;

typedef int feature_elem;

typedef struct H264FeatureVector {
  feature_elem min;
  feature_elem max;
  feature_elem *v;
} H264FeatureVector;

typedef struct H264FeatureContext {
  int N;
  FILE *file; 
  H264FeatureVector *vec;
} H264FeatureContext;

void addCounts(H264FeatureContext *feature_context, int *level, int total_coeff);
void storeFeatures(H264FeatureContext *feature_context);
void refreshFeatures(H264FeatureContext *feature_context);

#endif /* AVCODEC_H264_EXTRACT */