#include "h264_extract.h"

void addCounts(H264FeatureContext *fc, int *level, int total_coeff) {
  int i;
  
  for (i = 0; i < total_coeff; i++) {
    if (level[i] < fc->vec->min) fc->vec->min = level[i];
    if (level[i] > fc->vec->max) fc->vec->max = level[i];
    if (level[i] >= -50 && level[i] < 50) {
      fc->vec->v[level[i]+50] += 1;
    }
  }
  fc->N++;
}

void storeFeatures(H264FeatureContext *fc) {
  int i;
  
  fprintf(fc->file, "N = %i, min = %i, max = %i, v = ", fc->N, fc->vec->min, fc->vec->max);
  for (i = 0; i < 100; i++) {
    fprintf(fc->file, " %i ", fc->vec->v[i]);
  }
  fprintf(fc->file, "\n");
  fflush(fc->file);
}

void refreshFeatures(H264FeatureContext* feature_context) {
  feature_context->vec->min = 0;
  feature_context->vec->max = 0;
  feature_context->N = 0;
  for (int i = 0; i < 100; i++) {
    feature_context->vec->v[i] = 0;
  }
}


