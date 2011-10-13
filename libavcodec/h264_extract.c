#include "h264_extract.h"

void addCounts(H264FeatureContext *fc, int *level, int *run_before, int total_coeff, int totalZeros, int mb_type, int mb_xy, int qp, int n) {
  int i;
  
  if (mb_type&0x0002)
    fprintf(fc->file, "%i: is intra 16x16! ", mb_xy);
  else
    fprintf(fc->file, "%i:                 ", mb_xy);
  
  fprintf(fc->file, "n = %i, total_coeff = %i, totalZeros = %i", n, total_coeff, totalZeros);
  for (i = 0; i < total_coeff; i++) {
    if (level[i] < fc->vec->min) fc->vec->min = level[i];
    if (level[i] > fc->vec->max) fc->vec->max = level[i];
    if (level[i] >= -50 && level[i] < 50) {
      fc->vec->v[level[i]+50] += 1;
    }
    fprintf(fc->file, "(%i, %i)", level[i], run_before[i]);
  }
  fprintf(fc->file, "\n");
  fc->N++;
  
  if (mb_xy < 396) {
    fc->vec->mb_t[mb_xy] = mb_type;
  } else {
    printf("Etwas stimmt hiuer nicht! %d \n", mb_xy);
  }
  
  fc->vec->qp[qp]++;
}

void storeFeatures(H264FeatureContext *fc) {
  int i;
  
  fprintf(fc->file, "N = %i, min = %i, max = %i, v = ", fc->N, fc->vec->min, fc->vec->max);
  for (i = 0; i < 100; i++) {
    fprintf(fc->file, " %i ", fc->vec->v[i]);
  }
  fprintf(fc->file, "\n");
  
  // store mb types
  for (i = 0; i < 396; i++) {
    fprintf(fc->file, " %u ", fc->vec->mb_t[i]);
  }//*/
  fprintf(fc->file, "\n");
  
  // store qp histogram
  for (i = 0; i < 50; i++) {
    if (i == 20) 
      fprintf(fc->file, " [%u] ", fc->vec->qp[i]);
    else 
      fprintf(fc->file, " %u ", fc->vec->qp[i]);
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
  for (int i = 0; i < 396; i++) {
    feature_context->vec->mb_t[i] = 0;
  }
  for (int i = 0; i < 50; i++) {
    feature_context->vec->qp[i] = 0;
  }
}


