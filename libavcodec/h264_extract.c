#include "h264_extract.h"

void addCounts(H264FeatureContext *fc, int level[16], int total_coeff) {
  int i;
  
  for (i = 0; i < total_coeff; i++) {
    if (level[i] >= -50 && level[i] < 50) {
      fc->vec->v[level[i]+50] += 1.;
    }
  }
  fc->N++;
}