#include "h264_extract.h"

void addCounts(H264FeatureContext *fc, int *level, int total_coeff) {
  int i = 0;
  int y = 0;
  
  for (i = 0; i < total_coeff; i++) {
    y = level[i];
    if (level[i] >= -50 && level[i] < 50) {
      fc->vec->v[level[i]+50] += 1.;
    }
  }
  fc->N++;
}
