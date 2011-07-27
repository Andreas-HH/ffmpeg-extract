#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#define FEATURE_DIMENSION 100;

typedef struct H264FeatureVector {
  double *v;
} H264FeatureVector;

typedef struct H264FeatureContext {
  int N;
  H264FeatureVector *vec;
} H264FeatureContext;

void addCounts(H264FeatureContext *featureVec, int level[16], int total_coeff);

#endif /* AVCODEC_H264_EXTRACT */