#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#include <stdlib.h>
#include <stdio.h>
// #include "mpegvideo.h"

#define FEATURE_DIMENSION    100
#define QP_RANGE              10
#define QP_OFFSET             14
#define USED_PIXELS           54
#define TYPE_I_SLICE           0
#define TYPE_P_SLICE           1
#define TYPE_B_SLICE           2
#define LUMA_RANGE            10
#define CHROMA_DC_RANGE        4
#define CHROMA_AC_RANGE        5
#define PROB_SKIP              0.
#define PROB_KEEP              0.5
#define PROB_INCREASE          0.5

typedef int feature_elem;

static const int num_coefs[3] = {16, 4, 15}; // Luma, Cr DC, Cb DC, Cr AC, Cb AC
static const int luma_ranges[16] =      {14,  10, 10,  8, 8, 8,  5, 5, 5, 5,  3, 3, 3,  2, 2,  1};
static const int chroma_dc_ranges[4] =  {8,  6, 6,  4};
static const int chroma_ac_ranges[15] = {4,  4,  3, 3, 3,  2, 2, 2, 2,  1, 1, 1,  1, 1,  1};

/*static const int histogram_range[5][16] = { 
  {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20},
  {10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0},
  {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0}
};//*/

typedef struct H264FeatureVector {
//   int num_blocks;
  int slice_type;
  int vector_num;
  feature_elem min;
  feature_elem max;
  int *mb_t;
  int *qp;
  int ***N;                     //[qp][block][coef]
  feature_elem ****histogram;   //[qp][block][coef][element]
} H264FeatureVector;

typedef struct H264FeatureContext {
//   int ldc_range, lac_range, cdc_range, cac_range, qp_range;
  int refreshed;
  int *tape;
  int i_slices;
  int p_slices;
  int b_slices;
  int **histogram_ranges;  // entry n specifies range -n..n, excluding 0
  FILE *file;
  FILE *p_hist;
  FILE *p_pair;
  FILE *b_hist;
  FILE *b_pair;
  H264FeatureVector *vec;
} H264FeatureContext;

void   simulate_hiding_plusminus(int *coefs, int size);
// double getRandNr();
void   setup_ranges(int **ranges, int luma, int chroma_dc, int chroma_ac);
int    get_block_index(int n);
void   constructProperCoefArray(int *result, int *level, int *run_before, int total_coeff, int totalZeros, int blocknum,  H264FeatureVector *vec);
void   addCounts(H264FeatureContext *feature_context, int *level, int *run_before, int total_coeff, int totalZeros, int mb_type, int mb_xy, int qp, int n);
void   storeCounts(H264FeatureContext *feature_context);
void   refreshFeatures(H264FeatureContext *feature_context);

#endif /* AVCODEC_H264_EXTRACT */