#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#define QP_RANGE              20
#define QP_OFFSET             10
#define TYPE_I_SLICE           2
#define TYPE_P_SLICE           0
#define TYPE_B_SLICE           1
#define PROB_KEEP              0.5
#define PROB_INCREASE          0.5
#define ACCEPT_L               1
#define ACCEPT_LC              7
#define ACCEPT_C               6
#define MAX_RATE               0.1
#define NUM_BINS               50

typedef int feature_elem;

// static const int MAX_RATE = 0.1;
// static const int NUM_BINS = 10;
static const char *blockstrings[8] = {"clean", "L", "C_dc", "LC_dc", "C_ac", "LC_ac", "C", "LC"};
static const int num_coefs[3]         =  {16, 4, 15}; // Luma, Cr DC, Cb DC, Cr AC, Cb AC
// static const int luma_ranges[16]      =  {14,  10, 10,  8, 8, 8,  5, 5, 5, 5,  3, 3, 3,  2, 2,  1};
// static const int chroma_dc_ranges[4]  =  {8,  6, 6,  4};
// static const int chroma_ac_ranges[15] =  {4,  4,  3, 3, 3,  2, 2, 2, 2,  1, 1, 1,  1, 1,  1};
static const int ranges[3][16]        =  {{14,  10, 10,  8, 8, 8,  5, 5, 5, 5,  3, 3, 3,  0, 0,  0}, 
                                          {8,  6, 6,  4,   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, 
					  {4,  4,  3, 3, 3,  0, 0, 0, 0,  0, 0, 0,  0, 0,  0,   -1}};

typedef struct H264FeatureVector {
  int vector_num;
  int min;
  int max;
//   int *mb_t;
  int *qp;
//   int ***N;                     // [slice_type][qp][block][coef]
  int *****histograms;          // [slice_type][qp][block][coef][element]
  int *****pairs;               // [slice_type][qp][block][element_left][element_right]
  
  int vector_histograms_dim;
  int vector_pairs_dim;
  double *vector_histograms;
  double *vector_pairs;
} H264FeatureVector;

typedef struct H264FeatureContext {
//   int ldc_range, lac_range, cdc_range, cac_range, qp_range;
  int *tape;
  int extract_rate;
  
  // vector stats:
  int slice_type;
  int refreshed;
  int i_slices;
  int p_slices;
  int b_slices;
  int hidden_bits_b;
  int hidden_bits_p;
  int num_4x4_blocks_b;
  int num_4x4_blocks_p;
  
  // simulation parameters:
  int current_qp;
  int *proper_coefs;
  int blocknum;
  int accept_blocks;  // bit at position blocknum tells if that block type is accepted or not  (accept & (1 << blocknum))
  double p_hide;
  int locked_by_init;
//   pthread_mutex_t *thread_mutex;
//   pthread_mutex_t *main_mutex;
//   pthread_attr_t *thread_attr;
//   pthread_t *thread;
  
//   int **histogram_ranges;  // entry n specifies range -n..n, excluding 0 ::: [block][coef]
  FILE *logfile;
  FILE **files_hist;  // [SLICE_TYPE]
  FILE **files_pair;
  FILE ****rate_bins_hist; // [accept][sl][rate_idx]
  FILE ****rate_bins_pair; // [accept][sl][rate_idx]
//   FILE *b_hist;
//   FILE *b_pair;
  H264FeatureVector *vec;
} H264FeatureContext;

void myprint(char *text);

H264FeatureContext* init_features(char* method_name, int accept_blocks, double p_hide, FILE**** bins_h, FILE**** bins_p, int extract_rate);
void close_features(H264FeatureContext *fc);
FILE**** init_rate_bins(char* type, char* method_name);
void close_rate_bins(FILE ****bins);
int get_block_index(int n);
int get_rate_index(double rate);
// void perform_hiding_plusminus(void *vfc);
void simulate_hiding_plusminus(H264FeatureContext *fc);
// void wait_for_simulation(H264FeatureContext *fc);
// void setup_ranges(int **ranges, int luma, int chroma_dc, int chroma_ac);
void constructProperCoefArray(int *result, int *level, int *run_before, int total_coeff, int totalZeros, int blocknum,  H264FeatureVector *vec);
void addCounts(H264FeatureContext* fc, int qp, int blocknum);
void storeCounts(H264FeatureContext *feature_context);
void storeFeatureVectors(H264FeatureContext *feature_context);
void refreshFeatures(H264FeatureContext *feature_context);

#endif /* AVCODEC_H264_EXTRACT */