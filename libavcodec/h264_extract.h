#ifndef AVCODEC_H264_EXTRACT
#define AVCODEC_H264_EXTRACT

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define METHOD                 44 // 0 = clean, 1=pm
#define QP_RANGE               1
#define QP_OFFSET              16
#define QP_DELTA               4
#define QP_JUMPS               4
#define TYPE_I_SLICE           2
#define TYPE_P_SLICE           0
#define TYPE_B_SLICE           1
#define PROB_INCREASE          0.5
#define ACCEPT_L               1
#define ACCEPT_LC              7
#define ACCEPT_C               6
// #define MAX_RATE               0.24
// #define NUM_BINS               12
#define FRAC                   0.4
// #define FRAC_C                 1.
#define PROB_DELTA             0.001
#define STEGF                  10
#define THRESHOLD              0
#define MIN_COEF               1   // this counts from 1, coef 1 is DC coef
// #define ACCEPT                 1

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )

typedef int feature_elem;
typedef uint32_t store_elem;

static const char *blockstrings[8]        = {"clean", "L", "C_dc", "LC_dc", "C_ac", "LC_ac", "C", "LC"};
static const char num_coefs[3]            =  {16, 4, 15}; // Luma, Cr DC, Cb DC, Cr AC, Cb AC
static const unsigned char ranges[3][16]  =  {{15,  11, 10,  8, 8, 8,  5, 5, 3, 3,  0, 0, 0,  0, 0,  0}, 
                                              {8,  6, 6,  4,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					      {4,  4,  3, 3, 3,  0, 0, 0, 0,  0, 0, 0,  0, 0,  0,  0}};

typedef struct H264FeatureVector {
  int vector_num;
  store_elem *****histograms;          // [slice_type][qp][block][coef][element]
  store_elem *****pairs;               // [slice_type][qp][block][element_left][element_right]
  store_elem *****uvsv;                // [slice_type][qp][coef][element_u][element_v]  | coef=0 is DC, coef in [1..16) AC
  
  int vector_histograms_dim;
  int vector_pairs_dim;
  int vector_uvsv_dim;
  store_elem *vector_histograms;
  store_elem *vector_pairs;
  store_elem *vector_uvsv;
} H264FeatureVector;

typedef struct H264FeatureContext {
  int qp;
  int *tape;
  int slice_type;
  int refreshed;
  uint64_t hidden_bits_b;
  uint64_t hidden_bits_p;
  uint64_t num_coefs_b;
  uint64_t num_coefs_p;
  uint64_t num_vectors_b;
  uint64_t num_vectors_p;
  double bpnc_b;
  double bpnc_p;
  char *logName;
  
  // simulation parameters:
  int current_qp;
  int *proper_coefs;
  int accept_blocks;  // bit at position blocknum tells if that block type is accepted or not  (accept & (1 << blocknum))
  double p_hide;

  int *seenUs; // [0] = DC, otheriwse = AC, ranges from 0 to 4
  int **lastUs;
  int ux, uy, x, y;
  FILE *logfile;
  FILE **files_hist;  // [SLICE_TYPE]
  H264FeatureVector *vec;
} H264FeatureContext;

void myprint(char *text);
void myRandom(double* r);

H264FeatureContext* init_features(char* method_name, int accept_blocks, double p_hide, int qp);
void close_features(H264FeatureContext *fc);
void writeHeader(FILE *file, char pair, char slice_type, char method, double prob, char accept, int qp);
int get_block_index(int n);
void simulate_hiding_plusminus(H264FeatureContext* fc, int blocknum, int thresh);
void simulate_hiding_f5(H264FeatureContext* fc, int blocknum, int thresh);
void simulate_hiding_lsb(H264FeatureContext* fc, int blocknum, int thresh);
void constructProperCoefArray(int* result, int* level, int* run_before, int total_coeff, int totalZeros);
void addCounts(H264FeatureContext* fc, int qp, int n, int len);
void storeFeatureVectors(H264FeatureContext *feature_context);
void refreshFeatures(H264FeatureContext *feature_context);

#endif /* AVCODEC_H264_EXTRACT */