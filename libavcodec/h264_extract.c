#include "h264_extract.h"

/*
 (Luma DC:   48)
 Luma AC:   0 - 15
 Chroma DC: 49, 50
 Chroma AC: 16 - 19, 32 - 35
 */


void myprint(char *text) {
   printf("%s", text);
}

// When it comes to other embedding methods, be sure to change the method in header (init_rate_bins and init_features)
void simulate_hiding_plusminus(H264FeatureContext *fc) {
  int i;
  int *coefs = fc->tape;
  double r;
  int sl = fc->slice_type;
  
  if (sl == TYPE_I_SLICE) return;
  if (!(fc->accept_blocks & (1 << fc->blocknum))) {
    return;
  }

  for (i = 0; i < num_coefs[fc->blocknum]; i++) {
    if (coefs[i]<2 && coefs[i]>-2) continue;
    r = (((double) rand()) / ((double) RAND_MAX));
    if (r < fc->p_hide) {
      switch (fc->slice_type) {
	case TYPE_P_SLICE:
	  fc->hidden_bits_p++;
	  break;
	case TYPE_B_SLICE:
	  fc->hidden_bits_b++;
	  break;
      }
      if (coefs[i] == 2) {
	coefs[i]++;
      } else if (coefs[i] == -2) {
	coefs[i]--;
      } else{
	r = (((double) rand()) / ((double) RAND_MAX));
	if (r < PROB_INCREASE) {  // if changing, increase or decrease?
	  coefs[i] += 1;
	} else {
	  coefs[i] -= 1;
	}
      }
    }
  }
}

int get_block_index(int n) {
  int r = -1;
  
  if (n >= 0  && n <= 15)                             r = 0; // Luma
  if (n == 49 || n == 50)                             r = 1; // Chroma DC
  if ((n >= 16 && n <= 19) || (n >= 32 && n <= 35))   r = 2; // Chroma AC
  
  return r;
}

int get_rate_index(double rate) {
  int idx = (int) ((rate/MAX_RATE)*((double) NUM_BINS) + 0.5);
  if (idx >= NUM_BINS) return -1;
  return idx;
}

void constructProperCoefArray(int *result, int *level, int *run_before, int total_coeff, int totalZeros, int blocknum, H264FeatureVector *vec) {
  int i;
  int pos;
  int zerosToSee;
  
  for (i = 0; i < 16; i++)
    result[i] = 0;
  
  pos = total_coeff + totalZeros - 1;
  zerosToSee = totalZeros;
  result[pos] = level[0];
  pos--;
  for (i = 1; i < total_coeff; i++) {
    if (zerosToSee > 0) {
      zerosToSee -= run_before[i];
      pos -= run_before[i];
    }
    result[pos] = level[i];
    pos--;
  }
}

void addCounts(H264FeatureContext *fc, int qp, int n, int len) {
  int i, l, r;
  int coef_index;
  int *tape = fc->tape;
  int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  int sl = fc->slice_type;
  
  if (blocknum == -1 || qp_index < 0 || qp_index >= QP_RANGE)
    return;

//   if (sl != TYPE_P_SLICE) return;  // add only P-Slices 
  if (sl == TYPE_I_SLICE) return;
  // probably want to count bpnc later!
  
    for (i = 0; i < len; i++) { // num_coefs[blocknum]
      coef_index = tape[i];
      if (coef_index <= 0)  {
	coef_index = coef_index + ranges[blocknum][i]; // fc->histogram_
	if (coef_index < 0) continue;
      }
      else if (coef_index > 0)  {
	coef_index = coef_index + ranges[blocknum][i]; // -1
	if (coef_index > 2*ranges[blocknum][i]) continue; // -1
      }
      fc->vec->histograms[sl][qp_index][blocknum][i][coef_index]++;
    }
  //pairs
    for (i = 0; i < len-1; i++) { // num_coefs[blocknum]
      l = tape[i] + ranges[blocknum][0];                // we are allowed to exceed the local range here
      r = tape[i+1] + ranges[blocknum][1];              // space is limited by ranges of first two coefs
      if (l < 0 || l > 2*ranges[blocknum][0]) continue;
      if (r < 0 || r > 2*ranges[blocknum][1]) continue;
      fc->vec->pairs[sl][qp_index][blocknum][l][r]++;
    }
  // UvsV
  if (n == 49) {
    memcpy(fc->lastUs[0], tape, num_coefs[1]*sizeof(int));
    fc->seenUs[0] = 1; // it could happen that we skip one frame and align perfectly, 
  }                // VERY unlikely though. Comparing x,y with ux,uy should be enough
  else if (n == 50 && fc->seenUs[0] && fc->ux == fc->x && fc->uy == fc->y) {
    for (i = 0; i < num_coefs[1]; i++) {
      l = fc->lastUs[0][i] + ranges[1][0];
      r = tape[i+1] + ranges[1][0];
      if (l < 0 || l > 2*ranges[1][0]) continue;
      if (r < 0 || r > 2*ranges[1][0]) continue;
      fc->vec->uvsv[sl][qp_index][0][l][r]++;
    }
    fc->seenUs[0] = 0;
  }
  if (n >= 16 && n <= 19) {
    memcpy(fc->lastUs[n-15], tape, num_coefs[2]*sizeof(int));
    fc->seenUs[n-15] = 1;
  } else if ((n >= 32 && n <= 35) && fc->seenUs[n-31] && fc->ux == fc->x && fc->uy == fc->y) {
    for (i = 0; i < num_coefs[2]; i++) {
      l = fc->lastUs[n-31][i] + ranges[2][i];
      r = tape[i+1] + ranges[2][i];
      if (l < 0 || l > 2*ranges[2][i]) continue;
      if (r < 0 || r > 2*ranges[2][i]) continue;
      fc->vec->uvsv[sl][qp_index][i+1][l][r]++;
    }
    fc->seenUs[n-31] = 0;
  }
}

void storeFeatureVectors(H264FeatureContext* fc) {
  int sl, i, j, k, l;
  double scale;
  uint64_t count_h, count_p, count_u;
  uint64_t N_h, N_p, N_u;

  if (!fc->refreshed) {
    return;
  }
  
  for (sl = 0; sl < 2; sl++) {  // DROPPING B-FRAMES !!!
    count_h = 0ll;
    count_p = 0ll;
    count_u = 0ll;
    N_h = 0ll;
    N_p = 0ll;
    N_u = 0ll;
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
	for (k = 0; k < num_coefs[j]; k++) {
	  for (l = 0; l < 2*ranges[j][k]+1; l++) {
	    fc->vec->vector_histograms[count_h] = (double) fc->vec->histograms[sl][i][j][k][l]; // (double) 
	    count_h++;
	    N_h += fc->vec->histograms[sl][i][j][k][l];
	  }
	}
        // pairs
	for (k = 0; k < 2*ranges[j][0]+1; k++) {
	  for (l = 0; l < 2*ranges[j][1]+1; l++) {
// 	    if (k == ranges[j][0] && l == ranges[j][1]) continue; // 00: &&; 0x x0: ||
  //            printf("(%i, %i)\n", count_h, count_p);
	    fc->vec->vector_pairs[count_p] = (double) fc->vec->pairs[sl][i][j][k][l]; // (double) 
	    count_p++;
	    N_p += fc->vec->pairs[sl][i][j][k][l];
	  }
	}
      }
      // UvsV
      for (k = 0; k < 2*ranges[1][0]+1; k++) {
	for (l = 0; l < 2*ranges[1][0]+1; l++) {
	  fc->vec->vector_uvsv[count_u] = (double) fc->vec->uvsv[sl][i][0][k][l];
	  count_u++;
	  N_u += fc->vec->uvsv[sl][i][0][k][l];
	}
      }
      for (j = 1; j < 16; j++) {
	if (ranges[2][j-1] == 0) break;
	for (k = 0; k < 2*ranges[2][j-1]+1; k++) {
	  for (l = 0; l < 2*ranges[2][j-1]+1; l++) {
	    fc->vec->vector_uvsv[count_u] = (double) fc->vec->uvsv[sl][i][j-1][k][l];
	    count_u++;
	    N_u += fc->vec->uvsv[sl][i][j-1][k][l];
	  }
        }
      }
    }
    
    scale = 1./(double) N_h;
    for (i = 0; i < fc->vec->vector_histograms_dim; i++) {
      fc->vec->vector_histograms[i] *= scale;
    }
    scale = 1./(double) N_p;
    for (i = 0; i < fc->vec->vector_pairs_dim; i++) {
      fc->vec->vector_pairs[i] *= scale;
    }
    if (N_u != 0ll) {
      scale = 1./(double) N_u;
      for (i = 0; i < fc->vec->vector_uvsv_dim; i++) {
	fc->vec->vector_uvsv[i] *= scale;
      }
//       printf("N_u != 0, count_u = %i \n", fc->vec->vector_uvsv_dim);
    }

    // N_u = 0 is fairly common, chroma channels don't necessarily align
    if (N_h == 0 || N_p == 0) continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!

    fwrite(fc->vec->vector_histograms, sizeof(double), fc->vec->vector_histograms_dim, fc->files_hist[sl]); // double
    fwrite(fc->vec->vector_pairs, sizeof(double), fc->vec->vector_pairs_dim, fc->files_hist[sl]); // double
    fwrite(fc->vec->vector_uvsv, sizeof(double), fc->vec->vector_uvsv_dim, fc->files_hist[sl]); // double
  }
  
  fc->refreshed = 0;
}


void refreshFeatures(H264FeatureContext* fc) {
  int i, j, k, l, sl;
  
  fc->hidden_bits_b = 0ull;
  fc->hidden_bits_p = 0ull;

  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
//           feature_context->vec->N[i][j][k] = 0;
          for (l = 0; l < 2*ranges[j][k]; l++) {
            fc->vec->histograms[sl][i][j][k][l] = 0ull;
          }
        }
        // pairs
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          for (l = 0; l < 2*ranges[j][1]+1; l++) {
            fc->vec->pairs[sl][i][j][k][l] = 0ull;
          }
        }
      }
      // UvsV
      for (k = 0; k < 2*ranges[1][0]+1; k++) {
        for (l = 0; l < 2*ranges[1][0]+1; l++) {
          fc->vec->uvsv[sl][i][0][k][l] = 0ull;
        }
      }
      for (j = 1; j < 16; j++) {
	for (k = 0; k < 2*ranges[2][j-1]+1; k++) {
          for (l = 0; l < 2*ranges[2][j-1]+1; l++) {
            fc->vec->uvsv[sl][i][j][k][l] = 0ull;
          }
        }
      }
    }
  }
  
  fc->refreshed = 1;
}

H264FeatureContext* init_features(char* method_name, int accept_blocks, double p_hide, FILE**** bins_h, FILE**** bins_p, int extract_rate) { 
  int i, j, k, sl;
  H264FeatureContext *fc;
  H264FeatureVector *fv;
  long *****v;
  long *****w;
  long *****u;
  int hist_dim = 0;
  int pair_dim = 0;
  int uvsv_dim = 0;
  int *tape;
  char method;
  char b_h_path[512];
  char b_p_path[512];
  char p_h_path[512];
  char p_p_path[512];
  char *blockstring = blockstrings[accept_blocks];//"L";
    
  v = malloc(2*sizeof(uint64_t****));
  w = malloc(2*sizeof(uint64_t****));
  u = malloc(2*sizeof(uint64_t****));
  for (sl = 0; sl < 2; sl++) {
    v[sl] = malloc(QP_RANGE*sizeof(uint64_t***));
    w[sl] = malloc(QP_RANGE*sizeof(uint64_t***));
    u[sl] = malloc(QP_RANGE*sizeof(uint64_t***));
    for (i = 0; i < QP_RANGE; i++) {
      v[sl][i] = malloc(3*sizeof(uint64_t**));   // 3 types of blocks (Luma + Chroma DC/AC)
      w[sl][i] = malloc(3*sizeof(uint64_t**));
      for (j = 0; j < 3; j++) {
        // setup histograms
        v[sl][i][j] = malloc(num_coefs[j]*sizeof(uint64_t*));
        for (k = 0; k < num_coefs[j]; k++) {
          v[sl][i][j][k] = malloc((2*ranges[j][k]+1)*sizeof(uint64_t));
        }
        // setup pairs
        w[sl][i][j] = malloc((2*ranges[j][0]+1)*sizeof(uint64_t*)); // +1 for the zero
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          w[sl][i][j][k] = malloc((2*ranges[j][1]+1)*sizeof(uint64_t));
        }
      }
      // setup uvsv
      u[sl][i] = malloc(16*sizeof(uint64_t**));
      u[sl][i][0] = malloc((2*ranges[1][0]+1)*sizeof(uint64_t*));
      for (k = 0; k < 2*ranges[1][0]+1; k++) {
        u[sl][i][0][k] = malloc((2*ranges[1][0]+1)*sizeof(uint64_t));
      }
      for (j = 1; j < 16; j++) { // there are 16 chroma ac coefs
        u[sl][i][j] = malloc((2*ranges[2][j-1]+1)*sizeof(uint64_t*));
        for (k = 0; k < 2*ranges[2][j-1]+1; k++) {
          u[sl][i][j][k] = malloc((2*ranges[2][j-1]+1)*sizeof(uint64_t));
        }
      }
    }
  }

  for (i = 0; i < 3; i++) { // scan through blocks
    for (j = 0; j < num_coefs[i]; j++) {
      hist_dim += 2*ranges[i][j] + 1;
    }
    pair_dim += (2*ranges[i][0]+1) * (2*ranges[i][1]+1); // (0,0) included, other wise - 1
  }
  uvsv_dim += (2*ranges[1][0]+1)*(2*ranges[1][0]+1);
  for (j = 1; j < 16; j++) {
    if (2*ranges[2][j-1] > 0) uvsv_dim += (2*ranges[2][j-1]+1)*(2*ranges[2][j-1]+1);
  }

  fc    =   malloc(sizeof(H264FeatureContext));
  fv    =   malloc(sizeof(H264FeatureVector));
  tape  =   malloc(16*sizeof(int));  
  
  fv->vector_histograms_dim = QP_RANGE*hist_dim;
  fv->vector_pairs_dim      = QP_RANGE*pair_dim;
  fv->vector_uvsv_dim       = QP_RANGE*uvsv_dim;
  fv->vector_histograms     = malloc(fv->vector_histograms_dim*sizeof(double));
  fv->vector_pairs          = malloc(fv->vector_pairs_dim*sizeof(double));
  fv->vector_uvsv           = malloc(fv->vector_uvsv_dim*sizeof(double));
  fv->histograms     = v;
  fv->pairs          = w;
  fv->uvsv           = u;
  fc->vec            = fv;
  fc->files_hist     = malloc(2*sizeof(FILE*));
  if (p_hide < 0) {
    fc->files_hist[0]  = fopen("clean/p_clean.fv", "a");
    fc->files_hist[1] = fopen("clean/b_clean.fv", "a");
  } else {
      sprintf(p_h_path, "%s/p_%s_%i.fv", method_name, blockstring, (int) (p_hide*1000.+0.5));
      sprintf(b_h_path, "%s/b_%s_%i.fv", method_name, blockstring, (int) (p_hide*1000.+0.5));
      fc->files_hist[0]  = fopen(p_h_path, "a");
      fc->files_hist[1]  = fopen(b_h_path, "a");
  }
  
  if (p_hide < 0) {
    method = 0;
  } else {
    method = METHOD;
  }
  writeHeader(fc->files_hist[0], 0, TYPE_P_SLICE, method, 0, p_hide, accept_blocks);
  writeHeader(fc->files_hist[1], 0, TYPE_B_SLICE, method, 0, p_hide, accept_blocks);

  fc->tape = tape;
  fc->accept_blocks = accept_blocks;
  fc->p_hide = p_hide;
  fc->refreshed   = 0;
  fc->seenUs = (int*) malloc(5*sizeof(int));
  fc->lastUs = (int**) malloc(5*sizeof(int*));
  fc->lastUs[0] = (int*) malloc(num_coefs[1]*sizeof(int));
  for (i = 1; i < 5; i++) {
    fc->lastUs[i] = (int*) malloc(num_coefs[2]*sizeof(int));
  }

  return fc;
}

void close_features(H264FeatureContext* fc) {
//   printf("closing features \n");
  int i, j, k, sl;
  storeFeatureVectors(fc);
  
  if (fc->files_hist[0] != NULL) {
    fclose(fc->files_hist[0]);
  }
  if (fc->files_hist[1] != NULL) {
    fclose(fc->files_hist[1]);
  }

  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < num_coefs[j]; k++) {
          free(fc->vec->histograms[sl][i][j][k]);
        }
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          free(fc->vec->pairs[sl][i][j][k]);
        }
        free(fc->vec->histograms[sl][i][j]);
        free(fc->vec->pairs[sl][i][j]);
      }
      // free uvsv
      for (k = 1; k < 2*ranges[1][0]+1; k++) {
        free(fc->vec->uvsv[sl][i][0][k]);
      }
      free(fc->vec->uvsv[sl][i][0]);
      for (j = 1; j < 16; j++) { // there are 16 chroma coefs
        for (k = 0; k < 2*ranges[2][j-1]+1; k++) {
          free(fc->vec->uvsv[sl][i][j][k]);
        }
        free(fc->vec->uvsv[sl][i][j]);
      }
      free(fc->vec->histograms[sl][i]);
      free(fc->vec->pairs[sl][i]);
      free(fc->vec->uvsv[sl][i]);
    }
    free(fc->vec->histograms[sl]);
    free(fc->vec->pairs[sl]);
    free(fc->vec->uvsv[sl]);
  }
  free(fc->vec->histograms);
  free(fc->vec->pairs);
  free(fc->vec->uvsv);
  free(fc->vec->vector_histograms);
  free(fc->vec->vector_pairs);
  
  for (i = 0; i < 5; i++) {
    free(fc->lastUs[i]);
  }
  free(fc->lastUs);
  free(fc->seenUs);
  
  free(fc->tape);
  free(fc->vec);
  free(fc);
}

void writeHeader(FILE *file, char pair, char slice_type, char method, char using_rate, double rate, char accept){
  if (file == NULL || ftell(file) != 0L) return;
  int i, j;
  char qp_offset = QP_OFFSET;
  char qp_range = QP_RANGE;
   
  fwrite(&pair, sizeof(char), 1, file);
  fwrite(&slice_type, sizeof(char), 1, file);
  fwrite(&method, sizeof(char), 1, file);
  if (method != 0) {  // don't write a rate/accept for clean features
    fwrite(&using_rate, sizeof(char), 1, file);
    fwrite(&rate, sizeof(double), 1, file);
    fwrite(&accept, sizeof(char), 1, file);
  }
  fwrite(&qp_offset, sizeof(char), 1, file);
  fwrite(&qp_range, sizeof(char), 1, file);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < num_coefs[i]; j++) {
      fwrite(&(ranges[i][j]), sizeof(unsigned char), 1, file);
    }
  }
}
