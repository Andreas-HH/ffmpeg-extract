#include "h264_extract.h"

/*
 (Luma DC:   48)
 Luma AC:   0 - 15
 Chroma DC: 49, 50
 Chroma AC: 16 - 19, 32 - 35
 */

void simulate_hiding_plusminus(int* coefs, int size) {
  int i;
  double r;
  
  for (i = 0; i < size; i++) {
     if (coefs[i]<2) continue;
    r = (((double) rand()) / ((double) RAND_MAX));
//     printf("%f \n ", r);
    if (r > PROB_SKIP) {
      r = (((double) rand()) / ((double) RAND_MAX));
      if (r < PROB_KEEP) {
	r = (((double) rand()) / ((double) RAND_MAX));
	if (r < PROB_INCREASE) {
	  coefs[i] += 1;
	} else {
	  coefs[i] -= 1;
	}
      }
    }
  }
}


void setup_ranges(int **ranges, int luma, int chroma_dc, int chroma_ac) {
  int i;
  
  for (i = 0; i < 16; i++) {
    if (luma_ranges[i] > 2)
      ranges[0][i] = luma_ranges[i];
    else
      ranges[0][i] = 0;
  }
  for (i = 0; i < 4; i++) {
    if (chroma_dc_ranges[i] > 2)
      ranges[1][i] = chroma_dc_ranges[i];
    else
      ranges[1][i] = 0;
  }
  for (i = 0; i < 15; i++) {
    if (chroma_ac_ranges[i] > 2)
      ranges[2][i] = chroma_ac_ranges[i];
    else
      ranges[2][i] = 0;
  }
}

int get_block_index(int n) {
  int r = -1;
  
  if (n >= 0  && n <= 15)                             r = 0; // Luma
  if (n == 49 || n == 50)                             r = 1; // Chroma DC
  if ((n >= 16 && n <= 19) || (n >= 32 && n <= 35))   r = 2; // Chroma AC
  
  return r;
}

void constructProperCoefArray(int *result, int *level, int *run_before, int total_coeff, int totalZeros, int blocknum, H264FeatureVector *vec) {
  int i;
  int pos;
  int zerosToSee;
  
  for (i = 0; i < 16; i++)
    result[i] = 0;
  
  if (blocknum == 2) // Chroma AC
    pos = total_coeff + totalZeros;
  else
    pos = total_coeff + totalZeros - 1;
  zerosToSee = totalZeros;
  result[pos] = level[0];
  pos--;
  for (i = 1; i < total_coeff; i++) {
    if (zerosToSee > 0) {
      zerosToSee -= run_before[i];
      pos -= run_before[i];
    }
    if (level[i] < vec->min) vec->min = level[i];
    if (level[i] > vec->max) vec->max = level[i];
    result[pos] = level[i];
    pos--;
  }
}

void addCounts(H264FeatureContext *fc, int qp, int n) {
  int i, l, r;
  int coef_index;
  int *tape = fc->tape;
  int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  int sl = fc->slice_type;
  
//   if (sl != TYPE_P_SLICE) return;  // add only P-Slices 
  if (sl == TYPE_I_SLICE) return;
  
  fc->vec->qp[qp]++;
  if (blocknum == -1) {
//     fprintf(fc->file, "blocknum is -1 o_O %i \n", n);
    return;
  }
  if (qp_index < 0 || qp_index >= QP_RANGE) {
//     fprintf(fc->file, "qp out of range o_O %i \n", qp);
    return;
  }
  
//   constructProperCoefArray(tape, level, run_before, total_coeff, totalZeros, blocknum, fc->vec);
  for (i = 0; i < num_coefs[blocknum]; i++) {
    coef_index = tape[i];
    if (coef_index == 0) continue;
    else if (coef_index < 0)  {
      coef_index = coef_index + fc->histogram_ranges[blocknum][i];
      if (coef_index < 0) continue;
    }
    else if (coef_index > 0)  {
      coef_index = coef_index + fc->histogram_ranges[blocknum][i]-1;
      if (coef_index > 2*fc->histogram_ranges[blocknum][i]-1) continue;
    }
//     fc->vec->N[qp_index][blocknum][i]++;
    fc->vec->histograms[sl][qp_index][blocknum][i][coef_index]++;
  }
  
  //pairs
  for (i = 0; i < num_coefs[blocknum]-1; i++) {
    l = tape[i] + fc->histogram_ranges[blocknum][0];                // we are allowed to exceed the local range here
    r = tape[i+1] + fc->histogram_ranges[blocknum][1];              // space is limited by ranges of first two coefs
    if (l < 0 || l > 2*fc->histogram_ranges[blocknum][0]) continue;
    if (r < 0 || r > 2*fc->histogram_ranges[blocknum][1]) continue;
    fc->vec->pairs[sl][qp_index][blocknum][l][r]++;
  }
}

void storeCounts(H264FeatureContext *fc) {
  int i, j, k, l;
//   if (fc->slice_type != TYPE_P_SLICE) return;
  
  if (!fc->refreshed) return;  // avoid zero-vector at the beginning
  
  fprintf(fc->logfile, "I-Slices: %i, P-Slices: %i, B-Slices: %i    ", fc->i_slices, fc->p_slices, fc->b_slices);
  
  fprintf(fc->logfile, "min = %i, max = %i, vector_num = %i, slice_type = %i, histograms:\n", fc->vec->min, fc->vec->max, fc->vec->vector_num, fc->slice_type);
    // store qp histogram
  for (i = 0; i < 50; i++) {
    if (i == 20) 
      fprintf(fc->logfile, " [%u] ", fc->vec->qp[i]);
    else 
      fprintf(fc->logfile, " %u ", fc->vec->qp[i]);
  }
  fprintf(fc->logfile, "\n");
  
  for (i = 0; i < QP_RANGE; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < num_coefs[j]; k++) {
	// don't save zero hestograms
 	if (fc->histogram_ranges[j][k] == 0 || fc->vec->histograms[0][i][j][k][fc->histogram_ranges[j][k]] == 0) continue;
// 	if (fc->vec->N[i][j][k] == 0) continue;
	
	fprintf(fc->logfile, "qp = %i, blockn = %i, --- ", QP_OFFSET + i, j);//, fc->vec->N[i][j][k]);
	for (l = 0; l < fc->histogram_ranges[j][k]; l++) {
	  fprintf(fc->logfile, "%i ", fc->vec->histograms[0][i][j][k][l]);
	}
	fprintf(fc->logfile, "[0] ");
	for (l = fc->histogram_ranges[j][k]; l < 2*fc->histogram_ranges[j][k]; l++) {
	  fprintf(fc->logfile, "%i ", fc->vec->histograms[0][i][j][k][l]);
	}
	fprintf(fc->logfile, "\n");
      }
    }
  }
  
  fc->refreshed = 0;
}

void storeFeatureVectors(H264FeatureContext* fc) {
  int sl, i, j, k, l;
  int count_h, count_p;
  int N_h, N_p;
  long pos;
  
  if (!fc->refreshed) return;
  
  for (sl = 0; sl < 2; sl++) {
    count_h = 0;
    count_p = 0;
    N_h = 0;
    N_p = 0;
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
          for (l = 0; l < 2*fc->histogram_ranges[j][k]; l++) {
//             printf("(%i, %i, %i, %i) -- (%i, %i)\n", i, j, k, l, count_h, count_p);
             fc->vec->vector_histograms[count_h] = (double) fc->vec->histograms[sl][i][j][k][l];
            count_h++;
             N_h += fc->vec->pairs[sl][i][j][k][l];
//             fprintf(fc->files_hist[sl], "%i ", fc->vec->histograms[sl][i][j][k][l]);
          }
        }
        // pairs
        for (k = 0; k < 2*fc->histogram_ranges[j][0]+1; k++) {
          for (l = 0; l < 2*fc->histogram_ranges[j][1]+1; l++) {
            if (k == fc->histogram_ranges[j][0] && l == fc->histogram_ranges[j][1]) continue; // 00: &&; 0x x0: ||
//             printf("(%i, %i)\n", count_h, count_p);
             fc->vec->vector_pairs[count_p] = (double) fc->vec->pairs[sl][i][j][k][l];
            count_p++;
             N_p += fc->vec->pairs[sl][i][j][k][l];
//             fprintf(fc->files_pair[sl], "%i ", fc->vec->pairs[sl][i][j][k][l]);
          }
        }
      }
    }
//     printf("(%i, %i)\n", count_h, count_p);
    
    for (i = 0; i < fc->vec->vector_histograms_dim; i++)
      fc->vec->vector_histograms[i] *= ((double) fc->vec->vector_histograms_dim)/((double) N_h);
    for (i = 0; i < fc->vec->vector_pairs_dim; i++)
      fc->vec->vector_pairs[i] *= ((double) fc->vec->vector_pairs_dim)/((double) N_p);
  
    pos = ftell(fc->files_hist[sl]);
    if (pos == 0L)
      fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_hist[sl]);
    fwrite(fc->vec->vector_histograms, sizeof(double), fc->vec->vector_histograms_dim, fc->files_hist[sl]);
    pos = ftell(fc->files_pair[sl]);
    if (pos == 0L)
      fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, fc->files_pair[sl]);
    fwrite(fc->vec->vector_pairs, sizeof(double), fc->vec->vector_pairs_dim, fc->files_pair[sl]);
//     fprintf(fc->files_hist[sl], "\n");
//     fprintf(fc->files_pair[sl], "\n");
  }
//   printf("fertig!\n");
  
  fc->refreshed = 0;
}


void refreshFeatures(H264FeatureContext* feature_context) {
  int i, j, k, l, sl;
  
  feature_context->i_slices = 0;
  feature_context->p_slices = 0;
  feature_context->b_slices = 0;
  feature_context->vec->vector_num++;
  feature_context->vec->min = 0;
  feature_context->vec->max = 0;

  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
//           feature_context->vec->N[i][j][k] = 0;
          for (l = 0; l < 2*feature_context->histogram_ranges[j][k]; l++) {
            feature_context->vec->histograms[sl][i][j][k][l] = 0;
          }
        }
        // pairs
        for (k = 0; k < 2*feature_context->histogram_ranges[j][0]; k++) {
          for (l = 0; l < 2*feature_context->histogram_ranges[j][1]; l++) {
            feature_context->vec->pairs[sl][i][j][k][l] = 0;
          }
        }
      }
    }
  }
//   for (i = 0; i < 396; i++) {
//     feature_context->vec->mb_t[i] = 0;
//   }
  for (i = 0; i < 50; i++) {
    feature_context->vec->qp[i] = 0;
  }
  
  feature_context->refreshed = 1;
}


