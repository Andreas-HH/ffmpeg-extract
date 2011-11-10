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
  
//   for (i = 0; i < luma; i++) {
//      ranges[0][i] = 30;
// //     ranges[0][i] = luma - i;
//   }
//   for (; i < 16; i++) {
//     ranges[0][i] = 30; // 0
//   }
//   
//   for (i = 0; i < 4; i++) {
// //     ranges[1][i] = chroma_dc;
//     ranges[1][i] = 20;
//   }
//   
//   for (i = 0; i < luma; i++) {
// //     ranges[2][i] = chroma_ac - i;
//     ranges[2][i] = 20;
//   }
//   for (; i < 15; i++) {
//     ranges[2][i] = 20; // 0
//   }
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

void addCounts(H264FeatureContext *fc, int *level, int *run_before, int total_coeff, int totalZeros, int mb_type, int mb_xy, int qp, int n) {
  int i;
  int coef_index;
  int *tape = fc->tape;
  int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  
  if (fc->vec->slice_type != TYPE_P_SLICE) return;  // add only P-Slices 
  
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
    fc->vec->N[qp_index][blocknum][i]++;
    fc->vec->histogram[qp_index][blocknum][i][coef_index]++;
  }
//   fc->vec->qp[qp]++;

//   if (mb_type&0x0002)
//     fprintf(fc->file, "%i: is intra 16x16! ", mb_xy);
//   else
//     fprintf(fc->file, "%i:                 ", mb_xy);
//   
//   fprintf(fc->file, "n = %i, total_coeff = %i, totalZeros = %i ", n, total_coeff, totalZeros);
//   for (i = 0; i < total_coeff; i++) {
//     if (level[i] < fc->vec->min) fc->vec->min = level[i];
//     if (level[i] > fc->vec->max) fc->vec->max = level[i];
//     if (level[i] >= -50 && level[i] < 50) {
//       fc->vec->v[level[i]+50] += 1;
//     }
//     fprintf(fc->file, "(%i, %i)", level[i], run_before[i]);
//   }
//   fprintf(fc->file, "coef = {");
//   for (i = 0; i < 16; i++) {
//     fprintf(fc->file, "%i, ", tape[i]);
//   }
//   fprintf(fc->file, "} ");
//   fprintf(fc->file, "\n");
// //   fc->vec->N++;
//   
//   if (mb_xy < 396) {
//     fc->vec->mb_t[mb_xy] = mb_type;
//   } else {
//     printf("Etwas stimmt hiuer nicht! %d \n", mb_xy);
//   }
  
//   if (0 <= n && n <= 15)   blocknum = n;   // Luma default
//   if (n == 49)             blocknum = 16;  // Cr DC
//   if (16 <= n && n <= 19)  blocknum = n;   // Cr AC
//   if (n == 50)             blocknum = 19;  // Cr DC
//   if (16 <= n && n <= 19)  blocknum = n;   // Cr AC

//   fc->vec->v[(qp-17)*1600 + blocknum*100 + 50]
}

void storeCounts(H264FeatureContext *fc) {
  int i, j, k, l;
  
  if (!fc->refreshed) return;  // avoid zero-vector at the beginning
  
  fprintf(fc->file, "I-Slices: %i, P-Slices: %i, B-Slices: %i    ", fc->i_slices, fc->p_slices, fc->b_slices);
  
  fprintf(fc->file, "min = %i, max = %i, vector_num = %i, slice_type = %i, histograms:\n", fc->vec->min, fc->vec->max, fc->vec->vector_num, fc->vec->slice_type);
    // store qp histogram
  for (i = 0; i < 50; i++) {
    if (i == 20) 
      fprintf(fc->file, " [%u] ", fc->vec->qp[i]);
    else 
      fprintf(fc->file, " %u ", fc->vec->qp[i]);
  }
  fprintf(fc->file, "\n");
//   for (i = 0; i < 100; i++) {
//     fprintf(fc->file, " %i ", fc->vec->v[i]);
//   }
//   fprintf(fc->file, "\n");
  for (i = 0; i < QP_RANGE; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < num_coefs[j]; k++) {
	// don't save zero hestograms
	if (fc->vec->N[i][j][k] == 0) continue;
	
	fprintf(fc->file, "qp = %i, blockn = %i, N = %i --- ", QP_OFFSET + i, j, fc->vec->N[i][j][k]);
	for (l = 0; l < fc->histogram_ranges[j][k]; l++) {
	  fprintf(fc->file, "%i ", fc->vec->histogram[i][j][k][l]);
	}
	fprintf(fc->file, "[0] ");
	for (l = fc->histogram_ranges[j][k]; l < 2*fc->histogram_ranges[j][k]; l++) {
	  fprintf(fc->file, "%i ", fc->vec->histogram[i][j][k][l]);
	}
	fprintf(fc->file, "\n");
      }
    }
  }
  
  fc->refreshed = 0;
  
//   // store mb types
//   for (i = 0; i < 396; i++) {
//     fprintf(fc->file, " %u ", fc->vec->mb_t[i]);
//   }//*/
//   fprintf(fc->file, "\n");
  
//   fflush(fc->file);
}

void refreshFeatures(H264FeatureContext* feature_context) {
  int i, j, k, l;
  
  feature_context->i_slices = 0;
  feature_context->p_slices = 0;
  feature_context->b_slices = 0;
  feature_context->vec->vector_num++;
  feature_context->vec->min = 0;
  feature_context->vec->max = 0;
//   feature_context->vec->vector_num = 0;
//   feature_context->vec->N = 0;
//   for (int i = 0; i < 100; i++) {
//     feature_context->vec->v[i] = 0;
//   }
//   for (i = 0; i < feature_context->vec->num_histograms; i++) {
//     feature_context->vec->N[i] = 0;
//   }
  for (i = 0; i < QP_RANGE; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < num_coefs[j]; k++) {
	feature_context->vec->N[i][j][k] = 0;
	for (l = 0; l < 2*feature_context->histogram_ranges[j][k]; l++) {
	  feature_context->vec->histogram[i][j][k][l] = 0;
	}
      }
    }
  }
  for (i = 0; i < 396; i++) {
    feature_context->vec->mb_t[i] = 0;
  }
  for (i = 0; i < 50; i++) {
    feature_context->vec->qp[i] = 0;
  }
  
  feature_context->refreshed = 1;
}


