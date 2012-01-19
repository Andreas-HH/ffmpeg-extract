#include "h264_extract.h"

/*
 (Luma DC:   48)
 Luma AC:   0 - 15
 Chroma DC: 49, 50
 Chroma AC: 16 - 19, 32 - 35
 */


/*void perform_hiding_plusminus(void *vfc) {
  int i;
  H264FeatureContext *fc = (H264FeatureContext*) vfc;
  int *coefs;// = fc->tape;
  double r;
  
  printf("thread: locking main mutex \n");
  pthread_mutex_lock(fc->main_mutex);    // tell main thread you are ready to do work
  
  while (1) { // need some variable in fc that tells thread to terminate
    printf("thread: locking thread mutex \n");
    pthread_mutex_lock(fc->thread_mutex);  // waiting for main to prepare coefs
    printf("thread: unlocking main mutex \n");
    pthread_mutex_unlock(fc->main_mutex);    // tell main thread you are ready to do work
//     printf("thread: locking main mutex \n");
//     pthread_mutex_lock(fc->main_mutex);  // stop main from doing stupid things
    printf("Hallo \n");
    if (fc->blocknum == -1) {
      printf("thread: unlocking thread mutex (blocknum) \n");
      pthread_mutex_unlock(fc->thread_mutex);  // tell main that I'm finished
      printf("thread: locking main mutex \n");
      pthread_mutex_unlock(fc->main_mutex);  // tell main that I'm finished
      continue;  // return
    }
    
    coefs = fc->tape;
    memcpy(fc->tape, fc->proper_coefs, 16*sizeof(int));
    for (i = 0; i < num_coefs[fc->blocknum]; i++) {
      if (coefs[i]<2 && coefs[i]>-2) continue;
      r = (((double) rand()) / ((double) RAND_MAX));
  //     printf("%f \n ", r);
      if (r < fc->p_hide) {
	switch (fc->slice_type) {
	  case TYPE_P_SLICE:
	    fc->hidden_bits_p++;
	    break;
	  case TYPE_B_SLICE:
	    fc->hidden_bits_b++;
	    break;
	}
	r = (((double) rand()) / ((double) RAND_MAX));
	if (r < PROB_KEEP) {   // decide to keep or change value
	  r = (((double) rand()) / ((double) RAND_MAX));
	  if (r < PROB_INCREASE) {  // if changing, increase or decrease?
	    coefs[i] += 1;
	  } else {
	    coefs[i] -= 1;
	  }
	}
      }
    }
    printf("Habe for überlebt! \n");
    addCounts(fc, fc->current_qp, fc->blocknum); 
    printf("thread: unlocking thread mutex \n");
    pthread_mutex_unlock(fc->thread_mutex);  // tell main that I'm finished
    printf("thread: locking main mutex \n");
    pthread_mutex_unlock(fc->main_mutex);  // tell main that I'm finished
  }
//   pthread_exit(NULL);
}*/

void myprint(char *text) {
   printf("%s", text);
}

// void simulate_hiding_plusminus(H264FeatureContext *fc) {
//   int rc;
// 
//   rc = pthread_create(fc->thread, fc->thread_attr, (void *) &perform_hiding_plusminus, (void *) fc);
//   if (rc) printf("ERROR: pthread_create failed! %i \n", rc);
// }

void simulate_hiding_plusminus(H264FeatureContext *fc) {
  int i;
  int *coefs = fc->tape;
  double r;
  int sl = fc->slice_type;
  
  if (sl == TYPE_I_SLICE) return;
//   if (sl == TYPE_B_SLICE) return;   // maybe wis hto change that later
  if (!(fc->accept_blocks & (1 << fc->blocknum))) {
//     if (fc->accept_blocks == ACCEPT_LC && fc->blocknum!=-1)
//     printf("rejected block of type %i, accept=%i \n", fc->blocknum, fc->accept_blocks);
    return;
  }
//   if (fc->accept_blocks == 7)
//     printf("accepted block %i \n", fc->blocknum);
  
//   memcpy(fc->tape, fc->proper_coefs, 16*sizeof(int));
  for (i = 0; i < num_coefs[fc->blocknum]; i++) {
    if (coefs[i]<2 && coefs[i]>-2) continue;
    r = (((double) rand()) / ((double) RAND_MAX));
//     printf("%f \n ", r);
    if (r < fc->p_hide) {
      switch (fc->slice_type) {
	case TYPE_P_SLICE:
	  fc->hidden_bits_p++;
	  break;
	case TYPE_B_SLICE:
	  fc->hidden_bits_b++;
	  break;
      }
      r = (((double) rand()) / ((double) RAND_MAX));
      if (r < PROB_KEEP) {   // decide to keep or change value
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

// void wait_for_simulation(H264FeatureContext* fc) {
//   int rc;
//   void *thread_status;
//   
//   rc = pthread_join(*(fc->thread), &thread_status);
//   if (rc) printf("ERROR: pthread_join failed: %i \n", rc);
// }


/*void setup_ranges(int **ranges, int luma, int chroma_dc, int chroma_ac) {
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
}*/

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
  
//   if (blocknum == 2) // Chroma AC
//     pos = total_coeff + totalZeros;
//   else
  pos = total_coeff + totalZeros - 1;
  zerosToSee = totalZeros;
  result[pos] = level[0];
  pos--;
  for (i = 1; i < total_coeff; i++) {
    if (zerosToSee > 0) {
      zerosToSee -= run_before[i];
      pos -= run_before[i];
    }
//     if (level[i] < vec->min) vec->min = level[i];
//     if (level[i] > vec->max) vec->max = level[i];
    result[pos] = level[i];
    pos--;
  }
}

void addCounts(H264FeatureContext *fc, int qp, int blocknum, int len) {
  // need len to have correct pair counts
//   fflush(stderr);
//   printf("adding counts \n");
  int i, l, r;
  int coef_index;
  int *tape = fc->tape;
//   int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  int sl = fc->slice_type;
  
  if (blocknum == -1 || qp_index < 0 || qp_index >= QP_RANGE)
    return;
//     if (qp_index < 0 || qp_index >= QP_RANGE) {
// //     fprintf(fc->file, "qp out of range o_O %i \n", qp);
//     return;
//   }
  
//   if (sl != TYPE_P_SLICE) return;  // add only P-Slices 
  if (sl == TYPE_I_SLICE) return;
  switch (sl) {
    case TYPE_P_SLICE:
      fc->num_4x4_blocks_p++;
      break;
    case TYPE_B_SLICE:
      fc->num_4x4_blocks_b++;
      break;  
  }
  
//   fc->vec->qp[qp]++;

//   if (blocknum == -1) {
// //     fprintf(fc->file, "blocknum is -1 o_O %i \n", n);
//     return;
//   }

  
//   printf("extracted qp histogram %i, %i \n", sl, blocknum);
  
//   constructProperCoefArray(tape, level, run_before, total_coeff, totalZeros, blocknum, fc->vec);
  if (fc->extract_rate || fc->files_hist[sl] != NULL) {  // only add counts if results get stored
    for (i = 0; i < len; i++) { // num_coefs[blocknum]
//       if (blocknum == 2 && i == 0) continue;            // We don't want to count the fake 0 of chroma ac
      coef_index = tape[i];
      if (coef_index == 0) continue;
      else if (coef_index < 0)  {
	coef_index = coef_index + ranges[blocknum][i]; // fc->histogram_
	if (coef_index < 0) continue;
      }
      else if (coef_index > 0)  {
	coef_index = coef_index + ranges[blocknum][i]-1; // fc->histogram_
	if (coef_index > 2*ranges[blocknum][i]-1) continue; // fc->histogram_
      }
//       if (blocknum == 2)
// 	fc->vec->histograms[sl][qp_index][blocknum][i-1][coef_index]++; // Chroma AC always has coef_0 = 0
//       else 
	fc->vec->histograms[sl][qp_index][blocknum][i][coef_index]++;
    }
  }
//   printf("extracted hist ");
  //pairs
  if (fc->extract_rate || fc->files_pair[sl] != NULL) {  // only add counts if results get stored
    for (i = 0; i < len-1; i++) { // num_coefs[blocknum]
//       if (blocknum == 2 && i == 0) continue;            // We don't want to count the fake 0 of chroma ac
      l = tape[i] + ranges[blocknum][0];                // we are allowed to exceed the local range here
      r = tape[i+1] + ranges[blocknum][1];              // space is limited by ranges of first two coefs
      if (l < 0 || l > 2*ranges[blocknum][0]) continue;
      if (r < 0 || r > 2*ranges[blocknum][1]) continue;
      fc->vec->pairs[sl][qp_index][blocknum][l][r]++;
    }
  }
  
//   printf("and pairs \n");
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
 	if (ranges[j][k] == 0 || fc->vec->histograms[0][i][j][k][ranges[j][k]] == 0) continue;
// 	if (fc->vec->N[i][j][k] == 0) continue;
	
	fprintf(fc->logfile, "qp = %i, blockn = %i, --- ", QP_OFFSET + i, j);//, fc->vec->N[i][j][k]);
	for (l = 0; l < ranges[j][k]; l++) {
	  fprintf(fc->logfile, "%i ", fc->vec->histograms[0][i][j][k][l]);
	}
	fprintf(fc->logfile, "[0] ");
	for (l = ranges[j][k]; l < 2*ranges[j][k]; l++) {
	  fprintf(fc->logfile, "%i ", fc->vec->histograms[0][i][j][k][l]);
	}
	fprintf(fc->logfile, "\n");
      }
    }
  }
  
  fc->refreshed = 0;
}

void storeFeatureVectors(H264FeatureContext* fc) {
  int bin;
  int sl, i, j, k, l;
  int count_h, count_p;
  int N_h, N_p;
  long pos;
  double sum_h, sum_p;
  double scale_h, scale_p;
  FILE *current_hist_file;
  FILE *current_pair_file;
  
//   if (fc->accept_blocks == 7)
//     printf("Storing blocks in ACCEÜT_LC \n");
  
//   printf("Storing feature vectors.. \n");
  
  if (!fc->refreshed) {
//     printf("not refreshed! \n");
    return;
  }
  
//   printf("I am fresh! \n");
  
  for (sl = 0; sl < 2; sl++) {  // DROPPING B-FRAMES !!!
    count_h = 0;
    count_p = 0;
    N_h = 0;
    N_p = 0;
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
//         if (fc->extract_rate || fc->files_hist[sl] != NULL) {
	  for (k = 0; k < num_coefs[j]; k++) {
	    for (l = 0; l < 2*ranges[j][k]; l++) {
  //             printf("(%i, %i, %i, %i) -- (%i, %i)\n", i, j, k, l, count_h, count_p);
	      fc->vec->vector_histograms[count_h] = (double) fc->vec->histograms[sl][i][j][k][l];
	      count_h++;
	      N_h += fc->vec->histograms[sl][i][j][k][l];
  //             fprintf(fc->files_hist[sl], "%i ", fc->vec->histograms[sl][i][j][k][l]);
	    }
	  }
// 	}
        // pairs
//         if (fc->extract_rate || fc->files_pair[sl] != NULL) {
	  for (k = 0; k < 2*ranges[j][0]+1; k++) {
	    for (l = 0; l < 2*ranges[j][1]+1; l++) {
	      if (k == ranges[j][0] && l == ranges[j][1]) continue; // 00: &&; 0x x0: ||
  //             printf("(%i, %i)\n", count_h, count_p);
	      fc->vec->vector_pairs[count_p] = (double) fc->vec->pairs[sl][i][j][k][l];
	      count_p++;
	      N_p += fc->vec->pairs[sl][i][j][k][l];
  //             fprintf(fc->files_pair[sl], "%i ", fc->vec->pairs[sl][i][j][k][l]);
	    }
	  }
// 	}
      }
    }
//     printf("found: (%i, %i) calculated: (%i, %i)\n", count_h, count_p, fc->vec->vector_histograms_dim, fc->vec->vector_pairs_dim);
    
    
//     printf("Just did the scaling \n");
  
    // check if vector sums orrectly, for debugging only
    /*sum_h = 0.;
    for (i = 0; i < fc->vec->vector_histograms_dim; i++) {
      sum_h += fc->vec->vector_histograms[i];
    }
    sum_p = 0.;
    for (i = 0; i < fc->vec->vector_pairs_dim; i++) {
      sum_p += fc->vec->vector_histograms[i];
    }*/
//     printf("dim_h = %i,dim_p = %i, sum_h = %f, sum_p = %f,  N_h = %i, N_p = %i \n", fc->vec->vector_histograms_dim, fc->vec->vector_pairs_dim, sum_h, sum_p, N_h, N_p);
//     if (fc->accept_blocks == ACCEPT_LC)
//       printf("LC: hidden bits: %i, numblocks: %i \n", fc->hidden_bits_p, fc->num_4x4_blocks_p);
    if (sl == 0) if (fc->num_4x4_blocks_p == 0) {
//       printf("N_h = 0, accept = %i, hidden bits_p: %i, numblocks: %i \n", fc->accept_blocks, fc->hidden_bits_p, fc->num_4x4_blocks_p);
      continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!
    }
    if (sl == 1) if (fc->num_4x4_blocks_b == 0) continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!
   
    scale_h = ((double) fc->vec->vector_histograms_dim)/((double) N_h);
    for (i = 0; i < fc->vec->vector_histograms_dim; i++)
      fc->vec->vector_histograms[i] *= scale_h;
    scale_p = ((double) fc->vec->vector_pairs_dim)/((double) N_p);
    for (i = 0; i < fc->vec->vector_pairs_dim; i++)
      fc->vec->vector_pairs[i] *= scale_p;
    
//       	if (fc->accept_blocks == ACCEPT_LC)
// 	  printf("LC: wtf_2 \n");
    if (fc->extract_rate) {
      if (sl == 0) {
// 	printf("%f: p-rate=%f ", fc->p_hide, ((double)fc->hidden_bits_p)/((double)fc->num_4x4_blocks_p));
	bin = get_rate_index(((double)fc->hidden_bits_p)/((double)fc->num_4x4_blocks_p));
      }
      if (sl == 1) 
	bin = get_rate_index(((double)fc->hidden_bits_b)/((double)fc->num_4x4_blocks_b));
//             printf("%f: extract rate! bin=%i ", fc->p_hide, bin);
//       printf("found bin: %i \n", bin);
//       	if (fc->accept_blocks == ACCEPT_LC)
// 	  printf("LC: looking for files, bin=%i \n", bin);
      if (bin >= 0) {
	current_hist_file = fc->rate_bins_hist[fc->accept_blocks][sl][bin];
        current_pair_file = fc->rate_bins_pair[fc->accept_blocks][sl][bin];
      }
      else {
	current_hist_file = NULL;
	current_pair_file = NULL;
      }
    } else {
//       printf("%f, extract prob ", fc->p_hide);
//       printf("not using bins! \n");
      current_hist_file = fc->files_hist[sl];
      current_pair_file = fc->files_pair[sl];
    }
    
//     printf("%f: Have current_files, start writing \n", fc->p_hide);
//     if (current_hist_file == NULL) printf("%f: hist is NULL! \n", fc->p_hide);
//     if (current_pair_file == NULL) printf("%f: pair is NULL! \n", fc->p_hide);
//     current_file = NULL;
//     printf("b-rate: %f, p-rate: %f \n", ((double)fc->hidden_bits_b)/((double)fc->num_4x4_blocks_b), ((double)fc->hidden_bits_p)/((double)fc->num_4x4_blocks_p));
//     if (fc->files_hist[sl] != NULL) {
    if (current_hist_file != NULL) {
//       if (fc->accept_blocks == ACCEPT_LC)
// 	printf("Storing some LC vector \n");
//       printf(" fc->files_hist[%i] != NULL \n", sl);
//       pos = ftell(current_hist_file);
// //       printf(" Having pos \n");
//       if (pos == 0L)
// 	fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, current_hist_file);
      if (N_h != 0) {
  
// 	printf(" Writing hist of slice %i \n", sl);
//         if (sl == 0) {
// 	bin = get_rate_index(((double)fc->hidden_bits_p)/((double)fc->num_4x4_blocks_p));
// 	if (bin >= 0) 
// 	  fwrite(fc->rate_bins_hist[fc->accept_blocks][sl][bin], sizeof(double), fc->vec->vector_histograms_dim, fc->files_hist[sl]);
// 	}
//         printf("writung %i doubles! \n", fc->vec->vector_histograms_dim);
// 	if (fc->p_hide < 0) {
// 	  for (i = 0; i < fc->vec->vector_histograms_dim; i++)
// 	    printf("%f ", fc->vec->vector_histograms[i]);
// 	  printf("\n");
// 	}
	fwrite(fc->vec->vector_histograms, sizeof(double), fc->vec->vector_histograms_dim, current_hist_file);
	
      }
    }
//     printf("Stored histograms \n");
    if (current_pair_file != NULL) {
      pos = ftell(current_pair_file);
      if (pos == 0L)
	fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, current_pair_file);
      if (N_p != 0) {
	fwrite(fc->vec->vector_pairs, sizeof(double), fc->vec->vector_pairs_dim, current_pair_file);
      }
    }
//     if (fc->files_pair[sl] != NULL) {
// //       printf(" fc->files_pair[%i] != NULL \n", sl);
//       pos = ftell(fc->files_pair[sl]);
//       if (pos == 0L)
// 	fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, fc->files_pair[sl]);
//       if (N_p != 0) {
// // 	printf(" Writing pair of slice %i \n", sl);
// 	fwrite(fc->vec->vector_pairs, sizeof(double), fc->vec->vector_pairs_dim, fc->files_pair[sl]);
//       }
//     }
  }
//   printf("fertig!\n");
  
  fc->refreshed = 0;
}


void refreshFeatures(H264FeatureContext* fc) {
  int i, j, k, l, sl;
  
//   printf("Refreshing features \n");
  
  fc->i_slices = 0;
  fc->p_slices = 0;
  fc->b_slices = 0;
  fc->vec->vector_num++;
  fc->vec->min = 0;
  fc->vec->max = 0;
  fc->num_4x4_blocks_b = 0;
  fc->num_4x4_blocks_p = 0;
  fc->hidden_bits_b = 0;
  fc->hidden_bits_p = 0;

  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
//           feature_context->vec->N[i][j][k] = 0;
          for (l = 0; l < 2*ranges[j][k]; l++) {
            fc->vec->histograms[sl][i][j][k][l] = 0;
          }
        }
        // pairs
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          for (l = 0; l < 2*ranges[j][1]+1; l++) {
            fc->vec->pairs[sl][i][j][k][l] = 0;
          }
        }
      }
    }
  }
//   for (i = 0; i < fc->vec->vector_histograms_dim; i++)
//     fc->vec->vector_histograms[i] = 0.;
//   for (i = 0; i < fc->vec->vector_pairs_dim; i++)
//     fc->vec->vector_pairs[i] = 0.;
//   for (i = 0; i < 396; i++) {
//     feature_context->vec->mb_t[i] = 0;
//   }
  for (i = 0; i < 50; i++) {
    fc->vec->qp[i] = 0;
  }
  
//   printf("refresh: fertif \n");
//   pthread_mutex_lock(fc->thread_mutex);
//   pthread_create(fc->thread, fc->thread_attr, (void *) &perform_hiding_plusminus, (void *) fc);
  
  fc->refreshed = 1;
}

FILE**** init_rate_bins(char* type, char *method_name, int dim) {
  int i, j;
  FILE ****result;
  char *blockstring;// = blockstrings[accept];
  char path[512];
  
  result = (FILE****) malloc(8*sizeof(FILE***)); // accept: [0..7]
  for (i = 0; i < 8; i++) {
    if (i != ACCEPT_L && i != ACCEPT_C && i != ACCEPT_LC) continue;
    blockstring = blockstrings[i];
//     printf("blockstring = %s \n", blockstring);
    result[i] = (FILE***) malloc(2*sizeof(FILE**)); // slice_type
    
    result[i][0] = (FILE**) malloc(NUM_BINS*sizeof(FILE*));
    for (j = 0; j < NUM_BINS; j++) {
//       printf("j=%i, %f \n", j, (((double) (j))/((double) NUM_BINS))*MAX_RATE);
      sprintf(path, "%s/%s/rate/p_%s/p_%s_rate_%i.fv", method_name, blockstring, type, type, (int) ((((double) (10000*j))/((double) NUM_BINS))*MAX_RATE));
//       printf("opening rate file: %s \n", path);
      result[i][0][j] = fopen(path, "a");
      if (result[i][0][j] != NULL && ftell(result[i][0][j]) == 0L) 
        fwrite(&dim, sizeof(int), 1, result[i][0][j]);
    }
    
    result[i][1] = (FILE**) malloc(NUM_BINS*sizeof(FILE*));
    for (j = 0; j < NUM_BINS; j++) {
      sprintf(path, "%s/%s/rate/b_%s/b_%s_rate_%i.fv", method_name, blockstring, type, type, (int) ((((double) (10000*j))/((double) NUM_BINS))*MAX_RATE));
//       result[i][1][j] = fopen(path, "a");
      result[i][1][j] = NULL;
      if (result[i][1][j] != NULL && ftell(result[i][1][j]) == 0L) 
        fwrite(&dim, sizeof(int), 1, result[i][1][j]);
    }
  }
  
  return result;
}

void close_rate_bins(FILE**** bins) {
  int i, j;
  
//   printf("Closing bins \n");
  
  for (i = 0; i < 8; i++) {
    if (i != ACCEPT_L && i != ACCEPT_C && i != ACCEPT_LC) continue;
    for (j = 0; j < NUM_BINS; j++) {
      if (bins[i][0][j] != NULL) {
	fflush(bins[i][0][j]);
	fclose(bins[i][0][j]);
      }
    }
//     printf("close: closed p files \n");
    free(bins[i][0]);
    
    for (j = 0; j < NUM_BINS; j++) {
      if (bins[i][1][j] != NULL) {
	fflush(bins[i][1][j]);
	fclose(bins[i][1][j]);
      }
    }
//     printf("close: closed b  files \n");
    free(bins[i][1]);
    free(bins[i]);
  }
  free(bins);
}



H264FeatureContext* init_features(char* method_name, int accept_blocks, double p_hide, FILE**** bins_h, FILE**** bins_p, int extract_rate) {
  int i, j, k, sl;
  H264FeatureContext *fc;
  H264FeatureVector *fv;
//   int **ranges;
//   int ***N;
  int *****v;
  int *****w;
  int hist_dim = 0;
  int pair_dim = 0;
//   int *dim = &hist_dim;
//   int *mb_t;
  int *qp;
  int *tape;
  char b_h_path[512];
  char b_p_path[512];
  char p_h_path[512];
  char p_p_path[512];
  char *blockstring = blockstrings[accept_blocks];//"L";
//   int num_histograms = 16+2*(4+15);
  
//   printf("Initialising some featureset... ");
  
//   ranges    = malloc(3*sizeof(int*));
//   ranges[0] = malloc(16*sizeof(int));
//   ranges[1] = malloc(4*sizeof(int));
//   ranges[2] = malloc(15*sizeof(int));
//   setup_ranges(ranges, LUMA_RANGE, CHROMA_DC_RANGE, CHROMA_AC_RANGE);
  
//   N = av_malloc(QP_RANGE*sizeof(int**));
  v = malloc(2*sizeof(int****));
  w = malloc(2*sizeof(int****));
  for (sl = 0; sl < 2; sl++) {
    v[sl] = malloc(QP_RANGE*sizeof(int***));
    w[sl] = malloc(QP_RANGE*sizeof(int***));
    for (i = 0; i < QP_RANGE; i++) {
  //     N[i] = av_malloc(3*sizeof(int*));
      // setup histograms
      v[sl][i] = malloc(3*sizeof(int**));   // 3 types of blocks (Luma + Chroma DC/AC)
      w[sl][i] = malloc(3*sizeof(int**));
      for (j = 0; j < 3; j++) {
  //       N[i][j] = av_malloc(num_coefs[j]*sizeof(int));
        // setup histograms
        v[sl][i][j] = malloc(num_coefs[j]*sizeof(int*));
        for (k = 0; k < num_coefs[j]; k++) {
          v[sl][i][j][k] = malloc(2*ranges[j][k]*sizeof(int));
        }
        // setup pairs
        w[sl][i][j] = malloc((2*ranges[j][0]+1)*sizeof(int*)); // +1 for the zero
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          w[sl][i][j][k] = malloc((2*ranges[j][1]+1)*sizeof(int));
        }
      }
    }
  }
  
  for (i = 0; i < 3; i++) { // scan through blocks
    for (j = 0; j < num_coefs[i]; j++) {
      hist_dim += 2*ranges[i][j];
    }
    pair_dim += (2*ranges[i][0]+1) * (2*ranges[i][1]+1) - 1; // (0,0) not included
  }
//   v     =   av_malloc(USED_PIXELS*FEATURE_DIMENSION*QP_RANGE*sizeof(feature_elem));
//   mb_t  =   av_malloc(396*sizeof(int));
  qp    =   malloc(50*sizeof(int));
  fc    =   malloc(sizeof(H264FeatureContext));
  fv    =   malloc(sizeof(H264FeatureVector));
  tape  =   malloc(16*sizeof(int));  
//   fv->N = N;
  
//   FILE *tmp = fopen("p.fv", "r");
//   fread(&hist_dim, sizeof(int), 1, tmp);
//   fclose(tmp);
//    *dim=sizeof(int);
  
  fv->vector_histograms_dim = QP_RANGE*hist_dim;
  fv->vector_pairs_dim      = QP_RANGE*pair_dim;
  fv->vector_histograms     = malloc(fv->vector_histograms_dim*sizeof(double));
  fv->vector_pairs          = malloc(fv->vector_pairs_dim*sizeof(double));
  fv->histograms     = v;
  fv->pairs          = w;
//   fv->mb_t         = mb_t;
  fv->qp             = qp;
  fc->vec            = fv;
//   fc->thread_mutex   = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t));  // unlock from main notifies all threads that try to lock this mutex
//   fc->main_mutex     = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t));  // used to notfy main thread that all work is finished in detached thread
//   fc->thread_attr    = (pthread_attr_t*) malloc(sizeof(pthread_attr_t));
//   fc->thread         = (pthread_t*) malloc(sizeof(pthread_t));
//   fc->logfile        = fopen("features.log", "w"); // fclose(file)
  fc->files_hist     = malloc(2*sizeof(FILE*));
  fc->files_pair     = malloc(2*sizeof(FILE*));
  fc->extract_rate   = extract_rate;
  if (p_hide < 0) {
    fc->files_hist[0]  = fopen("clean/p_histograms.fv", "a");
//     fc->files_hist[1]  = fopen("clean/b_histograms.fv", "a");
    fc->files_hist[1] = NULL;
    fc->files_pair[0]  = fopen("clean/p_pairs.fv", "a");
//     fc->files_pair[1]  = fopen("clean/b_pairs.fv", "a");
    fc->files_pair[1] = NULL;
  } else {
//     printf("opening files ");
    if (extract_rate) {
//       printf("using rate \n");
      fc->rate_bins_hist = bins_h;
      fc->rate_bins_pair = bins_p;
      
      fc->files_hist[0] = NULL;
      fc->files_hist[1] = NULL;
      fc->files_pair[0] = NULL;
      fc->files_pair[1] = NULL;
    } else {
//       printf("opening prob files \n");
      sprintf(p_h_path, "%s/%s/prob/p_hist/p_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(p_p_path, "%s/%s/prob/p_pair/p_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(b_h_path, "%s/%s/prob/b_hist/b_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(b_p_path, "%s/%s/prob/b_pair/b_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
  //     printf("%s", p_h_path);
      fc->files_hist[0]  = fopen(p_h_path, "a");
      fc->files_hist[1]  = fopen(b_h_path, "a");
//       fc->files_hist[1] = NULL;
      fc->files_pair[0]  = fopen(p_p_path, "a");
      fc->files_pair[1]  = fopen(b_p_path, "a");
  //     fc->files_pair[0] = NULL;
//       fc->files_pair[1] = NULL;
    }
  }
  
  if (fc->files_hist[0] != NULL && ftell(fc->files_hist[0]) == 0L) 
    fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_hist[0]);
  if (fc->files_hist[1] != NULL && ftell(fc->files_hist[1]) == 0L) 
    fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_hist[1]);
  if (fc->files_pair[0] != NULL && ftell(fc->files_pair[0]) == 0L) 
    fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_pair[0]);
  if (fc->files_pair[1] != NULL && ftell(fc->files_pair[1]) == 0L) 
    fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_pair[1]);
  
//   fc->p_hist = fopen("p_histograms.fv", "a");
//   fc->p_pair = fopen("p_pairs.fv", "a");
//   fc->b_hist = fopen("b_histograms.fv", "a");
//   fc->b_pair = fopen("b_pairs.fv", "a");
  fc->tape = tape;
//   fc->histogram_ranges = ranges;
  fc->accept_blocks = accept_blocks;
  fc->p_hide = p_hide;
//   fc->rate_bins_hist = (FILE****) malloc(8*sizeof(FILE***));
//   fc->rate_bins_pair = (FILE****) malloc(8*sizeof(FILE***));
  
  fc->i_slices    = 0;
  fc->p_slices    = 0;
  fc->b_slices    = 0;
  fv->vector_num  = 0;
  fc->refreshed   = 0;
  
//   printf("Done. \n");
//   if (!(p_hide < 0)) {
//     pthread_mutex_init(fc->thread_mutex, NULL);
//     pthread_mutex_init(fc->main_mutex, NULL);
//     pthread_attr_init(fc->thread_attr);
//     pthread_attr_setdetachstate(fc->thread_attr, PTHREAD_CREATE_DETACHED);
//     printf("init: locking thread mutex \n");
//     pthread_mutex_lock(fc->thread_mutex);
//     pthread_create(fc->thread, fc->thread_attr, (void *) &perform_hiding_plusminus, (void *) fc);
//     fc->locked_by_init = 1;
//   }
//   printf("init done! \n");
  
  return fc;
//   h->feature_context = fc;
}

void close_features(H264FeatureContext* fc) {
  int i, j, k, sl;
//   H264FeatureContext *fc = h->feature_context;
//   printf("Closing features... \n");
//   storeCounts(fc);
  storeFeatureVectors(fc);
  
//   printf("close: storing done \n");
//   fflush(fc->logfile);
//   fclose(fc->logfile);
//   printf("close: done logfile \n");
  if (fc->files_hist[0] != NULL) {
    fflush(fc->files_hist[0]);
    fclose(fc->files_hist[0]);
  }
//   printf("close: done p_hist \n");
  if (fc->files_hist[1] != NULL) {
    fflush(fc->files_hist[1]);
    fclose(fc->files_hist[1]);
  }
//   printf("close: done b_hist \n");
  if (fc->files_pair[0] != NULL) {
    fflush(fc->files_pair[0]);
    fclose(fc->files_pair[0]);
  }
//   printf("close: done p_pair \n");
  if (fc->files_pair[1] != NULL) {
    fflush(fc->files_pair[1]);
    fclose(fc->files_pair[1]);
  }
//   printf("close: done b_pair \n");

//   printf("close: closed all files \n");
  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < num_coefs[j]; k++) {
          free(fc->vec->histograms[sl][i][j][k]);
        }
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          free(fc->vec->pairs[sl][i][j][k]);
        }
  //       av_free(fc->vec->N[i][j]);
        free(fc->vec->histograms[sl][i][j]);
        free(fc->vec->pairs[sl][i][j]);
      }
  //     av_free(fc->vec->N[i]);
      free(fc->vec->histograms[sl][i]);
      free(fc->vec->pairs[sl][i]);
    }
  //   av_free(fc->vec->N);
    free(fc->vec->histograms[sl]);
    free(fc->vec->pairs[sl]);
  }
//   printf("close: fc->vec nearly down \n");
  free(fc->vec->histograms);
  free(fc->vec->pairs);
  free(fc->vec->vector_histograms);
  free(fc->vec->vector_pairs);
  
//   free(fc->histogram_ranges[0]);
//   free(fc->histogram_ranges[1]);
//   free(fc->histogram_ranges[2]);
//   free(fc->histogram_ranges);
//   printf("close: fast fertig \n");
//   free(fc->thread_attr);
//   free(fc->thread);
  free(fc->tape);
//   av_free(fc->vec->mb_t);
  free(fc->vec->qp);
  free(fc->vec);
  free(fc);
}

