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
      if (coefs[i] == 2) {
	coefs[i]++;
      } else if (coefs[i] == -2) {
	coefs[i]--;
      } else{
// 	r = (((double) rand()) / ((double) RAND_MAX));
// 	if (r < PROB_KEEP) {   // decide to keep or change value
	r = (((double) rand()) / ((double) RAND_MAX));
	if (r < PROB_INCREASE) {  // if changing, increase or decrease?
	  coefs[i] += 1;
	} else {
	  coefs[i] -= 1;
	}
// 	}
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
  
//   printf("constructing coef array on blocknum: %i \n", blocknum);
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

void addCounts(H264FeatureContext *fc, int qp, int n, int len) {
  // need len to have correct pair counts
//   fflush(stderr);
//   printf("adding counts \n");
  int i, l, r;
  int coef_index;
  int *tape = fc->tape;
  int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  int sl = fc->slice_type;
  
  if (blocknum == -1 || qp_index < 0 || qp_index >= QP_RANGE)
    return;
//     if (qp_index < 0 || qp_index >= QP_RANGE) {
// //     fprintf(fc->file, "qp out of range o_O %i \n", qp);
//     return;
//   }
  if (blocknum == 1) {
    printf("n = %i, x = %i, y = %i \n", n, fc->x, fc->y);
  }
  if (n == 49) {
    memcpy(fc->lastU, tape, num_coefs[1]*sizeof(int));
    fc->seenU = 1; // it could happen that we skip one frame and align perfectly, 
  }                // VERY unlikely though. Comparing x,y with ux,uy should be enough
  if (n == 50 && fc->seenU && fc->ux == fc->x && fc->uy == fc->y) {
    printf("I could count U vs V now! \n");
    for (i = 0; i < num_coefs[1]; i++) printf("%i, ", fc->lastU[i]);
    printf("\n");
    for (i = 0; i < num_coefs[1]; i++) printf("%i, ", tape[i]);
    printf("\n");
    fc->seenU = 0;
  }
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
  
//   constructProperCoefArray(tape, level, run_before, total_coeff, totalZeros, blocknum, fc->vec);
  if (fc->extract_rate || fc->files_hist[sl] != NULL) {  // only add counts if results get stored
    for (i = 0; i < len; i++) { // num_coefs[blocknum]
//       if (blocknum == 2 && i == 0) continue;            // We don't want to count the fake 0 of chroma ac
      coef_index = tape[i];
//       if (coef_index == 0) continue;
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
  }
  //pairs
  if (fc->extract_rate || fc->files_pair[sl] != NULL) {  // only add counts if results get stored
    for (i = 0; i < len-1; i++) { // num_coefs[blocknum]
      l = tape[i] + ranges[blocknum][0];                // we are allowed to exceed the local range here
      r = tape[i+1] + ranges[blocknum][1];              // space is limited by ranges of first two coefs
      if (l < 0 || l > 2*ranges[blocknum][0]) continue;
      if (r < 0 || r > 2*ranges[blocknum][1]) continue;
      fc->vec->pairs[sl][qp_index][blocknum][l][r]++;
    }
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
//   long pos;
//   double sum_h, sum_p;
//   double scale_h, scale_p;
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
	    for (l = 0; l < 2*ranges[j][k]+1; l++) {
  //             printf("(%i, %i, %i, %i) -- (%i, %i)\n", i, j, k, l, count_h, count_p);
	      fc->vec->vector_histograms[count_h] = fc->vec->histograms[sl][i][j][k][l]; // (double) 
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
// 	      if (k == ranges[j][0] && l == ranges[j][1]) continue; // 00: &&; 0x x0: ||
  //             printf("(%i, %i)\n", count_h, count_p);
	      fc->vec->vector_pairs[count_p] = fc->vec->pairs[sl][i][j][k][l]; // (double) 
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
    if (sl == 0) if (fc->num_4x4_blocks_p == 0) continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!
    if (sl == 1) if (fc->num_4x4_blocks_b == 0) continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!
   
    // leave all pre-processing to Stegosaurus
    /*scale_h = ((double) fc->vec->vector_histograms_dim)/((double) N_h);
    for (i = 0; i < fc->vec->vector_histograms_dim; i++)
      fc->vec->vector_histograms[i] *= scale_h;
    scale_p = ((double) fc->vec->vector_pairs_dim)/((double) N_p);
    for (i = 0; i < fc->vec->vector_pairs_dim; i++)
      fc->vec->vector_pairs[i] *= scale_p;*/
    
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
      if (N_h != 0) {
	fwrite(fc->vec->vector_histograms, sizeof(int), fc->vec->vector_histograms_dim, current_hist_file); // double
      }
    }
//     printf("Stored histograms \n");
    if (current_pair_file != NULL) {
//       pos = ftell(current_pair_file);
//       if (pos == 0L)
// 	fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, current_pair_file);
      if (N_p != 0) {
	fwrite(fc->vec->vector_pairs, sizeof(int), fc->vec->vector_pairs_dim, current_pair_file); // double
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

FILE**** init_rate_bins(char pair, char *method_name, int dim) {
  int i, j;
  FILE ****result;
  char *blockstring;// = blockstrings[accept];
  char path[512];
  char *types[2] = {"hist", "pair"};
  char *type = types[pair];   // a bit tricky but should work if we only have histograms and pairs, maybe giving the string and then a hastable on int is better style
  double rate;
  char using_rate = 1;
  
  result = (FILE****) malloc(8*sizeof(FILE***)); // accept: [0..7]
  for (i = 0; i < 8; i++) {
    if (i != ACCEPT_L && i != ACCEPT_C && i != ACCEPT_LC) continue;
    blockstring = blockstrings[i];
//     printf("blockstring = %s \n", blockstring);
    result[i] = (FILE***) malloc(2*sizeof(FILE**)); // slice_type
    
    result[i][0] = (FILE**) malloc(NUM_BINS*sizeof(FILE*));
    for (j = 0; j < NUM_BINS; j++) {
//       printf("j=%i, %f \n", j, (((double) (j))/((double) NUM_BINS))*MAX_RATE);
      rate = (((double) j)/((double) NUM_BINS))*MAX_RATE;
      sprintf(path, "%s/%s/rate/p_%s/p_%s_rate_%i.fv", method_name, blockstring, type, type, (int) (10000. * rate));
//       printf("opening rate file: %s \n", path);
      result[i][0][j] = fopen(path, "a");
      writeHeader(result[i][0][j], pair, TYPE_P_SLICE, METHOD, using_rate, rate, i);
//       if (result[i][0][j] != NULL && ftell(result[i][0][j]) == 0L) 
//         fwrite(&dim, sizeof(int), 1, result[i][0][j]);
    }
    
    result[i][1] = (FILE**) malloc(NUM_BINS*sizeof(FILE*));
    for (j = 0; j < NUM_BINS; j++) {
      rate = (((double) j)/((double) NUM_BINS))*MAX_RATE;
      sprintf(path, "%s/%s/rate/b_%s/b_%s_rate_%i.fv", method_name, blockstring, type, type, (int) (10000. * rate));
//       result[i][1][j] = fopen(path, "a");
      result[i][1][j] = NULL;
      writeHeader(result[i][1][j], pair, TYPE_B_SLICE, METHOD, using_rate, rate, i);
//       if (result[i][1][j] != NULL && ftell(result[i][1][j]) == 0L) 
//         fwrite(&dim, sizeof(int), 1, result[i][1][j]);
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
// 	fflush(bins[i][0][j]);
	fclose(bins[i][0][j]);
      }
    }
//     printf("close: closed p files \n");
    free(bins[i][0]);
    
    for (j = 0; j < NUM_BINS; j++) {
      if (bins[i][1][j] != NULL) {
// 	fflush(bins[i][1][j]);
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
  int *****v;
  int *****w;
  int *****u;
  int hist_dim = 0;
  int pair_dim = 0;
  int *qp;
  int *tape;
  char method;
  char b_h_path[512];
  char b_p_path[512];
  char p_h_path[512];
  char p_p_path[512];
  char *blockstring = blockstrings[accept_blocks];//"L";
  
  v = malloc(2*sizeof(int****));
  w = malloc(2*sizeof(int****));
  u = malloc(2*sizeof(int****));
  for (sl = 0; sl < 2; sl++) {
    v[sl] = malloc(QP_RANGE*sizeof(int***));
    w[sl] = malloc(QP_RANGE*sizeof(int***));
    u[sl] = malloc(QP_RANGE*sizeof(int***));
    for (i = 0; i < QP_RANGE; i++) {
      v[sl][i] = malloc(3*sizeof(int**));   // 3 types of blocks (Luma + Chroma DC/AC)
      w[sl][i] = malloc(3*sizeof(int**));
      for (j = 0; j < 3; j++) {
        // setup histograms
        v[sl][i][j] = malloc(num_coefs[j]*sizeof(int*));
        for (k = 0; k < num_coefs[j]; k++) {
          v[sl][i][j][k] = malloc((2*ranges[j][k]+1)*sizeof(int));
        }
        // setup pairs
        w[sl][i][j] = malloc((2*ranges[j][0]+1)*sizeof(int*)); // +1 for the zero
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          w[sl][i][j][k] = malloc((2*ranges[j][1]+1)*sizeof(int));
        }
      }
      // setup uvsv
      u[sl][i] = malloc((num_coefs[2]+1)*sizeof(int**));
      u[sl][i][0] = malloc((2*ranges[1][0]+1)*sizeof(int*));
      for (k = 1; k < 2*ranges[1][0]+1; k++) {
        u[sl][i][0][k] = malloc((2*ranges[1][0]+1)*sizeof(int));
      }
      for (j = 1; j < num_coefs[2]+1; j++) { // there are 16 chroma coefs
        u[sl][i][j] = malloc((2*ranges[2][0]+1)*sizeof(int*));
        for (k = 1; k < 2*ranges[2][0]+1; k++) {
          u[sl][i][j][k] = malloc((2*ranges[2][0]+1)*sizeof(int));
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
  qp    =   malloc(50*sizeof(int));
  fc    =   malloc(sizeof(H264FeatureContext));
  fv    =   malloc(sizeof(H264FeatureVector));
  tape  =   malloc(16*sizeof(int));  
  
  fv->vector_histograms_dim = QP_RANGE*hist_dim;
  fv->vector_pairs_dim      = QP_RANGE*pair_dim;
  fv->vector_histograms     = malloc(fv->vector_histograms_dim*sizeof(double));
  fv->vector_pairs          = malloc(fv->vector_pairs_dim*sizeof(double));
  fv->histograms     = v;
  fv->pairs          = w;
  fv->uvsv           = u;
  fv->qp             = qp;
  fc->vec            = fv;
  fc->files_hist     = malloc(2*sizeof(FILE*));
  fc->files_pair     = malloc(2*sizeof(FILE*));
  fc->extract_rate   = extract_rate;
  if (p_hide < 0) {
    fc->files_hist[0]  = fopen("clean/p_histograms.fv", "a");
    fc->files_hist[1] = NULL;
    fc->files_pair[0]  = fopen("clean/p_pairs.fv", "a");
    fc->files_pair[1] = NULL;
  } else {
    if (extract_rate) {
      fc->rate_bins_hist = bins_h;
      fc->rate_bins_pair = bins_p;
      
      fc->files_hist[0] = NULL;
      fc->files_hist[1] = NULL;
      fc->files_pair[0] = NULL;
      fc->files_pair[1] = NULL;
    } else {
      sprintf(p_h_path, "%s/%s/prob/p_hist/p_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(p_p_path, "%s/%s/prob/p_pair/p_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(b_h_path, "%s/%s/prob/b_hist/b_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      sprintf(b_p_path, "%s/%s/prob/b_pair/b_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
      fc->files_hist[0]  = fopen(p_h_path, "a");
      fc->files_hist[1]  = fopen(b_h_path, "a");
      fc->files_pair[0]  = fopen(p_p_path, "a");
      fc->files_pair[1]  = fopen(b_p_path, "a");
    }
  }
  
  if (p_hide < 0) {
    method = 0;
  } else {
    method = METHOD;
  }
  writeHeader(fc->files_hist[0], 0, TYPE_P_SLICE, method, 0, p_hide, accept_blocks);
  writeHeader(fc->files_hist[1], 0, TYPE_B_SLICE, method, 0, p_hide, accept_blocks);
  writeHeader(fc->files_pair[0], 1, TYPE_P_SLICE, method, 0, p_hide, accept_blocks);
  writeHeader(fc->files_pair[1], 1, TYPE_B_SLICE, method, 0, p_hide, accept_blocks);
//   if (fc->files_hist[0] != NULL && ftell(fc->files_hist[0]) == 0L) 
//     fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_hist[0]);
//   if (fc->files_hist[1] != NULL && ftell(fc->files_hist[1]) == 0L) 
//     fwrite(&(fc->vec->vector_histograms_dim), sizeof(int), 1, fc->files_hist[1]);
//   if (fc->files_pair[0] != NULL && ftell(fc->files_pair[0]) == 0L) 
//     fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, fc->files_pair[0]);
//   if (fc->files_pair[1] != NULL && ftell(fc->files_pair[1]) == 0L) 
//     fwrite(&(fc->vec->vector_pairs_dim), sizeof(int), 1, fc->files_pair[1]);
  
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
  fc->lastU = (int*) malloc(num_coefs[1]*sizeof(int));
  fc->seenUs = (int*) malloc(num_coefs[2]*sizeof(int));
  fc->lastUs = (int**) malloc(num_coefs[2]*sizeof(int*));
  for (i = 0; i < num_coefs[2]; i++) {
    fc->lastUs[i] = (ranges[2][i]*sizeof(int));
  }

  return fc;
}

void close_features(H264FeatureContext* fc) {
  int i, j, k, sl;
  storeFeatureVectors(fc);
  
  if (fc->files_hist[0] != NULL) {
    fclose(fc->files_hist[0]);
  }
  if (fc->files_hist[1] != NULL) {
    fclose(fc->files_hist[1]);
  }
  if (fc->files_pair[0] != NULL) {
    fclose(fc->files_pair[0]);
  }
  if (fc->files_pair[1] != NULL) {
    fclose(fc->files_pair[1]);
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
      for (j = 1; j < num_coefs[2]+1; j++) { // there are 16 chroma coefs
        for (k = 1; k < 2*ranges[2][0]+1; k++) {
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
  
  free(fc->tape);
  free(fc->lastU);
  free(fc->vec->qp);
  free(fc->vec);
  free(fc);
}

void writeHeader(FILE *file, char pair, char slice_type, char method, char using_rate, double rate, char accept){
//   printf("Writing Header \n");
  if (file == NULL || ftell(file) != 0L) return;
//   printf("and not returning immediately \n");
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
  if (pair) {
    for (i = 0; i < 3; i++) {
      fwrite(&(ranges[i][0]), sizeof(unsigned char), 1, file);
      fwrite(&(ranges[i][1]), sizeof(unsigned char), 1, file);
    }
  } else {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < num_coefs[i]; j++) {
	fwrite(&(ranges[i][j]), sizeof(unsigned char), 1, file);
      }
    }
  }
}
