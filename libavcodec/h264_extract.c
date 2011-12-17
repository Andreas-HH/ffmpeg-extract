#include "h264_extract.h"

/*
 (Luma DC:   48)
 Luma AC:   0 - 15
 Chroma DC: 49, 50
 Chroma AC: 16 - 19, 32 - 35
 */


void perform_hiding_plusminus(H264FeatureContext *fc) {
  int i;
//   H264FeatureContext *fc = (H264FeatureContext*) vfc;
  int *coefs = fc->tape;
  double r;
  
  if (fc->blocknum == -1) return;
  
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
  addCounts(fc, fc->current_qp, fc->blocknum);
//   pthread_exit(NULL);
}

void *PrintHello(void *threadid)
{
   long tid = 5;
   printf("Hello World! It's me, thread #%ld!\n", tid);
   pthread_exit(NULL);
}

void simulate_hiding_plusminus(H264FeatureContext *fc) {
  int rc;

  rc = pthread_create(fc->thread, fc->thread_attr, (void *) &perform_hiding_plusminus, (void *) fc);
  if (rc) printf("ERROR: pthread_create failed! %i \n", rc);
}

void wait_for_simulation(H264FeatureContext* fc) {
  int rc;
  void *thread_status;
  
  rc = pthread_join(*(fc->thread), &thread_status);
  if (rc) printf("ERROR: pthread_join failed: %i \n", rc);
}


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

void addCounts(H264FeatureContext *fc, int qp, int blocknum) {
  int i, l, r;
  int coef_index;
  int *tape = fc->tape;
//   int blocknum = get_block_index(n);
  int qp_index = qp - QP_OFFSET;
  int sl = fc->slice_type;
  
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
  
  fc->vec->qp[qp]++;
  if (blocknum == -1) {
//     fprintf(fc->file, "blocknum is -1 o_O %i \n", n);
    return;
  }
  if (qp_index < 0 || qp_index >= QP_RANGE) {
//     fprintf(fc->file, "qp out of range o_O %i \n", qp);
    return;
  }
  
//   printf("extracted qp histogram %i, %i \n", sl, blocknum);
  
//   constructProperCoefArray(tape, level, run_before, total_coeff, totalZeros, blocknum, fc->vec);
  for (i = 0; i < num_coefs[blocknum]; i++) {
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
    fc->vec->histograms[sl][qp_index][blocknum][i][coef_index]++;
  }
//   printf("extracted hist ");
  //pairs
  for (i = 0; i < num_coefs[blocknum]-1; i++) {
    l = tape[i] + ranges[blocknum][0];                // we are allowed to exceed the local range here
    r = tape[i+1] + ranges[blocknum][1];              // space is limited by ranges of first two coefs
    if (l < 0 || l > 2*ranges[blocknum][0]) continue;
    if (r < 0 || r > 2*ranges[blocknum][1]) continue;
    fc->vec->pairs[sl][qp_index][blocknum][l][r]++;
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
  int sl, i, j, k, l;
  int count_h, count_p;
  int N_h, N_p;
  long pos;
  double sum_h, sum_p;
  double scale_h, scale_p;
  
//   printf("Storing feature vectors.. ");
  
  if (!fc->refreshed) return;
  
//   printf("I am fresh! \n");
  
  for (sl = 0; sl < 2; sl++) {
    count_h = 0;
    count_p = 0;
    N_h = 0;
    N_p = 0;
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
          for (l = 0; l < 2*ranges[j][k]; l++) {
//             printf("(%i, %i, %i, %i) -- (%i, %i)\n", i, j, k, l, count_h, count_p);
            fc->vec->vector_histograms[count_h] = (double) fc->vec->histograms[sl][i][j][k][l];
            count_h++;
            N_h += fc->vec->histograms[sl][i][j][k][l];
//             fprintf(fc->files_hist[sl], "%i ", fc->vec->histograms[sl][i][j][k][l]);
          }
        }
        // pairs
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
      }
    }
//     printf("(%i, %i)\n", count_h, count_p);
    
    scale_h = ((double) fc->vec->vector_histograms_dim)/((double) N_h);
    for (i = 0; i < fc->vec->vector_histograms_dim; i++)
      fc->vec->vector_histograms[i] *= scale_h;
    scale_p = ((double) fc->vec->vector_pairs_dim)/((double) N_p);
    for (i = 0; i < fc->vec->vector_pairs_dim; i++)
      fc->vec->vector_pairs[i] *= scale_p;
  
    // check if vector sums orrectly, for debugging only
    /*sum_h = 0.;
    for (i = 0; i < fc->vec->vector_histograms_dim; i++) {
      sum_h += fc->vec->vector_histograms[i];
    }
    sum_p = 0.;
    for (i = 0; i < fc->vec->vector_pairs_dim; i++) {
      sum_p += fc->vec->vector_histograms[i];
    }*/
//     printf("dim_h = %i, dim_p = %i, sum_h = %f, sum_p = %f,  N_h = %i, N_p = %i \n", fc->vec->vector_histograms_dim, fc->vec->vector_pairs_dim, sum_h, sum_p, N_h, N_p);
    if (N_h == 0) continue;  // ffmpeg does some test-decoding in the beginning. We dont want to count this!
    
    printf("b-rate: %f, p-rate: %f \n", ((double)fc->hidden_bits_b)/((double)fc->num_4x4_blocks_b), ((double)fc->hidden_bits_p)/((double)fc->num_4x4_blocks_p));
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
  feature_context->num_4x4_blocks_b = 0;
  feature_context->num_4x4_blocks_p = 0;
  feature_context->hidden_bits_b = 0;
  feature_context->hidden_bits_p = 0;

  for (sl = 0; sl < 2; sl++) {
    for (i = 0; i < QP_RANGE; i++) {
      for (j = 0; j < 3; j++) {
        // histograms
        for (k = 0; k < num_coefs[j]; k++) {
//           feature_context->vec->N[i][j][k] = 0;
          for (l = 0; l < 2*ranges[j][k]; l++) {
            feature_context->vec->histograms[sl][i][j][k][l] = 0;
          }
        }
        // pairs
        for (k = 0; k < 2*ranges[j][0]+1; k++) {
          for (l = 0; l < 2*ranges[j][1]+1; l++) {
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

H264FeatureContext* init_features(char *method_name, int accept_blocks, double p_hide) {
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
  char *blockstring = "L";
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
  fc->thread_attr    = (pthread_attr_t*) malloc(sizeof(pthread_attr_t));
  fc->thread         = (pthread_t*) malloc(sizeof(pthread_t));
  fc->logfile        = fopen("features.log", "w"); // fclose(file)
  fc->files_hist     = malloc(2*sizeof(FILE*));
  fc->files_pair     = malloc(2*sizeof(FILE*));
  if (p_hide < 0) {
    fc->files_hist[0]  = fopen("clean/p_histograms.fv", "a");
    fc->files_hist[1]  = fopen("clean/b_histograms.fv", "a");
    fc->files_pair[0]  = fopen("clean/p_pairs.fv", "a");
    fc->files_pair[1]  = fopen("clean/b_pairs.fv", "a");
  } else {
//     printf("opening files ");
    sprintf(p_h_path, "%s/%s/prob/p_hist/p_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
    sprintf(p_p_path, "%s/%s/prob/p_pair/p_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
    sprintf(b_h_path, "%s/%s/prob/b_hist/b_hist_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
    sprintf(b_p_path, "%s/%s/prob/b_pair/b_pair_prob_%i.fv", method_name, blockstring, (int) (p_hide*1000.));
//     printf("%s", p_h_path);
    fc->files_hist[0]  = fopen(p_h_path, "a");
    fc->files_hist[1]  = fopen(b_h_path, "a");
    fc->files_pair[0]  = fopen(p_p_path, "a");
    fc->files_pair[1]  = fopen(b_p_path, "a");
  }
  
//   fc->p_hist = fopen("p_histograms.fv", "a");
//   fc->p_pair = fopen("p_pairs.fv", "a");
//   fc->b_hist = fopen("b_histograms.fv", "a");
//   fc->b_pair = fopen("b_pairs.fv", "a");
  fc->tape = tape;
//   fc->histogram_ranges = ranges;
  fc->accept_blocks = accept_blocks;
  fc->p_hide = p_hide;
  
  fc->i_slices    = 0;
  fc->p_slices    = 0;
  fc->b_slices    = 0;
  fv->vector_num  = 0;
  fc->refreshed   = 0;
  
//   printf("Done. \n");
  fc->blocknum = -1;
  pthread_attr_init(fc->thread_attr);
  pthread_attr_setdetachstate(fc->thread_attr, PTHREAD_CREATE_JOINABLE);
  
  return fc;
//   h->feature_context = fc;
}

void close_features(H264FeatureContext* fc) {
  int i, j, k, sl;
//   H264FeatureContext *fc = h->feature_context;
//   printf("Closing features... ");
//   storeCounts(fc);
  storeFeatureVectors(fc);
  
  fflush(fc->logfile);
  fflush(fc->files_hist[0]);
  fflush(fc->files_hist[1]);
  fflush(fc->files_pair[0]);
  fflush(fc->files_pair[1]);
  fclose(fc->logfile);
  fclose(fc->files_hist[0]);
  fclose(fc->files_hist[1]);
  fclose(fc->files_pair[0]);
  fclose(fc->files_pair[1]);

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
  free(fc->vec->histograms);
  free(fc->vec->pairs);
  free(fc->vec->vector_histograms);
  free(fc->vec->vector_pairs);
  
//   free(fc->histogram_ranges[0]);
//   free(fc->histogram_ranges[1]);
//   free(fc->histogram_ranges[2]);
//   free(fc->histogram_ranges);
  
  free(fc->thread_attr);
  free(fc->thread);
  free(fc->tape);
//   av_free(fc->vec->mb_t);
  free(fc->vec->qp);
  free(fc->vec);
  free(fc);
}

