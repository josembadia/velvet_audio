/*
###########################################################################
#  Copyright 2023-24 Jose M. Badia <barrachi@uji.es> and                  #
#                    German Leon <leon@uji.es>                            #
#                                                                         #
#  audio_omp.c is free software: you can redistribute it and/or modify    #
#  it under the terms of the GNU General Public License as published by   #
#  the Free Software Foundation; either version 3 of the License, or      #
#  (at your option) any later version.                                    #
#                                                                         #
#  This file is distributed in the hope that it will be useful, but       #
#  WITHOUT ANY WARRANTY; without even the implied warranty of             #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      #
#  General Public License for more details.                               #
#                                                                         #
#  You should have received a copy of the GNU General Public License      #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>   #
#                                                                         #
###########################################################################
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//#define INT16
#define INT32
//#define FLOAT
#ifdef INT16
  typedef short btype;
#else
   #ifdef INT32
     typedef int btype;
   #else // FLOAT
     typedef float btype;
   #endif
#endif

//#define DEBUG
#define GOLDEN_TEST // Tests the result with the filter of size 88000
#define LAST_NB     // Filters the last block with nb samples

// Counts the number of lines of the text file
int countLines(char nombre[]) {
    FILE * fp;
    int nlines=0;
    char buffer[80];

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error opening file: %s\n", nombre);
        return 0;

    }
    else {
        fscanf(fp, "%hd", buffer);
        nlines++;
        while ( !feof(fp) ) {
           fscanf(fp, "%hd", buffer);
           nlines++;
        }
    }

    fclose(fp);
    return nlines;
}

// Reads a vector of type btype from the text file
int ReadVector(char nombre[], btype *v, int size) {
    char buffer[20];
    FILE * fp;
    int iReadValues=0;

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error opening file: %s\n", nombre);
        return 0;

    }
    else {
#ifdef INT16
        fscanf(fp, "%hd", &v[iReadValues++]);
#else
   #ifdef INT32
        fscanf(fp, "%d", &v[iReadValues++]);
   #else // FLOAT
        fscanf(fp, "%f", &v[iReadValues++]);
   #endif
#endif
        while ( !feof(fp) && (iReadValues < size) ) {
#ifdef INT16
           fscanf(fp, "%hd", &v[iReadValues++]);
#else
   #ifdef INT32
           fscanf(fp, "%d", &v[iReadValues++]);
   #else // FLOAT
           fscanf(fp, "%f", &v[iReadValues++]);
   #endif
#endif
        }
    }

    fclose(fp);
    return iReadValues-1;
}

// Reads the vector with the non-zero indices
int ReadVectorIndices(char nombre[], int v[]) {

    FILE * fp;
    int iReadValues=0;

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error opening file: %s\n", nombre);
        return 0;

    }
    else {
        fscanf(fp, "%d", &v[iReadValues++]);
        v[iReadValues-1]--;
        while ( !feof(fp) ) {
            fscanf(fp, "%d", &v[iReadValues++]);
            v[iReadValues-1]--;
        }
    }

    fclose(fp);
    return iReadValues-1;
}


// Generates a vector with nb random numbers
void GenerateRandomB(int *b, int nb) {
  srand(time(NULL));
  for (int i = 0; i < nb; i++) {
     if (rand() > RAND_MAX/2) b[i] = 1;
     else b[i] = -1;
  }
}

void GenerateRandomX(int *x, int nx) {
  srand(time(NULL));
  for (int i = 0; i < nx; i++) {
     x[i] = rand() % 100;
//     printf("x[%d] = %d\n", i, x[i]);
  }
}

// Compares filtered result with the golden result generated with matlab
int compareGolden(btype *v, btype *golden, int size) {
  int i = 0;
  while ( (i < size) && ( v[i] == golden[i] ) ) {
     i++;
  }

    if (i == size) return 1;
  else {
#ifdef FLOAT
     printf("i: %d, v: %f, golden: %f\n", i, v[i], golden[i]);
#else
     printf("i: %d, v: %d, golden: %d\n", i, v[i], golden[i]);
#endif
     return 0;
  }
}

void buildVectorDense(btype *b, int nb, int *ind, int nind, btype *bdense) {
   for (int i = 0; i < nind; i++)
      bdense[ind[i]] = b[i];
}

// Version 1. Conditional filter (COF)
// For each sample loop over the non-zero values of b
// If b=+1 add sample to the element of y, if b=-1, substract sample. Avoids product
void processAudio_v1(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
  btype *ptr = x;
  btype elem;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+ind[nind-1];
#else
  int end = nsamples;
#endif
  btype aux;

  for (int is = 0; is < end; is += bsize) {  
     for (iy = is; iy < MIN(is+bsize, end); iy++) {  
        aux = 0;
        ib = 0;
        while ( (ind[ib] <= iy) && (ib < nind) )  { 

          elem = *(ptr - ind[ib]);
          if (b[ib] > 0)
             aux += elem;
          else
             aux -= elem;

          ib++;
       }
       y[iy] = aux;

       ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 2. Sparse multiplier filter (SMF)
// For each sample loop over the non-zero values of b
void processAudio_v2(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
  btype *ptr = x;
  btype elem, aux;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+ind[nind-1];
#else
  int end = nsamples;
#endif

  for (int is = 0; is < end; is += bsize) {  
     for (iy = is; iy < MIN(is+bsize, end); iy++) {  
	aux = 0;
        ib = 0;
        while ( (ind[ib] <= iy) && (ib < nind) )  { 


           elem = *(ptr - ind[ib]);
           aux += elem*b[ib];

           ib++;
       }
       y[iy] = aux;
       ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 3. Double-vector filter (DVF)
// For each sample loop over the +1 of b and then over the -1 of b
// Avoids 'if' in the inner loop and deals first with the positive elements of b and next with the negative elements.
// Performs only additions and substractions as in version v0.
void processAudio_v3(btype *x, int nsamples, int nb, int *indp, int np, int *indn, int nn, btype *y, int bsize) {
  btype *ptr = x;
  btype  aux;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+MAX(indp[np-1], indn[nn-1]);
#else
  int end = nsamples;
#endif


  for (int is = 0; is < end; is += bsize) {  
     for (int iy = is; iy < MIN(is+bsize, end); iy++) {  
   
        // First add elements of x associated to elements of b=+1
        aux = 0;
        ib = 0;
        while ( (indp[ib] <= iy) && (ib < np) )  { 
           aux += *(ptr - indp[ib]);


           ib++;
        }
        // Next add elements of x associated to elements of b=-1
        ib = 0;
        while ( (indn[ib] <= iy) && (ib < nn) )  { 

            aux -= *(ptr - indn[ib]);

            ib++;
        }
        y[iy] = aux;

        ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 4. Transposed double-vector filter (DVF Transpose)
// Inverts the two main loops reducing the number of memory accesses and better leveraging the cache memory
// For each positive element of b loops over the elements of y and add the corresponding sample.
// Next, for each negative element of b loops over the elements of y and substracts the corresponding sampmle
void processAudio_v4(btype *x, int nsamples, int nb, int *indp, int np, int *indn, int nn, btype *y, int bsize) {
   int is, ib, off, ini, end, endb, ith, iy;

   // Process first block of nb samples increasing the number of samples processed on each iteration
   for (is = indp[0]; is <= nb; is += bsize) {
      for (ib = 0; ib < np; ib++) {
         off = indp[ib];
	 ini = MAX(is, off);
	 end = MIN(is+bsize, nb);
         for (int iy = ini; iy < end; iy++) {
            y[iy] += x[iy-off];
         }
      }
   } 

   for (is = indn[0]; is <= nb; is += bsize) {
      for (ib = 0; ib < nn; ib++) {
         off = indn[ib];
	 ini = MAX(is, off);
	 end = MIN(is+bsize, nb);
         for (iy = ini; iy < end; iy++) {
            y[iy] -= x[iy-off];
         }
      } 
   } 

   for (int iblock = nb; iblock < nsamples; iblock += nb) {

      end = MIN(iblock+nb, nsamples);

      for (is = iblock; is <= end; is += bsize) {
         for (ib = 0; ib < np; ib++) {
            endb = MIN(is+bsize, end);
            off = indp[ib];
            for (iy = is; iy < endb; iy++) {
               y[iy] += x[iy-off];
           }
        }
      } 

      for (is = iblock; is <= end; is += bsize) {
         for (ib = 0; ib < nn; ib++) {
            endb = MIN(is+bsize, end);
            off = indn[ib];
            for (iy = is; iy < endb; iy++) {
               y[iy] -= x[iy-off];
            }
         }
      }

   } // end for iblock

#ifdef LAST_NB
   // Process the last block of samples until the pipeline is empty
   for (is = nsamples; is < nsamples+indp[np-1]; is += bsize) {
      for (ib = 0; ib < np; ib++) {
         off = indp[ib];
	 endb = MIN(is+bsize, nsamples+off);
         for (iy = is; iy < endb; iy++) {
            y[iy] += x[iy-off];
         }
      }
   } 

   for (is = nsamples; is < nsamples+indn[nn-1]; is += bsize) {
      for (ib = 0; ib < nn; ib++) {
         off = indn[ib];
	 endb = MIN(is+bsize, nsamples+off);
         for (iy = is; iy < endb; iy++) {
            y[iy] -= x[iy-off];
         }
      }
   }
#endif

}

// Version 5. Transposed sparse multiplier filter (SMF Tranpose)
// As in version 4, but avoids looping twice over the indices and elements of y by performing and additional product as in v1
void processAudio_v5(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
   int off, ith, iy, endb, end;
   btype elemb;

   // Process first block of nb samples increasing the number of samples processed on each iteration
   for (int is = ind[0]; is <= nb; is += bsize ) {
      for (int ib = 0; ib < nind; ib++) {
         endb = MIN(is+bsize, nb);
         off = ind[ib];
         elemb = b[ib];
         for (int iy = MAX(is, ind[ib]); iy < endb; iy++) {
            y[iy] += x[iy-off]*elemb;
         }
      } 
   } // end is

   for (int iblock = nb; iblock < nsamples; iblock += nb) {

      end = MIN(iblock+nb, nsamples);
      for (int is = iblock; is <= end; is += bsize ) {
         for (int ib = 0 ; ib < nind ; ib++) {
	    endb = MIN(is+bsize, end);
            off = ind[ib];
            elemb = b[ib];
            for (int iy = is; iy < endb; iy++) {
               y[iy] += x[iy-off]*elemb;
           }
         } 
     } // end is

   } // end for iblock

#ifdef LAST_NB
   // Process the last block of samples until the pipeline is empty
   for (int is = nsamples; is <= nsamples+ind[nind-1]; is += bsize ) {
      for (int ib = 0 ; ib < nind ; ib++) {
         off = ind[ib];
         elemb = b[ib];
	 endb = MIN(is+bsize, nsamples+off);
         for (int iy = is; iy < endb; iy++) {
            y[iy] += x[iy-off]*elemb;
         }
      } 
   } // end is
#endif

}

int main(int argc, char *argv[]) {

  const int MAXINDICES = 100000;
  int nb = 88000;    // size of vector b containing the coeficients of the long filter
  int nind;          // number of non_zero elements in b
  int nsamples, nres;

  btype **x;        // Sample vector
  btype **y;        // Result vector with the filtered samples
  btype *b = (btype *) malloc(MAXINDICES*sizeof(btype));   // Store only non-zero {-1,+1} values
#ifdef GOLDEN_TEST
  btype *golden;   // resultado en matlab
#endif
  int *ind;  // indexes of all the non-zero elements of b
  int *indp; // indexes of the +1 elements of b
  int *indn; // indexes of the -1 elements of b

  int i, ch;

  int ntimes = 3; // Number of times that we execute the filtering process

  int np, nn;
  int ct;
  int nth;       // number of threads in the OpenMP parallel version
  int nchannels; // number of channels
  int nblocks;   // number of blocks
                     // -1: Read samples from file
		     // 0: Compute 10*nb samples
		     // 1: Compute only the first block of samples
  int bsize;     // block size
  char *type;    // data type: s (short), i (int), f (float)
  int version; // Filter version: 1, 2, 3, 4 or 5

  if (argc < 7) {
     printf("FORMAT: audio_omp <nth> <nch> <nblocks> <bsize> <nb> <type> <version>\n");
     return 1;
  } else {
     nth = atoi(argv[1]);
     nchannels = atoi(argv[2]);
     nblocks = atoi(argv[3]);
     bsize = atoi(argv[4]);
     nb = atoi(argv[5]);
     type = argv[6];
     version = atoi(argv[7]);
  }
#ifdef DEBUG
  printf("Number of threads: %d\n", nth); fflush(stdout);
#endif


  // Build files names from paramenters
  char xname[40];
  char bname[40]; // = "../data/b_88000_100.txt";
  char goldename[40] = "../data/goldenOnes.txt";
  char indname[40]; // = "../data/ind_88000_100.txt";
  char datadir[20] = "../data/";
  sprintf(xname, "%sx_%s.txt", datadir, type);
#ifdef DEBUG
  printf("xname: %s\n", xname);
#endif
  sprintf(bname, "%sb_%d.txt", datadir, nb);
#ifdef DEBUG
  printf("bname: %s\n", bname);
#endif
  sprintf(goldename, "%sgolden_%d_%s.txt", datadir, nb, type);
#ifdef DEBUG
  printf("goldename: %s\n", goldename);
#endif
  sprintf(indname, "%sind_%d.txt", datadir, nb);
#ifdef DEBUG
  printf("indname: %s\n", indname);
#endif


  nind = ReadVector(bname, b, MAXINDICES);
#ifdef DEBUG
  printf("Read %d values of b\n", nind); fflush(stdout);
#endif

  ind = (int *) malloc(nind*sizeof(int));
  indp = (int *) malloc(nind*sizeof(int));
  indn = (int *) malloc(nind*sizeof(int));

  ct = ReadVectorIndices(indname, ind);
#ifdef DEBUG
  printf("Read %d indices\n", ct); fflush(stdout);
#endif

  if (nblocks == -1) {  // Read samples from file
     nsamples = countLines(xname);
#ifdef DEBUG
     printf("Number of samples on file: %d\n", nsamples); fflush(stdout);
#endif
  }
  else if (nblocks == 1) { // Compute only the cost of processing the first block of samples
     nsamples = nb = bsize;
     ntimes = 1000;
  }
  else {
//     nsamples = MAX(2*nb, nblocks*bsize); // To compute the bandwidth
     nsamples = 10*nb;
     nblocks = MAX(nsamples/bsize, 1);
  }

//  bdense = (btype *) calloc(nb, sizeof(btype) );
//  buildVectorDense(b, nb, ind, nind, bdense);
//
  nres = nsamples+nb-1;

  // Allocate space for the samples
  x = (btype **) malloc(nchannels * sizeof(btype *));
  for (int ch = 0; ch < nchannels; ch++)
      x[ch] = (btype *) calloc(nres, sizeof(btype) );

  if (nblocks == -1) { // Read samples from file
     nsamples = ReadVector(xname, x[0], nsamples);
#ifdef DEBUG
     printf("Read %d samples.\n", nsamples); fflush(stdout);
#endif
     nblocks = nsamples/bsize;
  }
  else { // Init the samples to one, except the last nb-size block with zeros
     for (int i = 0; i < nsamples; i++)
        x[0][i] = 1;
  }

  // Copy the first channel to the rest
  for (int i = 1; i < nchannels; i++)
      memcpy(x[i], x[0], nsamples*sizeof(btype));


#ifdef GOLDEN_TEST
  golden = (btype *) malloc(nres*sizeof(btype));
  ct = ReadVector(goldename, golden, nres);
//  ct = ReadVector("../data/goldenOnes.txt", golden, nres);
#endif

  // Allocate space for the results
  y = (btype **) malloc(nchannels * sizeof(btype *));
  for (int i = 0; i < nchannels; i++)
      y[i] = (btype *) calloc(nres, sizeof(btype) );

  // Separate indexes of elements +1 and -1 of vector b
  np = nn = 0;
  for (int i = 0; i < nind; i++) { 
    ind[i]; // To convert matlab (from 1) to C (from 0) indexing
    if (b[i] > 0) { 
       indp[np++] = ind[i];
    }
    else { 
       indn[nn++] = ind[i];
    }
  }
#ifdef GOLDEN_TEST
  printf("Number of +1: %d; Number of -1: %d\n", np, nn);
#endif

  double start, end;
  omp_set_num_threads(nth);

//  for (ch = 0; ch < nchannels; ch++)
//     memset(y[ch], 0, nres*sizeof(btype));

  start = omp_get_wtime();

  switch (version) {
    case 1: 
      for (i = 0; i < ntimes; i++)
        #pragma omp parallel for num_threads(nth)
        for (ch = 0; ch < nchannels; ch++)
          processAudio_v1(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);

      break;

    case 2:
      for (i = 0; i < ntimes; i++)
        #pragma omp parallel for num_threads(nth)
        for (ch = 0; ch < nchannels; ch++)
          processAudio_v2(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);

    case 3:
      for (i = 0; i < ntimes; i++)
        #pragma omp parallel for num_threads(nth)
        for (ch = 0; ch < nchannels; ch++)
          processAudio_v3(x[ch], nsamples, nb, indp, np, indn, nn, y[ch], bsize);

      break;

    case 4:
      for (i = 0; i < ntimes; i++)
        #pragma omp parallel for num_threads(nth)
        for (ch = 0; ch < nchannels; ch++)
          processAudio_v4(x[ch], nsamples, nb, indp, np, indn, nn, y[ch], bsize);

      break;

    case 5:
      for (i = 0; i < ntimes; i++)
        #pragma omp parallel for num_threads(nth)
        for (ch = 0; ch < nchannels; ch++)
          processAudio_v5(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);
    } // end switch

  end = omp_get_wtime();

  printf("v%d #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
        version, nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );

#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif


  // Free memory
  for (ch = 0; ch < nchannels; ch++) {
     free(x[ch]);
     free(y[ch]);
  }
  free(x);
  free(y);
  free(b);
  free(ind);
  free(indp);
  free(indn);
#ifdef GOLDEN_TEST
  free(golden);
#endif

}
