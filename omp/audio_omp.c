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

//#define GOLDEN_TEST // Determina si se comprueba el resultado obtenido
//#define LAST_NB     // Determina si se computa la salidad el último bloque de nb muestras

int OpenWav(unsigned char **pwav, char nombre[])
{

  FILE *wav_fp;
  int iReadMues=0;
  unsigned char *speechdat;                     //Usada para leer los bytes del fichero
  unsigned char *header;                                //Cabecera de fichero
        //unsigned char *pInit;                         //puntero al principio del array
  int iSize;                                                    //Numero de bytes leidos del wav.
  int numread;                                          //Elementos leidos
  int j=0;

        if ((wav_fp = fopen(nombre, "rb")) == NULL)
        {
            printf("Can't open directory listing file '%s'\n", nombre);
        exit(0);
        }
    // obtain file size:
  fseek (wav_fp , 0 , SEEK_END);
  iSize = ftell (wav_fp);
  rewind (wav_fp);
  //Then I copy the first 44 bytes of the wav file (i.e. header) to a structure:
  header = (unsigned char *) calloc(44,1);
  fread(header, 1, 44, wav_fp);
  //Now, I write the samples (i.e. the data after the header) into an array:
  iSize=iSize-44; //le restamos la cabecera
  //Reservamos memoria para las muestras musicales del archivo wav
  speechdat = (unsigned char *) calloc(iSize,1);
  //leemos muestras del fichero
  numread = fread(speechdat, 1, iSize, wav_fp);

  *pwav = speechdat;

   fclose(wav_fp);      //Cerramos el fichero wav abierto antes
   free(header);
   // -- free(speechdat);


  return iSize;

}

int OpenWavConvert32(unsigned char **pwav, char nombre[])
{

   int j=0;
   int iSize=0;
   unsigned char * speechdat;
   unsigned char * pwavTmp;

   iSize = OpenWav(&speechdat, nombre);
   pwavTmp = (unsigned char *) calloc(iSize*4,1);
   //outwav[iorden] = (unsigned char *) calloc(iSize*4*2,1);
   memset(pwavTmp,0,iSize*4);
   //memset(outwav[iorden],0,iSize*4*2);
   //Ahora leemos el wav y ponemos los bytes en el orden que toca
   for(j=0;j<iSize;j=j+2)
   {
             pwavTmp[2*j+2]=speechdat[j];
       pwavTmp[2*j+3]=speechdat[j+1];
   }

   *pwav = pwavTmp;

   free(speechdat);
   iSize = iSize/2;
   return iSize;

}

int countLines(char nombre[]) {
    FILE * fp;
    int nlines=0;
    char buffer[80];

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error al abrir el archivo \n");
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

int ReadVector(char nombre[], btype *v, int size) {
    char buffer[20];
    FILE * fp;
    int iReadValues=0;
    // Lectura de los coeficientes del filtro, h1, es un camino con 220 coeficientes.
//    printf("%s\n",nombre);

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error al abrir el archivo \n");
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

int ReadVectorIndices(char nombre[], int v[]) {

    FILE * fp;
    int iReadValues=0;
    // Lectura de los coeficientes del filtro, h1, es un camino con 220 coeficientes.
//    printf("%s\n",nombre);

    fp = fopen(nombre, "r" );
    if (fp==NULL) {
        printf("Error al abrir el archivo \n");
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


void GenerateRandomB(int *b, int nb) {
  srand(time(NULL));
  for (int i = 0; i < nb; i++) {
     if (rand() > RAND_MAX/2) b[i] = 1;
     else b[i] = -1;
//     printf("b[%d] = %d\n", i, b[i]);
  }
}

void GenerateRandomX(int *x, int nx) {
  srand(time(NULL));
  for (int i = 0; i < nx; i++) {
     x[i] = rand() % 100;
//     printf("x[%d] = %d\n", i, x[i]);
  }
}

void PrintVectorInt(int *v, int size, char *name) {
   for (int i = 0; i < size; i++)
      printf("%s[%d] = %d\n", name, i, v[i]);
}

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

void  buildVectorDense(btype *b, int nb, int *ind, int nind, btype *bdense) {
   for (int i = 0; i < nind; i++)
      bdense[ind[i]] = b[i];
}

// Version 0
// For each sample loop over the non-zero values of b
// If b=+1 add sample to the element of y, if b=-1, substract sample. Avoids product
void processAudio_v0(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
  btype *ptr = x;
  btype elem;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+ind[nind-1];
#else
  int end = nsamples;
#endif
  btype aux;

  // Cada iteracion de este bucle es independiente y la podria hacer una hebra distinta.
//  #pragma omp parallel for schedule(static, 16) firstprivate(ptr) private(ib)
//  for (int is = 0; is < nsamples+nb-1; is += bsize) {  
  for (int is = 0; is < end; is += bsize) {  
//     printf("is: %d\n", is);
     for (iy = is; iy < MIN(is+bsize, end); iy++) {  
        aux = 0;
//        y[iy] = 0;
        ib = 0;
        while ( (ind[ib] <= iy) && (ib < nind) )  { 

//      printf("ind[%d] = %d, b[%d] = %d, %d\n", ib, ind[ib], idx, b[idx], is ); fflush(stdout);

          elem = *(ptr - ind[ib]);
          if (b[ib] > 0)
             aux += elem;
//             y[iy] += elem;
          else
             aux -= elem;
//             y[iy] -= elem;
//          (b[ib] > 0) ? (y[iy] += elem) : (y[iy] -= elem);
//      printf("y[%d] += %d --> %d\n", is, *(ptr - indp[ib]), y[is]); fflush(stdout);

          ib++;
       }
       y[iy] = aux;

       ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 1
// For each sample loop over the non-zero values of b
void processAudio_v1(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
  btype *ptr = x;
  btype elem, aux;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+ind[nind-1];
#else
  int end = nsamples;
#endif

  // Cada iteracion de este bucle es independiente y la podria hacer una hebra distinta.
//  #pragma omp parallel for schedule(static, 16) firstprivate(ptr) private(ib)
//  for (int is = 0; is < nsamples+nb-1; is += bsize) {  
  for (int is = 0; is < end; is += bsize) {  
//     printf("is: %d\n", is);
     for (iy = is; iy < MIN(is+bsize, end); iy++) {  
	aux = 0;
//        y[iy] = 0;
        ib = 0;
        while ( (ind[ib] <= iy) && (ib < nind) )  { 

//      printf("ind[%d] = %d, b[%d] = %d, %d\n", ib, ind[ib], idx, b[idx], is ); fflush(stdout);

           elem = *(ptr - ind[ib]);
           aux += elem*b[ib];
//           y[iy] += elem*b[ib];
//      printf("y[%d] += %d --> %d\n", is, *(ptr - indp[ib]), y[is]); fflush(stdout);

           ib++;
       }
       y[iy] = aux;
       ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 2
// For each sample loop over the +1 of b and then over the -1 of b
// Avoids 'if' in the inner loop and deals first with the positive elements of b and next with the negative elements.
// Performs only additions and substractions as in version v0.
void processAudio_v2(btype *x, int nsamples, int nb, int *indp, int np, int *indn, int nn, btype *y, int bsize) {
  btype *ptr = x;
  btype  aux;
  int ib, iy;
#ifdef LAST_NB
  int end = nsamples+MAX(indp[np-1], indn[nn-1]);
#else
  int end = nsamples;
#endif


  // Cada iteracion de este bucle es independiente y la podria hacer una hebra
  // distinta.
//  #pragma omp parallel for schedule(static, 16) firstprivate(ptr) private(ib)
  for (int is = 0; is < end; is += bsize) {  
//     printf("is: %d\n", is);
     for (int iy = is; iy < MIN(is+bsize, end); iy++) {  
   
        // First add elements of x associated to elements of b=+1
        aux = 0;
//        y[iy] = 0;
        ib = 0;
        while ( (indp[ib] <= iy) && (ib < np) )  { 
           aux += *(ptr - indp[ib]);
//           y[iy] += *(ptr - indp[ib]);
//      printf("y[%d] += %d --> %d\n", is, *(ptr - indp[ib]), y[is]); fflush(stdout);

//      printf("indp[%d] = %d, y[%d] += %d\n", ib,  indp[ib], is, y[is] ); fflush(stdout);

           ib++;
        }
        // Next add elements of x associated to elements of b=-1
        ib = 0;
        while ( (indn[ib] <= iy) && (ib < nn) )  { 

            aux -= *(ptr - indn[ib]);
//            y[iy] -= *(ptr - indn[ib]);
//      printf("y[%d] -= %d --> %d\n", is, *(ptr - indn[ib]), y[is]); fflush(stdout);

//      printf("indn[%d] = %d, y[%d] -= %d\n", ib,  indn[ib], is, y[is] ); fflush(stdout);

            ib++;
        }
        y[iy] = aux;

//    printf("%d\n", y[is]); fflush(stdout);

        ptr++; // Move to the next sample

     } // end for iy

  } // end for is

}

// Version 3
// Inverst the two main loops reducing the number of memory accesses and better leveraging the cache memory
// For each positive element of b loops over the elements of y and add the corresponding sample.
// Next, for each negative element of b loops over the elements of y and substracts the corresponding sampmle
// Parallelism over y. 
// For each non-zero element of b we can update in parallel all the elements of y with the corresponding sample
void processAudio_v3(btype *x, int nsamples, int nb, int *indp, int np, int *indn, int nn, btype *y, int bsize) {
   int is, ib, off, ini, end, endb, ith, iy;

//#pragma omp parallel private(is, ib, off, ini, end, endb, ith, iy)
//{

//   printf("start: %d, end: %d\n", 0, nb);
   // Process first block of nb samples increasing the number of samples processed on each iteration
   for (is = indp[0]; is <= nb; is += bsize) {
      for (ib = 0; ib < np; ib++) {
         off = indp[ib];
	 ini = MAX(is, off);
	 end = MIN(is+bsize, nb);
//#pragma omp for
         for (int iy = ini; iy < end; iy++) {
            y[iy] += x[iy-off];
//	 printf("[%d]: y[%d] += %d --> %d\n", ith, is, x[is],  y[is]); fflush(stdout);
         }
      }
   } 
//   printf("-----------\n"); fflush(stdout);

   for (is = indn[0]; is <= nb; is += bsize) {
      for (ib = 0; ib < nn; ib++) {
         off = indn[ib];
	 ini = MAX(is, off);
	 end = MIN(is+bsize, nb);
//#pragma omp for
         for (iy = ini; iy < end; iy++) {
            y[iy] -= x[iy-off];
//	 printf("y[%d] -= %d --> %d\n", is, *ptrx, y[is]); fflush(stdout);
         }
      } 
   } 
//   printf("-----------\n"); fflush(stdout);

   for (int iblock = nb; iblock < nsamples; iblock += nb) {

      end = MIN(iblock+nb, nsamples);
//      printf("start: %d, end: %d\n", iblock, end);

      for (is = iblock; is <= end; is += bsize) {
         for (ib = 0; ib < np; ib++) {
            endb = MIN(is+bsize, end);
            off = indp[ib];
//#pragma omp for
            for (iy = is; iy < endb; iy++) {
               y[iy] += x[iy-off];
           }
        }
      } 
//      printf("-----------\n"); fflush(stdout);

      for (is = iblock; is <= end; is += bsize) {
         for (ib = 0; ib < nn; ib++) {
            endb = MIN(is+bsize, end);
            off = indn[ib];
//#pragma omp for
            for (iy = is; iy < endb; iy++) {
               y[iy] -= x[iy-off];
            }
         }
      }
//      printf("-----------\n"); fflush(stdout);

   } // end for iblock

#ifdef LAST_NB
   // Process the last block of samples until the pipeline is empty
//   printf("last start: %d, end: %d\n", nsamples, nsamples+nb-1);
   for (is = nsamples; is < nsamples+indp[np-1]; is += bsize) {
      for (ib = 0; ib < np; ib++) {
         off = indp[ib];
	 endb = MIN(is+bsize, nsamples+off);
//#pragma omp for
         for (iy = is; iy < endb; iy++) {
            y[iy] += x[iy-off];
         }
      }
   } 
//   printf("-----------\n"); fflush(stdout);

   for (is = nsamples; is < nsamples+indn[nn-1]; is += bsize) {
      for (ib = 0; ib < nn; ib++) {
         off = indn[ib];
	 endb = MIN(is+bsize, nsamples+off);
//#pragma omp for
         for (iy = is; iy < endb; iy++) {
            y[iy] -= x[iy-off];
         }
      }
   }
//}  // end pragma parallel
#endif
//   printf("========================\n"); fflush(stdout);

}

// Version 4.
// As in version 3, but avoids looping twice over the indices and elements of y by performing and additional product as in v1
// Parallelism over y
// For every nonzero element of b all values of a n-size block of y can be updated in parallel 
// using the same sample
void processAudio_v4(btype *x, int nsamples, btype *b, int nb, int *ind, int nind, btype *y, int bsize) {
   int off, ith, iy, endb, end;
   btype elemb;

//#pragma omp parallel private(ith)
//{ 
//   ith = omp_get_thread_num();
   // Process first block of nb samples increasing the number of samples processed on each iteration
//   printf("start: %d, end: %d\n", 0, nb);
   for (int is = ind[0]; is <= nb; is += bsize ) {
      for (int ib = 0; ib < nind; ib++) {
         endb = MIN(is+bsize, nb);
         off = ind[ib];
         elemb = b[ib];
//#pragma omp for 
         for (int iy = MAX(is, ind[ib]); iy < endb; iy++) {
            y[iy] += x[iy-off]*elemb;
//  	    printf("[%d]: is: %d, ib: %d, y[%d] += x[%d]*b[%d]\n", ith, is, ib, iy, iy-off, ind[ib]); fflush(stdout);
         }
      } 
   } // end is
//   printf("-----------\n"); fflush(stdout);


   for (int iblock = nb; iblock < nsamples; iblock += nb) {

      end = MIN(iblock+nb, nsamples);
//      printf("start: %d, end: %d\n", iblock, MIN(iblock+nb, nsamples));
      for (int is = iblock; is <= end; is += bsize ) {
         for (int ib = 0 ; ib < nind ; ib++) {
	    endb = MIN(is+bsize, end);
            off = ind[ib];
            elemb = b[ib];
//#pragma omp for 
            for (int iy = is; iy < endb; iy++) {
               y[iy] += x[iy-off]*elemb;
//  	       printf("[%d]: is: %d, ib: %d, y[%d] += x[%d]*b[%d]\n", ith, is, ib, iy, iy-off, ind[ib]); fflush(stdout);
           }
         } 
     } // end is
//     printf("-----------\n"); fflush(stdout);

   } // end for iblock

#ifdef LAST_NB
   // Process the last block of samples until the pipeline is empty
//   printf("last start: %d, end: %d\n", nsamples, nsamples+nb-1);
   for (int is = nsamples; is <= nsamples+ind[nind-1]; is += bsize ) {
      for (int ib = 0 ; ib < nind ; ib++) {
         off = ind[ib];
         elemb = b[ib];
	 endb = MIN(is+bsize, nsamples+off);
//#pragma omp for 
         for (int iy = is; iy < endb; iy++) {
            y[iy] += x[iy-off]*elemb;
//  	    printf("[%d]: is: %d, ib: %d, y[%d] += x[%d]*b[%d]\n", ith, is, ib, iy, iy-off, ind[ib]); fflush(stdout);
         }
      } 
   } // end is
//}  // end pragma parallel
#endif
//   printf("========================\n"); fflush(stdout);

}

// Version 5
// Parallelism over the indices,  reduction on y. 
// The problem of the reduction clause when applied to arrays is that it cannot be applied to 
// very large arrays due to memory issues: every thread has to allocate a copy of the reduced sections of the array
void processAudio_v5(btype *x, int nsamples, int nb, int *indp, int np, int *indn, int nn, btype *y) {
   int ix, is, end;

#pragma omp parallel private(ix, is)
{

   // Process first block of nb samples increasing the number of samples processed on each iteration
   printf("start: %d, end: %d\n", 0, nb);
#pragma omp for reduction(+:y[0:nb])
   for (int ib = 0; ib < np; ib++) {
      ix = indp[ib];
      for (is = indp[ib]; is < nb; is++) {
         y[is] += x[ix++];
//	 printf("[%d]: y[%d] += %d --> %d\n", ith, is, x[is],  y[is]); fflush(stdout);
      }
   } 
   printf("-----------\n"); fflush(stdout);

#pragma omp for reduction(+:y[0:nb])
   for (int ib = 0; ib < nn; ib++) {
      ix = indn[ib];
      for (is = indn[ib]; is < nb; is++) {
         y[is] -= x[ix++];
//	 printf("y[%d] -= %d --> %d\n", is, *ptrx, y[is]); fflush(stdout);
      }
   } 
   printf("-----------\n"); fflush(stdout);

   for (int iblock = nb; iblock < nsamples; iblock += nb) {
      end = MIN(iblock+nb, nsamples);

      printf("start: %d, end: %d\n", iblock, end);
//#pragma omp for reduction(+:y[iblock:iblock+nb])
#pragma omp for reduction(+:y[iblock:nb])
      for (int ib = np-1 ; ib >=0 ; ib--) {
	 ix = iblock-indp[ib];
         for (is = iblock; is < end; is++) {
            y[is] += x[ix++];
        }
      } 
      printf("-----------\n"); fflush(stdout);

//#pragma omp for reduction(+:y[iblock:iblock+nb])
#pragma omp for reduction(+:y[iblock:nb])
      for (int ib = nn-1 ; ib >=0 ; ib--) {
         ix = iblock - indn[ib];
         for (is = iblock; is < end; is++) {
            y[is] -= x[ix++];
         }
      }
      printf("-----------\n"); fflush(stdout);

   } // end for iblock


#ifdef LAST_NB
   // Process the last block of samples until the pipeline is empty
   printf("last start: %d, end: %d\n", nsamples, nsamples+nb-1);
//#pragma omp for reduction(+:y[nsamples:nsamples+nb-1])
#pragma omp for reduction(+:y[nsamples:nb-1])
   for (int ib = np-1 ; ib >= 0 ; ib--) {
      ix = nsamples - indp[ib];
      for (is = nsamples; ix < nsamples; is++) {
         y[is] += x[ix++];
      }
   } 
   printf("-----------\n"); fflush(stdout);

//#pragma omp for reduction(+:y[nsamples:nsamples+nb-1])
#pragma omp for reduction(+:y[nsamples:nb-1])
   for (int ib = nn-1 ; ib >= 0 ; ib--) {
      ix = nsamples - indn[ib];
      for (is = nsamples; ix < nsamples; is++) {
         y[is] -= x[ix++];
      }
   }
#endif
}  // end pragma parallel

   printf("========================\n"); fflush(stdout);

}

int main(int argc, char *argv[]) {

  const int MAXSAMPLES = 5000000;
  const int MAXINDICES = 100000;
  int nb = 88000; // size of vector b
//  int nb = 100034; // size of vector b
  int nind; // number of non_zero elements in b
  int nsamples, nres;

  // Por ahora asumo que los datos de entrada son numeros enteros
  btype **x;
  btype **y;        // resultado a calcular
  btype *b = (btype *) malloc(MAXINDICES*sizeof(btype));   // guardamos solo los elementos no nulos (+1,-1) de b
//  btype *bdense;
#ifdef GOLDEN_TEST
  btype *golden;   // resultado en matlab
#endif
  int *ind; // indices a los elementos no nulos de b
  int *indp; // indices a los elementos no nulos de b
  int *indn; // indices a los elementos no nulos de b

  int i, ch;

  int ntimes = 3;

  int np, nn;
  int ct;
  int nth, nchannels, nblocks, bsize;
  char *type;

  if (argc < 6) {
     printf("FORMAT: audio_omp <nth> <nch> <nblocks> <bsize> <nb> <type>\n");
     return 1;
  } else {
     nth = atoi(argv[1]);
     nchannels = atoi(argv[2]);
     nblocks = atoi(argv[3]);
     bsize = atoi(argv[4]);
     nb = atoi(argv[5]);
     type = argv[6];
  }
//  printf("Number of threads: %d\n", nth); fflush(stdout);


  // Build files names from paramenters
  char xname[40];
  char bname[40]; // = "../data/b_88000_100.txt";
  char goldename[40] = "../data/goldenOnes.txt";
  char indname[40]; // = "../data/ind_88000_100.txt";
  char datadir[20] = "../data/";
  sprintf(xname, "%sx_%s.txt", datadir, type);
//  printf("xname: %s\n", xname);
  sprintf(bname, "%sb_%d.txt", datadir, nb);
//  printf("bname: %s\n", bname);
  sprintf(goldename, "%sgolden_%d_%s.txt", datadir, nb, type);
//  printf("goldename: %s\n", goldename);
  sprintf(indname, "%sind_%d.txt", datadir, nb);
//  printf("indname: %s\n", indname);



  //----------------Datos de entrada: un wav con la señal de piano 
/* 
//  int *x;
  int channel_number = 1;
  unsigned char **filewav = (unsigned char **) malloc (channel_number*sizeof(unsigned char *));
  nsamples = OpenWavConvert32(&filewav[0],"../data/001_piano.wav");
//  x = (int *)filewav[0];
  printf("Read %d audio samples \n", nsamples);
//  PrintVectorInt(x, nx, "x");
  golden = (int *) malloc(nsamples*sizeof(int));

  nsamples = nsamples + nb; // Add one zero block to finish processing the last block of samples
  printf("resized nsamples: %d\n", nsamples); fflush(stdout);
  // increases the size of the buffer to complete the last block of size nb and fills the last samples with zeros
  x = (int *) calloc( nsamples, sizeof(int) );
  memcpy(x, filewav[0], nsamples*sizeof(int));
//  PrintVectorInt(x, nsamples, "x");
 
  y = (int *) malloc(nsamples*sizeof(int));
*/

//  nsamples = 2*nb;
//  for (int i = 0; i < nsamples; i++) x[i] = 1;

/*
#ifdef INT32
  nsamples = ReadVector("../data/x_i.txt", x);
#else
   #ifdef FLOAT
      nsamples = ReadVector("../data/x_f.txt", x);
   #else
      nsamples = ReadVector("../data/x_i.txt", x);
   #endif
#endif
  printf("Read %d samples.\n", nsamples); fflush(stdout);
*/

//  PrintVectorInt(x, nsamples, "x");


  nind = ReadVector(bname, b, MAXINDICES);
//  printf("Read %d values of b\n", nind); fflush(stdout);

  ind = (int *) malloc(nind*sizeof(int));
  indp = (int *) malloc(nind*sizeof(int));
  indn = (int *) malloc(nind*sizeof(int));

  ct = ReadVectorIndices(indname, ind);
//  printf("Read %d indices\n", ct); fflush(stdout);

  if (nblocks == -1) {  // Read samples from file
     nsamples = countLines(xname);
//     printf("Number of samples on file: %d\n", nsamples); fflush(stdout);
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
//      x[ch] = (btype *) malloc(nsamples*sizeof(btype) );

  if (nblocks == -1) { // Read samples from file
     nsamples = ReadVector(xname, x[0], nsamples);
//     printf("Read %d samples.\n", nsamples); fflush(stdout);
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

/*
#ifdef INT32
  ct = ReadVector("../data/golden_100034_i.txt", golden);
#else
   #ifdef FLOAT
       ct = ReadVector("../data/golden_100034_f.txt", golden);
   #else
      ct = ReadVector("../data/golden_100034_s.txt", golden);
   #endif
#endif
  printf("Read %d values of golden result.\n", ct); fflush(stdout);
*/

#endif

  y = (btype **) malloc(nchannels * sizeof(btype *));
  for (int i = 0; i < nchannels; i++)
      y[i] = (btype *) calloc(nres, sizeof(btype) );

  /*
  nsamples = 16;
  nb = 8;
  nind = 3;
  ind[0]=2;
  ind[1]=5;
  ind[2]=8;
  b[0]=-1;
  b[1]=1;
  b[2]=1;
  */

  // Separate indexes of elements +1 and -1 of vector b
  np = nn = 0;
  for (int i = 0; i < nind; i++) { 
    ind[i]; // To convert matlab (from 1) to C (from 0) indexing
    if (b[i] > 0) { 
       indp[np++] = ind[i];
//       printf("b[%d] = %d, indp[%d] = %d\n", i, b[i], np-1, indp[np-1]);
    }
    else { 
       indn[nn++] = ind[i];
//       printf("b[%d] = %d, indn[%d] = %d\n", i, b[i], nn-1, indn[nn-1]);
    }
  }
//  printf("Number of +1: %d; Number of -1: %d\n", np, nn);

  double start, end;
  omp_set_num_threads(nth);

//  printf("#th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d\n",
//              nth, nchannels, nblocks, bsize, nb, nsamples );

  start = omp_get_wtime();
 
  for (i = 0; i < ntimes; i++)
#pragma omp parallel for num_threads(nth)
     for (ch = 0; ch < nchannels; ch++)
        processAudio_v0(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);

  end = omp_get_wtime();
  printf("v0 #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
              nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );
#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif

  printf("********************\n");

  for (ch = 0; ch < nchannels; ch++)
     memset(y[ch], 0, nres*sizeof(btype));

  start = omp_get_wtime();

  for (i = 0; i < ntimes; i++)
#pragma omp parallel for num_threads(nth)
     for (ch = 0; ch < nchannels; ch++)
        processAudio_v1(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);

  end = omp_get_wtime();
  printf("v1 #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
              nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );
#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif

  printf("********************\n");

  for (ch = 0; ch < nchannels; ch++)
     memset(y[ch], 0, nres*sizeof(btype));

  start = omp_get_wtime();

  for (i = 0; i < ntimes; i++)
#pragma omp parallel for num_threads(nth)
     for (ch = 0; ch < nchannels; ch++)
        processAudio_v2(x[ch], nsamples, nb, indp, np, indn, nn, y[ch], bsize);

  end = omp_get_wtime();
  printf("v2 #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
              nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );
#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif

  printf("********************\n");

  for (ch = 0; ch < nchannels; ch++)
     memset(y[ch], 0, nres*sizeof(btype));

  start = omp_get_wtime();

  for (i = 0; i < ntimes; i++)
#pragma omp parallel for num_threads(nth)
     for (ch = 0; ch < nchannels; ch++)
        processAudio_v3(x[ch], nsamples, nb, indp, np, indn, nn, y[ch], bsize);

  end = omp_get_wtime();
  printf("v3 #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
              nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );
#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif

  printf("********************\n");

  for (ch = 0; ch < nchannels; ch++)
     memset(y[ch], 0, nres*sizeof(btype));

  start = omp_get_wtime();

  for (i = 0; i < ntimes; i++)
#pragma omp parallel for num_threads(nth)
     for (ch = 0; ch < nchannels; ch++)
        processAudio_v4(x[ch], nsamples, b, nb, ind, nind, y[ch], bsize);

  end = omp_get_wtime();
  printf("v4 #th: %d ; #ch: %d ; #blocks: %d ; bsize: %d ; nb: %d; #samples: %d ; Tot time: %8.5f s. ; Time/block: %6.2f ms.\n",
              nth, nchannels, nblocks, bsize, nb, nsamples, (end-start)/ntimes, 1000*(end-start)/(nblocks*ntimes) );
#ifdef GOLDEN_TEST
  for (ch = 0; ch < nchannels; ch++)
     printf("compareGolden: %d\n", compareGolden(y[ch], golden, nres) );
#endif


  // Free memory
//  free(filewav[0]);
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
