#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <sys/time.h>
typedef struct timeval pca_time;

#define ACTSINGLE

#ifdef ACTSINGLE
typedef float actData;
#else
typedef double actData;
#endif


/*--------------------------------------------------------------------------------*/

void tick(pca_time *tt)
{
  gettimeofday(tt,NULL);
}
/*--------------------------------------------------------------------------------*/
void tock(pca_time *tt)
{
  pca_time tnow;
  gettimeofday(&tnow,NULL);
  double   dt=(tnow.tv_usec-tt->tv_usec)/1.0e6+(tnow.tv_sec-tt->tv_sec);  
  fprintf(stdout,"Tock registers %14.4e seconds.\n",dt);
  
}
/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/

actData **matrix(long n,long m)
{
  actData *data, **ptrvec;
  data=(actData *)malloc(sizeof(actData)*n*m);
  assert(data!=NULL);
  ptrvec=(actData **)malloc(sizeof(actData *)*n);
  assert(ptrvec!=NULL);
   for (int i=0;i<n;i++)
     ptrvec[i]=&data[i*m];
   
   return ptrvec;
}
/*--------------------------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"
{
#endif


void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc, int alen, int blen);

void sgemm_(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc, int alen, int blen);

#ifdef __cplusplus
}
#endif


/*--------------------------------------------------------------------------------*/
void actgemm(char transa, char transb, int m, int n, int k, actData alpha, actData *a, int lda, actData *b, int ldb, actData beta, actData *c, int ldc)
{
#ifdef ACTSINGLE
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1,1);
#else
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1,1);

#endif
}

/*--------------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
  int ndata=atoi(argv[1]);
  int ndet=atoi(argv[2]);
  int nvec=atoi(argv[3]);
  
  
  actData **mat=matrix(ndet,ndata);
  actData **cols=matrix(ndata,nvec);
  actData **rows=matrix(nvec,ndet);




  pca_time tt;
  tick(&tt);
#pragma omp parallel for shared(ndet,ndata,nvec,cols,rows,mat) default(none)  
  for (int i=0;i<ndet;i++)
    for (int j=0;j<ndata;j++)
      for (int k=0;k<nvec;k++) {
	mat[i][j]+=cols[j][k]*rows[k][i];
      }
  tock(&tt);

  if (nvec==1) {
    printf("Doing single vec.\n");
    tick(&tt);
#pragma omp parallel for shared(ndet,ndata,nvec,cols,rows,mat) default(none)  
    for (int i=0;i<ndet;i++)
      for (int j=0;j<ndata;j++)
	mat[i][j]+=cols[0][j]*rows[0][i];

    tock(&tt);
  }
  
  return 0;
  tick(&tt);
  actgemm('n','n',ndet,ndata,nvec,1.0,rows[0],ndet,cols[0],nvec,0.0,mat[0],ndata);
  tock(&tt);
  
}
