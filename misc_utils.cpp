#include <octave/oct.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif




//Random utilities that may be useful






int get_int_value(octave_value val)
{
  NDArray myptr=val.array_value();
  int myval=(int)myptr(0,0);
  return myval;

}



/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
int *ivector(int n)
{
  int *ivec=(int *)malloc(sizeof(int)*n);
  if (ivec)
    return ivec;
  else {
    printf("Malloc error in ivector of %d\n",n);
    return NULL;
  }
  return NULL;
}
/*--------------------------------------------------------------------------------*/
int **imatrix(int n, int m)
{
  int *ivec=(int *)malloc(sizeof(int)*n*m);
  if (!ivec) {
    printf("Malloc failure in imatrix.\n");
    return NULL;
  }
  int **imat=(int **)malloc(sizeof(int *)*n);
  if (!imat) {
    printf("Malloc failure in imatrix.\n");
    return NULL;
  }

  for (int i=0;i<n;i++)
    imat[i]=ivec+i*m;
  return imat;
}

/*--------------------------------------------------------------------------------*/
double *dvector(int n)
{
  double *dvec=(double *)malloc(sizeof(double)*n);
  if (dvec)
    return dvec;
  else {
    printf("Malloc error in dvector of %d\n",n);
    return NULL;
  }
  return NULL;
}
/*--------------------------------------------------------------------------------*/


int isprime(int n)
//slow prime tester, but not a big deal                                                                                                                                                                    
{
  int myrt=sqrt(n);
  for (int i=2;i<=myrt;i++)
    if (n%i==0)
      return 0;
  return 1;
}

/*--------------------------------------------------------------------------------*/

int nprimes(int n)
//find number of primes less than or equal to n                                                                                                                                                            
{
  int nn=0;
  for (int i=2;i<=n;i++)
    if (isprime(i))
      nn++;
  return nn;
}

/*--------------------------------------------------------------------------------*/
int fill_in_logs(double lognmax, double *plogs, int nprime, int **facs, int icur)
{

  if (nprime==1) {
    int tmp=0;
    for (int i=0;i<=lognmax/plogs[0];i++) {
      facs[0][icur+tmp]=i;
      tmp++;
    }
    return tmp;
  }
  else {
    int tot=0; //accumulate how many factors we have at this depth                                                                                                                                         
    for (int i=0;i<=lognmax/plogs[nprime-1];i++) {
      int itmp=fill_in_logs(lognmax-plogs[nprime-1]*i,plogs,nprime-1,facs,icur+tot);

      for (int j=0;j<itmp;j++)
        facs[nprime-1][j+tot+icur]=i;

      tot+=itmp;
      //abort for testing                                                                                                                                                                                  
      //if (nprime==3)                                                                                                                                                                                     
      //return tot;                                                                                                                                                                                        
    }
    return tot;
  }
  assert(1==0);  //never get here                                                                                                                                                                          
}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (find_good_fft_lens_fast_c, args, nargout, "Find easily factorizable numbers up to some value.  args are (max_val,max_prime).\n")
{
  int nmax=1000;
  int maxprime=7;
  if (args.length()>0)
    nmax=get_int_value(args(0));

  if (args.length()>1)
    maxprime=get_int_value(args(1));


  int nprime=nprimes(maxprime);
  int *primes=ivector(nprime);
  int icur=0;
  for (int i=2;i<=maxprime;i++) {
    if (isprime(i)) {
      primes[icur]=i;
      icur++;
    }
  }


  double nlog=log(nmax);
  double *plogs=dvector(nprime);
  for (int i=0;i<nprime;i++)
    plogs[i]=log(primes[i]);

  //find gross upper bound on number of possible factors                                                                                                                                                   
  int maxfacs=1;
  for (int i=0;i<nprime;i++) {
    int tmp=1+nlog/plogs[i];
    maxfacs*=tmp;
  }
  for (int i=2;i<nprime;i++)
    maxfacs/=2;
  maxfacs+=100;  //pad by 100 to make sure we don't run off the end.                                                                                                                                        printf("maximum possible is %d\n",maxfacs); 


  int **myfacs=imatrix(nprime,maxfacs);
  for (int i=0;i<nprime;i++)
    for (int j=0;j<maxfacs;j++)
      myfacs[i][j]=-1;


  int nfacs=fill_in_logs(nlog,plogs,nprime,myfacs,0);

  int *allnums=ivector(nfacs);
  for (int i=0;i<nfacs;i++) {
    double tot=0;
    for (int j=0;j<nprime;j++)
      tot+=myfacs[j][i]*plogs[j];
    allnums[i]=(int)(exp(tot)+0.01);
  }


  ColumnVector ans(nfacs);
  for (int i=0;i<nfacs;i++)
    ans(i)=allnums[i];
  
  free(allnums);
  free(myfacs[0]);
  free(myfacs);
  free(primes);
  free(plogs);

  return octave_value(ans);
  
}

