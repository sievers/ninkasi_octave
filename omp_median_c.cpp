#include <octave/oct.h>
#include <iostream>
#include <stdio.h>


/*--------------------------------------------------------------------------------*/

#define SSWAP(a,b) temp=(a);(a)=(b);(b)=temp;

/*!
 * Returns the kth smallest value in the input array arr[1...n].
 * Notice the twisted deal where the first valid array value is arr[1], not arr[0].  The caller has
 * to make sure not to screw that up.
 * The input array will be rearranged such that the return value lives in arr[k], and all values at
 * lower index are <=arr[k] while all values at higher index are >=arr[k].
 * If this routine looks kind of like it came from the Second Edition of a certain book of
 * programming recipes, section 8.5, then you are imagining things.  That would be unethical.
 * \param k   The index of the value that we want (counting into a sorted array).
 * \param n   Length of the data array.
 * \param arr The data array to select from.  WILL BE MODIFIED (partially sorted) by this operation.
 */

double sselect(unsigned long k, unsigned long n, double *arr)
{
  unsigned long i,ir,j,l,mid;
  double a,temp;

  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SSWAP(arr[l],arr[ir]);
      }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SSWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
        SSWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SSWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SSWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SSWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef SSWAP


/*---------------------------------------------------------------------------------------------------------*/


DEFUN_DLD (omp_median_c, args, nargout, "Get the median of an array along columns using openmp.\n")
{
  Matrix dat=args(0).matrix_value();
  dim_vector dm=dat.dims();
  unsigned long n=dm.elem(0);
  unsigned long m=dm.elem(1);

  Matrix meds(1,m);
  
#pragma omp parallel for shared(dat,m,n,meds) default(none)  
  for (long i=0;i<m;i++) {
    ColumnVector c=dat.column(i);
    double *d=c.fortran_vec();
    
    meds(i)=sselect((n+1)/2,n,d-1);
    //meds(i)= sselect(n/2,  n, d-1);


  }

  return octave_value(meds);

}

/*---------------------------------------------------------------------------------------------------------*/


DEFUN_DLD (omp_median_r, args, nargout, "Get the median of an array along rows using openmp.\n")
{
  Matrix dat=args(0).matrix_value();
  dim_vector dm=dat.dims();
  unsigned long n=dm.elem(0);
  unsigned long m=dm.elem(1);

  Matrix meds(n,1);
  
#pragma omp parallel for shared(dat,m,n,meds) default(none)  
  for (long i=0;i<n;i++) {
    RowVector r=dat.row(i);
    double *d=r.fortran_vec();
    meds(i)=sselect((m+1)/2,m,d-1);
  }
  
  return octave_value(meds);
  
}
