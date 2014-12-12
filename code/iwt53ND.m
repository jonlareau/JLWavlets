function x = iwt53ND(x,J,D)
nd = 1:ndims(x);
if nargin < 2
    J = 3;
end
if nargin < 3
    D = ndims(x):-1:1;
end
for i = D
    n  = nd;
    n(i) = 1;
    n(1) = i;
    x = permute(iwt53(permute(x,n),J),n);
end
%{
/**
 *  dwt53.c - Fast discrete biorthogonal CDF 5/3 wavelet forward and inverse transform (lifting implementation)
 *  
 *  This code is provided "as is" and is given for educational purposes.
 *  2007 - Gregoire Pau - gregoire.pau@ebi.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *tempbank=0;

/**
 *  fwt53 - Forward biorthogonal 5/3 wavelet transform (lifting implementation)
 *
 *  x is an input signal, which will be replaced by its output transform.
 *  n is the length of the signal, and must be a power of 2.
 *
 *  The first half part of the output signal contains the approximation coefficients.
 *  The second half part contains the detail coefficients (aka. the wavelets coefficients).
 *
 *  See also iwt53.
 */
void fwt53(double* x,int n) {
  double a;
  int i;

  // Predict 1
  a=-0.5;
  for (i=1;i<n-2;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  } 
  x[n-1]+=2*a*x[n-2];

  // Update 1
  a=0.25;
  for (i=2;i<n;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[0]+=2*a*x[1];

  // Scale
  a=sqrt(2.0);
  for (i=0;i<n;i++) {
    if (i%2) x[i]*=a;
    else x[i]/=a;
  }

  // Pack
  if (tempbank==0) tempbank=(double *)malloc(n*sizeof(double));
  for (i=0;i<n;i++) {
    if (i%2==0) tempbank[i/2]=x[i];
    else tempbank[n/2+i/2]=x[i];
  }
  for (i=0;i<n;i++) x[i]=tempbank[i];
}

/**
 *  iwt53 - Inverse biorthogonal 5/3 wavelet transform
 *
 *  This is the inverse of fwt53 so that iwt53(fwt53(x,n),n)=x for every signal x of length n.
 *
 *  See also fwt53.
 */
void iwt53(double* x,int n) {
  double a;
  int i;

  // Unpack
  if (tempbank==0) tempbank=(double *)malloc(n*sizeof(double));
  for (i=0;i<n/2;i++) {
    tempbank[i*2]=x[i];
    tempbank[i*2+1]=x[i+n/2];
  }
  for (i=0;i<n;i++) x[i]=tempbank[i];

  // Undo scale
  a=1/sqrt(2.0);
  for (i=0;i<n;i++) {
    if (i%2) x[i]*=a;
    else x[i]/=a;
  }

  // Undo update 1
  a=-0.25;
  for (i=2;i<n;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[0]+=2*a*x[1];

  // Undo predict 1
  a=0.5;
  for (i=1;i<n-2;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[n-1]+=2*a*x[n-2];
}

int main() {
  double x[32];
  int i;

  // Makes a fancy cubic signal
  for (i=0;i<32;i++) x[i]=i; //5+i+0.4*i*i-0.02*i*i*i;
  
  // Prints original sigal x
  printf("Original signal:\n");
  for (i=0;i<32;i++) printf("x[%d]=%f\n",i,x[i]);
  printf("\n");

  // Do the forward 5/3 transform
  fwt53(x,32);
  
  // Prints the wavelet coefficients
  printf("Wavelets coefficients:\n");
  for (i=0;i<32;i++) printf("wc[%d]=%f\n",i,x[i]);
  printf("\n");

  // Do the inverse 5/3 transform
  iwt53(x,32); 

  // Prints the reconstructed signal 
  printf("Reconstructed signal:\n");
  for (i=0;i<32;i++) printf("xx[%d]=%f\n",i,x[i]);
}



%}