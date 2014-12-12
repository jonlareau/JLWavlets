function x = iwt53(x,J)
%Adapted from http://www.embl.de/~gpau/misc/dwt53.c
if nargin < 2
    J = 8;
end
if nargin ==0
    x0 = double(imread('football.jpg'));
    %x0 = x0-mean(x0(:));
    %x0 = x0/max(abs(x0(:)));
    x0 = (x0-128)/256;
    figure(1)
    imagesc((x0-min(x0(:)))/(max(x0(:))-min(x0(:))));
    title('input image');
    figure(2)
    x = fwt53ND(x0,J,1:2);
    x = permute(x,[2,1,3]);
    imagesc(x(:,:,1)); colorbar;
    imagesc((x-min(x(:)))/(max(x(:))-min(x(:))));
    title('test image');
    drawnow;
    figure(3)
end

n = size(x,1);
sz = size(x);

%Reshape...
x = x(1:n,:);

for n = n./(2.^(J-1:-1:0))
    %unpack
    x([1:2:n,2:2:n]',:) = x(1:n,:);
    
    %// Undo scale
    a=1/sqrt(2.0);
    for (i=1:n)
        if mod(i-1,2)~=0
            x(i,:) = (x(i,:)/a);
        else
            x(i,:)=(x(i,:)*a);
        end
    end

    %// Undo update 2
    a=-0.25;
    for (i=3:2:n)
        x(i,:)=x(i,:)+(a*(x(i-1,:)+x(i+1,:)));
    end
    x(1,:)=x(1,:)+(2*a*x(2,:));
  
    %// Undo predict 2
    a=.5;
    for (i=2:2:n-1)
        x(i,:)=x(i,:)+(a*(x(i-1,:)+x(i+1,:)));
        %x[i]+=a*(x[i-1]+x[i+1]);
    end
    x(n,:)=x(n,:)+(2*a*x(n-1,:));
end

x = reshape(x,sz);

if nargin==0
    x1 = iwt53(permute(x,[2,1,3]),J);
    x1 = x1(1:size(x0,1),1:size(x0,2),1:size(x0,3));
    imagesc((x1-min(x1(:)))/(max(x1(:))-min(x1(:))));
    title('resultant image');
    
    mean(abs(x0(:)-x1(:)))
    %imagesc(abs(x1-x0)/max(abs(x1(:)-x0(:))));
    x = mean(abs(x0(:)-x1(:)));
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