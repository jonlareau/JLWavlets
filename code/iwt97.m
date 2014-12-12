function x = iwt97(x,J)

if nargin < 2
    szx = size(x);
    J = min(floor(log2(szx(1))));
end

n = size(x,1);
sz = size(x);

%Reshape...
x = x(1:n,:);

for n = floor(n./(2.^(J-1:-1:0)))
    %unpack
    x([1:2:n,2:2:n]',:) = x(1:n,:);
    
    %// Undo scale
    a=1.149604398;
    i = (1:n)';
    mi = 0~=mod(i-1,2);
    x(mi,:) = x(mi,:)*a;
    x(~mi,:) = x(~mi,:)/a;
    
    %// Undo update 2
    a=-0.4435068522;
    i = (3:2:n-1)';
    x(i,:)=x(i,:)+a*(x(i-1,:)+x(i+1,:));
    x(1,:)=x(1,:)+2*a*x(2,:);
    
    %// Undo predict 2
    a=-0.8829110762;
    i = (2:2:n-1)';
    x(i,:)=x(i,:)+a*(x(i-1,:)+x(i+1,:));
    x(n,:)=x(n,:)+2*a*x(n-1,:);
    
    %// Undo update 1
    a=0.05298011854;
    i = (3:2:n-1)';
    x(i,:)=x(i,:)+a*(x(i-1,:)+x(i+1,:));
    x(1,:)=x(1,:)+2*a*x(2,:);
    
    %// Undo predict 1
    a=1.586134342;
    i = (2:2:n-1)';
    x(i,:)=x(i,:)+a*(x(i-1,:)+x(i+1,:));
    x(n,:)=x(n,:)+2*a*x(n-1,:);
    
end

x = reshape(x,sz);



%{
/**
 *  iwt97 - Inverse biorthogonal 9/7 wavelet transform
 *
 *  This is the inverse of fwt97 so that iwt97(fwt97(x,n),n)=x for every signal x of length n.
 *
 *  See also fwt97.
 */
void iwt97(double* x,int n) {
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
  a=1.149604398;
  for (i=0;i<n;i++) {
    if (i%2) x[i]*=a;
    else x[i]/=a;
  }

  // Undo update 2
  a=-0.4435068522;
  for (i=2;i<n;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[0]+=2*a*x[1];

  // Undo predict 2
  a=-0.8829110762;
  for (i=1;i<n-2;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[n-1]+=2*a*x[n-2];

  // Undo update 1
  a=0.05298011854;
  for (i=2;i<n;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[0]+=2*a*x[1];

  // Undo predict 1
  a=1.586134342;
  for (i=1;i<n-2;i+=2) {
    x[i]+=a*(x[i-1]+x[i+1]);
  }
  x[n-1]+=2*a*x[n-2];
}

int main() {
  double x[32];
  int i;

  // Makes a fancy cubic signal
  for (i=0;i<32;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  
  // Prints original sigal x
  printf("Original signal:\n");
  for (i=0;i<32;i++) printf("x[%d]=%f\n",i,x[i]);
  printf("\n");

  // Do the forward 9/7 transform
  fwt97(x,32);
  
  // Prints the wavelet coefficients
  printf("Wavelets coefficients:\n");
  for (i=0;i<32;i++) printf("wc[%d]=%f\n",i,x[i]);
  printf("\n");

  // Do the inverse 9/7 transform
  iwt97(x,32);

  // Prints the reconstructed signal
  printf("Reconstructed signal:\n");
  for (i=0;i<32;i++) printf("xx[%d]=%f\n",i,x[i]);
}

%}