function [o] = JJL_afb(x, af)

% Analysis filter bank
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

N = length(x);
L = length(af)/2;
x = cshift(x,-L);

% lowpass filter
lo = upfirdn(x, af(:,1), 1, 2);
lo(1:L) = lo(N/2+[1:L]) + lo(1:L);
lo = lo(1:N/2);

% highpass filter
hi = upfirdn(x, af(:,2), 1, 2);
hi(1:L) = hi(N/2+[1:L]) + hi(1:L);
hi = hi(1:N/2);

o = zeros(size(x));
o(1:end/2) = lo;
o((end/2+1):end) = hi;
