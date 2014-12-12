function out = JJL_afb2D(x, af1, af2)

% 2D Analysis Filter Bank
%
%Adapted by Jonathan Lareau From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
   af2 = af1;
end

% filter along columns
[L, H] = JJL_afb2D_A(x, af1, 1);

out = zeros(size(x));
% filter along rows
[out(1:end/2,1:end/2), out((end/2+1):end,1:end/2)] = JJL_afb2D_A(L, af2, 2);
[out(1:end/2,(end/2+1):end), out((end/2+1):end,(end/2+1):end)] = JJL_afb2D_A(H, af2, 2);

