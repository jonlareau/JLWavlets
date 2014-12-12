function [y]= JJL_sfb3D(lo, hi, sf1, sf2, sf3)

% 3D Synthesis Filter Bank
%
% USAGE:
%   y = sfb3D(lo, hi, sf1, sf2, sf3);
% INPUT:
%   lo, hi - lowpass subbands
%   sfi - synthesis filters for dimension i
% OUPUT:
%   y - output array
% See afb3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 4
    sf2 = sf1;
    sf3 = sf1;
end

if isempty(hi)
    
    LLL = lo(1:end/2,1:end/2,1:end/2);
    LLH = lo((end/2+1):end,1:end/2,1:end/2);
    LHL = lo(1:end/2,(end/2+1):end,1:end/2);
    LHH = lo((end/2+1):end,(end/2+1):end,1:end/2);

    HLL = lo(1:end/2,1:end/2,(end/2+1):end);
    HLH = lo((end/2+1):end,1:end/2,(end/2+1):end);
    HHL = lo(1:end/2,(end/2+1):end,(end/2+1):end);
    HHH = lo((end/2+1):end,(end/2+1):end,(end/2+1):end);
else
    LLL = lo;
    LLH = hi{1};
    LHL = hi{2};
    LHH = hi{3};
    HLL = hi{4};
    HLH = hi{5};
    HHL = hi{6};
    HHH = hi{7};
end
% filter along dimension 3
LL = JJL_sfb3D_A(LLL, LLH, sf3, 3);
LH = JJL_sfb3D_A(LHL, LHH, sf3, 3);
HL = JJL_sfb3D_A(HLL, HLH, sf3, 3);
HH = JJL_sfb3D_A(HHL, HHH, sf3, 3);

% filter along dimension 3
L = JJL_sfb3D_A(LL, LH, sf2, 2);
H = JJL_sfb3D_A(HL, HH, sf2, 2);

% filter along dimension 1
y = JJL_sfb3D_A(L, H, sf1, 1);


% LOCAL FUNCTION

function y = JJL_sfb3D_A(lo, hi, sf, d)

% 3D Synthesis Filter Bank
% (along single dimension only)
%
% y = sfb3D_A(lo, hi, sf, d);
% sf - synthesis filters
% d  - dimension of filtering
% see afb2D_A

lpf = sf(:, 1);     % lowpass filter
hpf = sf(:, 2);     % highpass filter

% permute dimensions of lo and hi so that dimension d is first.
p = mod(d-1+[0:2], 3) + 1;
lo = permute(lo, p);
hi = permute(hi, p);

[N1, N2, N3] = size(lo);
N = 2*N1;
L = length(sf);
y = zeros(N+L-2, N2, N3);

for k = 1:N3
    y(:, :, k) = upfirdn(lo(:, :, k), lpf, 2, 1) + upfirdn(hi(:, :, k), hpf, 2, 1);
end
y(1:L-2, :, :) = y(1:L-2, :, :) + y(N+[1:L-2], :, :);
y = y(1:N, :, :);
y = cshift3D(y, 1-L/2, 1);

% permute dimensions of y (inverse permutation)
y = ipermute(y, p);


