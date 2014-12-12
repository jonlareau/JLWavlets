function out = JJL_dwt2D(x, J, af)

% discrete 2-D wavelet transform
%
% USAGE:
%   w = dwt2D(x, stages, af)
% INPUT:
%   x - N by M matrix
%       1) M, N are both even
%       2) min(M,N) >= 2^(J-1)*length(af)
%   J - number of stages
%   af - analysis filters
% OUPUT:
%   w - Matrix of wavelet coefficients
% EXAMPLE:
%   [af, sf] = farras;
%   x = rand(128,64);
%   J = 3;
%   w = dwt2D(x,J,af);
%   y = idwt2D(w,J,sf);
%   err = x - y; 
%   max(max(abs(err)))
%
%Adapted by Jonathan Lareau from:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
    [af, sf] = farras;
end

s = size(x);
out = x;
for k = 1:J
    if iscell(af)
        [sr, sc] = size(af);
        if (sr <= 1) || (sc <= 1)
            afk1 = af{min(k,length(af))};
            afk2 = af{min(k,length(af))};
        else
            afk1 = af{min(k,sr),1};
            afk2 = af{min(k,sr),2};
        end
    else
        afk1 = af; afk2 = af;
    end
        [out(1:s(1),1:s(2))] = JJL_afb2D(out(1:s(1),1:s(2)), afk1, afk2);
    s = s/2;
    if any(mod(s,2)~=0)
        'stoping recursion'
        break;
    end
end

