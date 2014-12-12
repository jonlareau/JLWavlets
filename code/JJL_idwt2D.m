function y = JJL_idwt2D(w, J, sf)

% Inverse 2-D Discrete Wavelet Transform
%
% USAGE:
%   y = idwt(w, J, sf)
% INPUT:
%   w - wavelet coefficients
%   J  - number of stages
%   sf - synthesis filters
% OUTPUT:
%   y - output array
% See dwt2D
%
%Adapted by Jonathan Lareau from:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
    [af, sf] = farras;
end

s = size(w);
y = w;
for k = J:-1:1
    if iscell(sf)
        [sr, sc] = size(sf);
        if (sr <= 1) || (sc <= 1)
            sfk1 = sf{min(k,length(sf))};
            sfk2 = sf{min(k,length(sf))};
        else
            sfk1 = sf{min(k,sr),1};
            sfk2 = sf{min(k,sr),2};
        end
    else
        sfk1 = sf; sfk2 = sf;
    end
    j = s ./ (2^(k-1));
    [y(1:j(1),1:j(2))] = JJL_sfb2D(y(1:j(1),1:j(2)), [], sfk1, sfk2);
end

