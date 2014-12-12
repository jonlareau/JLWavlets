function y = JJL_idwt3D(w,J,sf)

% Inverse 3-D Discrete Wavelet Transform
%
% USAGE:
%   y = idwt3D(w, J, sf)
% INPUT:
%   w - wavelet coefficient
%   J  - number of stages
%   sf - synthesis filters
% OUTPUT:
%   y - output array
% See: dwt3D
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

%
    s = size(w);
    y = w;
    for k = J:-1:1
    if iscell(sf)
        [sr, sc] = size(sf);
        if (sr <= 1) || (sc <= 1)
            sfk1 = sf{min(k,length(sf))};
            sfk2 = sf{min(k,length(sf))};
            sfk3 = sf{min(k,length(sf))};
        else
            sfk1 = sf{min(k,sr),1};
            sfk2 = sf{min(k,sr),2};
            sfk3 = sf{min(k,sr),3};
        end
    else
        sfk1 = sf; sfk2 = sf; sfk3 = sf;
    end
        j = s ./ (2^(k-1));
        [y(1:j(1),1:j(2),1:j(3))] = JJL_sfb3D(y(1:j(1),1:j(2),1:j(3)), [], sfk1, sfk2, sfk3);
    end

