function y = JJL_idwt(w, J, sf)
% Inverse Discrete 1-D Wavelet Transform
%
% USAGE:
%    y = idwt(w, J, sf)
% INPUT:
%    w - wavelet coefficients
%    J - number of stages
%    sf - synthesis filters
% OUTPUT:
%    y - output signal
% See also dwt
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
    s = size(w);
    y = w;
    for k = J:-1:1
        if iscell(sf)
            sfk = sf{min(k,length(sf))};
        else
            sfk = sf;
        end
        j = s ./ (2^(k-1));
        [y(1:j(1))] = JJL_sfb(y(1:j(1)), [], sfk);
    end
    

