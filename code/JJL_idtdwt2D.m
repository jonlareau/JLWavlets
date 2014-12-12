function [y] = JJL_idtdwt2D(wh,wg,J,Fsf,sf)
% Inverse 2-D Dual-Tree Discrete Wavelet Transform
% 
% USAGE:
%   y = idtdwt2D(w, J, hsf, gsf)
%   y = idtdwt2D(wh,wg,J,hsf,gsf);
% INPUT:
%   w{1} = wh, the output from the top dwt2D filter bank 
%   w{2} = wg, the output from the bottom dwt2D filter bank
%       w can also be a matrix of complex numbers, in which case wg =
%       imag(w), and wh = real(w);
%   J - number of stages
%   Fsf - First stage synthesis filters 
%   sf -  Subsequent stage synthesis filters 
% OUPUT:
%   y - output array
% See idualtree2D
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 4
    sf = Fsf;
    Fsf = J;
    J = wg;
    if iscell(wh)
        wg = wh{2};
        wh = wh{1};
    else
        wg = imag(wh);
        wh = real(wh);
    end
elseif nargin == 2
    J = wg;
    if iscell(wh)
        wg = wh{2};
        wh = wh{1};
    else
        wg = imag(wh);
        wh = real(wh);
    end
end
    
if nargin < 4
    [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J);
else
    [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J,[],[],Fsf,sf);
end

s0  = size(wh)/(2^(J));

%sum and difference terms...
w1 = (wh+wg)/sqrt(2);
w2 = (wh-wg)/sqrt(2);
w1(1:s0(1),1:s0(2)) = wh(1:s0(1),1:s0(2));
w2(1:s0(1),1:s0(2)) = wg(1:s0(1),1:s0(2));

% Tree 1
y1 = JJL_idwt2D(w1,J,hsf);

% Tree 2
y2 = JJL_idwt2D(w2,J,gsf);

% normalization
y = (y1 + y2)/sqrt(2);
%y = y/sqrt(2);
