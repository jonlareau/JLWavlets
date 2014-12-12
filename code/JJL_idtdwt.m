function [y] = JJL_idtdwt(wh,wg,J,Fsf,sf)
% Inverse Dual-Tree Discrete Wavelet Transform
% 
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 4
    sf = Fsf;
    Fsf = J;
    J = wg;
    wg = wh{2};
    wh = wh{1};
elseif nargin == 2
    J = wg;
    wg = wh{2};
    wh = wh{1};
end
    
if nargin < 4
    [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J);
else
    [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J,[],[],Fsf,sf);
end
    

% Tree 1
y1 = JJL_idwt2D(wh,J,hsf);

% Tree 2
y2 = JJL_idwt2D(wg,J,gsf);

% normalization
y = (y1 + y2)/sqrt(2);
%y = y/sqrt(2);
