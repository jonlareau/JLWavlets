function varargout = JJL_dtdwt(x,J,Faf,af,mode)
% 1D Dual-Tree Discrete Wavelet Transform
%
%USAGE:
%   [wh, wg] = JJL_dtdwt(x,J,Faf,af)
%   [wh, wg] = JJL_dtdwt(x,J)
%   [w] = JJL_dtdwt(...)
%
%   x - Input data vector
%   J - # scales to evaluate
%   Faf - First stage analysis filters
%   af - subsequent stage analysis filters
%
%   w = {wh,wg}
%   wh - wavelet coefficients from 1st filterbank (~Real Part)
%   wg - wavelet coefficients from 2nd filterbank (~Imag Part)
%
%Adapted from:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 0
    x = rand(64,1);
    J = 4;
end
if nargin <5
    mode = 'reg';
end

if nargin >= 4
    [haf,gaf] = JJL_makeDTFilters(J,Faf,af);
else
    [haf,gaf] = JJL_makeDTFilters(J);
end

% normalization
x = x/sqrt(2);
wh = JJL_dwt(x,J,haf);
wg = JJL_dwt(x,J,gaf);

if nargout == 1
    if strcmpi(mode,'combine')
        varargout{1} = wh+sqrt(-1)*wg;
    else
        varargout{1} = {wh,wg};
    end
else
    varargout{1} = wh;
    varargout{2} = wg;
end
