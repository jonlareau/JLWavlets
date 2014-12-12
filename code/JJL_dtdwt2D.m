function varargout = JJL_dtdwt2D(x,J,Faf,af,mode,ALTERNATE)
%Real 2D Dual-Tree Discrete Wavelet Transform
%
%The Real 2D DT-DWT is computed by using two seperable 2D DWT's in
%parrallel and then computing their sum and difference terms.
%
%USAGE:
%   [ws,wd] = JJL_dtdwt2D(x,J)
%   [ws,wd] = JJL_dtdwt2D(x,J,Faf,af)
%   [ws,wd] = JJL_dtdwt2D(x,J,Faf,af,OutMode)
%
%   ws - Sum Terms
%   wd - Difference Terms
%
%Note: this 
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 0
    x = rand(64,64);
    J = 4;
end

if nargin < 5 || isempty(mode)
    mode = 'reg';
end

if nargin >= 4
    [haf,gaf] = JJL_makeDTFilters(J,Faf,af);
elseif nargin < 4
    [haf,gaf] = JJL_makeDTFilters(J);
end

% normalization
x = x/sqrt(2);

s0  = size(x)/(2^(J));

w1 = JJL_dwt2D(x,J,haf);
w2 = JJL_dwt2D(x,J,gaf);

%sum and difference terms...
wh = (w1+w2)/sqrt(2);
wg = (w1-w2)/sqrt(2);
wh(1:s0(1),1:s0(2)) = w1(1:s0(1),1:s0(2));
wg(1:s0(1),1:s0(2)) = w2(1:s0(1),1:s0(2));

if nargout == 1
    if strcmpi(mode,'combine')
        varargout = wh+sqrt(-1)*wg;
    else
        varargout{1} = {wh,wg};
    end
else
    varargout{1} = wh;
    varargout{2} = wg;
end