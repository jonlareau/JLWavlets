function [varargout] = JJL_cplxdtdwt2D(x, J, Faf, af)
% Dual-Tree Complex 2D Discrete Wavelet Transform
%
%USAGE:
%   [W] = JJL_cplxdtdwt2D(x, J)
%   [W1,W2] = JJL_cplxdtdwt2D(x, J)
%   [Wr1,Wi1,Wr2,Wi2] = JJL_cplxdtdwt2D(x, J)
%
%   W - {Wr1 Wr2; Wi1 Wi2}
%   W1 - Wr1 + sqrt(-1)*Wi1
%   W2 - Wr2 + sqrt(-1)*Wi2
%   Wr1,Wi1,Wr2,Wi2 - the real and imaginary wavelet coefficients for each
%       direction
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
    [Faf, Fsf] = FSfarras;
    [af, sf] = dualfilt1;
end

% normalization
x = x/2;
AF = [];
for m = 1:2
    for n = 1:2
        AF{1,1} = Faf{m};
        AF{1,2} = Faf{n};
        AF{2,1} = af{m};
        AF{2,2} = af{n};
        w{m,n} = JJL_dwt2D(x,J,AF);
    end
end

s0 = size(x) / 2^(J);
[wr1 wi2] = pm(w{1,1},w{2,2});
[wr2 wi1] = pm(w{1,2},w{2,1});
wr1(1:s0(1),1:s0(2)) = w{1,1}(1:s0(1),1:s0(2));
wi1(1:s0(1),1:s0(2)) = w{2,1}(1:s0(1),1:s0(2));
wr2(1:s0(1),1:s0(2)) = w{1,2}(1:s0(1),1:s0(2));
wi2(1:s0(1),1:s0(2)) = w{2,2}(1:s0(1),1:s0(2));

%Format the output accordingly...
if nargout == 1
    varargout{1} = {wr1 wr2;wi1 wi2};
elseif nargout == 2
    varargout{1} = wr1 + j * wi1;
    varargout{2} = wr2 + j * wi2;
elseif nargout == 4
    varargout{1} = wr1;
    varargout{2} = wi1;
    varargout{3} = wr2;
    varargout{4} = wi2;
end
 
    

