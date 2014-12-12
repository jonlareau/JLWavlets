function y = JJL_sfb2D(lo, hi, sf1, sf2)

% 2D Synthesis Filter Bank
%
% USAGE:
%   y = sfb2D(lo, hi, sf1, sf2); (For cell wavelet coeffs)
%   y = sfb2D(w, [], sf1, sf2); (For matrix wavelet coeffs)
% INPUT:
%   lo, hi - lowpass, highpass subbands
%   sf1 - synthesis filters for the columns
%   sf2 - synthesis filters for the rows
% OUTPUT:
%   y - output array
% See afb2D
%
%Adapted by Jonathan Lareau From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/


if nargin < 4
    sf2 = sf1;
end

if isempty(hi)
    hi{1} = lo((end/2+1):end,1:end/2);
    hi{2} = lo(1:end/2,(end/2+1):end);
    hi{3} = lo((end/2+1):end,(end/2+1):end);
    lo = lo(1:end/2,1:end/2);
end

% filter along rows
lo = JJL_sfb2D_A(lo,    hi{1}, sf2, 2);
hi = JJL_sfb2D_A(hi{2}, hi{3}, sf2, 2);

% filter along columns
y = JJL_sfb2D_A(lo, hi, sf1, 1);


