function w = JJL_dwt(x, J, af)

% Discrete 1-D Wavelet Transform
%
% USAGE:
%    w = dwt(x, J, af)
% INPUT:
%    x - N-point vector, where
%            1) N is divisible by 2^J
%            2) N >= 2^(J-1)*length(af)
%    J - number of stages
%    af - analysis filters
%    af(:, 1) - lowpass filter (even length)
%    af(:, 2) - highpass filter (evenlength)

s = size(x);
w = x;
for k = 1:J
    [w(1:s(1))] = JJL_afb(w(1:s(1)), af);
    s = s/2;
end
