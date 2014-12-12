function [w1 w2 w3 w4] = JJL_dtdwt3D(x, J, Faf, af)
%Reak 3D Dual-Tree Discrete Wavelet Transform
%
%The Real 3D DT-DWT is computed by using four seperable 3D DWT's in
%parrallel and then computing their subbands.
%
%USAGE:
%   [u,v,q,r] = JJL_dtdwt3D(x,J,Faf,af)
%   [w] = JJL_dtdwt3D(x,J,Faf,af)
%
%   w = {u,v,q,r}
%   u = (w1 + w2 - w3 - w4)/2;
%   v = (w1 - w2 + w3 + w4)/2; 
%   q = (w1 + w2 - w3 + w4)/2;
%   r = (w1 + w2 + w3 - w4)/2;
%
%   w1,w2,w3,w4 - The four separate sets of wavelet coefficients...
%
%   NOTE: the low frequency information for each subband is not altered:
%       ie// u(LLL) = w1(LLL)
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin == 0
    %For debugging purposes only
    J = 3;
    x = rand(64,64,64);
    [Faf,Fsf] = FSfarras;
    [af,sf] = dualfilt1;
end

% normalization
x = x/2;

M = [
    1 1 1
    2 2 1
    2 1 2
    1 2 2
    ];

w = cell(1,4);
for i = 1:4
    %Set up each of the channel filters...
    f1 = M(i,1);
    f2 = M(i,2);
    f3 = M(i,3);
    AF = cell(2,3);

    %Set up the filter alternations for each scale...
    [AF{1,1},AF{1,2},AF{1,3}] = deal(Faf{f1}, Faf{f2}, Faf{f3}); %First stage
    [AF{2,1},AF{2,2},AF{2,3}] = deal(af{f1}, af{f2}, af{f3}); %Subsequent stages

    %Perform the wavelet transform for this channel
    w{i} = JJL_dwt3D(x,J,AF);
end

%Sum and difference terms...
[w1 w2 w3 w4] = pm4(w{1}, w{2}, w{3}, w{4});

%Keep the Low frequency data untouched...
s0  = size(x)/(2^(J));
w1(1:s0(1),1:s0(2),1:s0(3)) = w{1}(1:s0(1),1:s0(2),1:s0(3));
w2(1:s0(1),1:s0(2),1:s0(3)) = w{2}(1:s0(1),1:s0(2),1:s0(3));
w3(1:s0(1),1:s0(2),1:s0(3)) = w{3}(1:s0(1),1:s0(2),1:s0(3));
w4(1:s0(1),1:s0(2),1:s0(3)) = w{4}(1:s0(1),1:s0(2),1:s0(3));

if nargout == 1
    w1 = {w1 w2 w3 w4};
    w2 = []; w3 = []; w4 = [];
end