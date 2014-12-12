function y = JJL_idtdwt3D(w1,w2,w3,w4,J,Fsf,sf)
%y = JJL_idtdwt3D(w1,w2,w3,w4,J,Fsf,sf)
%y = JJL_idtdwt3D(w,J,Fsf,sf)
% Inverse 3D Dual-Tree Discrete Wavelet Transform
%
%ADAPTED From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 4
    sf = w4; Fsf = w3; J = w2;
    w2 = w1{2}; w3 = w1{3}; w4 = w1{4}; w1 = w1{1};
end

%Invert the sum and difference erms
[w{1} w{2} w{3} w{4}] = pm4inv(w1,w2,w3,w4);

%...while keeping the low frequency information 
s0  = size(w1)/(2^(J));
w{1}(1:s0(1),1:s0(2),1:s0(3)) = w1(1:s0(1),1:s0(2),1:s0(3));
w{2}(1:s0(1),1:s0(2),1:s0(3)) = w2(1:s0(1),1:s0(2),1:s0(3));
w{3}(1:s0(1),1:s0(2),1:s0(3)) = w3(1:s0(1),1:s0(2),1:s0(3));
w{4}(1:s0(1),1:s0(2),1:s0(3)) = w4(1:s0(1),1:s0(2),1:s0(3));

M = [
    1 1 1
    2 2 1
    2 1 2
    1 2 2
];

% initialize output array
y = zeros(size(w{1}));

for i = 1:4
    f1 = M(i,1);
    f2 = M(i,2);
    f3 = M(i,3);
    SF = cell(2,3);
    [SF{1,1},SF{1,2},SF{1,3}] = deal(Fsf{f1}, Fsf{f2}, Fsf{f3});
    [SF{2,1},SF{2,2},SF{2,3}] = deal(sf{f1}, sf{f2}, sf{f3});
    
    y = y + JJL_idwt3D(w{i},J,SF);
end

% normalization
y = y/2;


