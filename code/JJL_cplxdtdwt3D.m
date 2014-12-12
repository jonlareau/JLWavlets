function w = JJL_cplxdtdwt3D(x, J, Faf, af)

% 3D Complex Dual-Tree Discrete Wavelet Transform
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 0
    J = 4;
    x = ones(64,64,64); 
   [Faf, Fsf] = FSfarras;
   [af, sf] = dualfilt1;
end
    
% normalization
x = x/sqrt(8);

%Variable initialization
w = cell(2,2,2);
AF = cell(2,3);
for m = 1:2
    for n = 1:2
        for p = 1:2
            w{m,n,p} = zeros(size(x));
        end
    end
end

%For each orientation get the wavelet coefficients
for m = 1:2
    for n = 1:2
        for p = 1:2
            AF{1,1} = Faf{m};
            AF{1,2} = Faf{n};
            AF{1,3} = Faf{p};
            AF{2,1} = af{m};
            AF{2,2} = af{n};
            AF{2,3} = af{p};
            
            w{m,n,p}(:) = JJL_dwt3D(x,J,AF);
            
        end
    end
end

%Now get the sum and difference terms
s0 = size(x)/2^(J);
msk = ones(size(x));
msk(1:s0(1),1:s0(2),1:s0(3)) = 0;
msk = logical(msk);
[w{1,1,1}(msk) w{2,2,1}(msk) w{2,1,2}(msk) w{1,2,2}(msk)] = ...
    pm4(w{1,1,1}(msk), w{2,2,1}(msk), w{2,1,2}(msk), w{1,2,2}(msk));
[w{2,2,2}(msk) w{1,1,2}(msk) w{1,2,1}(msk) w{2,1,1}(msk)] = ...
    pm4(w{2,2,2}(msk), w{1,1,2}(msk), w{1,2,1}(msk), w{2,1,1}(msk));
    
%{
if nargin == 0
    %Show that this implementation is equivalent to the Poly implementation
    w1 = cplxdual3D(x*sqrt(8),J,Faf,af);
    y1 = icplxdual3D(w1,J,Fsf,sf);
    y = JJL_icplxdtdwt3D(w,J,Fsf,sf);
    
    m = abs(x*sqrt(8)-y);
    MSQE = max(m(:))
    m = abs(x*sqrt(8)-y1);
    MSQE = max(m(:))
end
%}



