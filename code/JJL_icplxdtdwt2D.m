function y = JJL_icplxdtdwt2D(wr1,wi1,wr2,wi2, J, Fsf, sf)

% Inverse Dual-Tree Complex 2D Discrete Wavelet Transform
% 
%y = JJL_icplxdtdwt2D(wr1,wi1,wr2,wi2, J, Fsf, sf)
%y = JJL_icplxdtdwt2D(wr1,wi1,wr2,wi2, J)
%y = JJL_icplxdtdwt2D(w, J, Fsf, sf)
%y = JJL_icplxdtdwt2D(w, J)
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
if nargin == 2
    J = wi1;
    wi1 = wr1{2,1};
    wr2 = wr1{1,2};
    wi2 = wr1{2,2};
    wr1 = wr1{1,1};
    
    [Faf, Fsf] = FSfarras;
    [af, sf] = dualfilt1;
    
elseif nargin == 4
    J = wi1;
    Fsf = wr2;
    sf = wi2;
    wi1 = wr1{2,1};
    wr2 = wr1{1,2};
    wi2 = wr1{2,2};
    wr1 = wr1{1,1};
elseif nargin < 6
    
    [Faf, Fsf] = FSfarras;
    [af, sf] = dualfilt1;
end

%for j = 1:J
%   for m = 1:3
%        [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
%        [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
%    end
%end

s0 = size(wr1) / 2^(J);
[w{1,1},w{2,2}] = pm(wr1, wi2);
[w{1,2},w{2,1}] = pm(wr2, wi1);
w{1,1}(1:s0(1),1:s0(2)) = wr1(1:s0(1),1:s0(2));
w{2,1}(1:s0(1),1:s0(2)) = wi1(1:s0(1),1:s0(2));
w{1,2}(1:s0(1),1:s0(2)) = wr2(1:s0(1),1:s0(2));
w{2,2}(1:s0(1),1:s0(2)) = wi2(1:s0(1),1:s0(2));

y  = zeros(size(wr1));
%yc = y;
SF = cell(2,2);
for m = 1:2
    for n = 1:2
        
        SF{1,1} = Fsf{m};
        SF{1,2} = Fsf{n};
        SF{2,1} = sf{m};
        SF{2,2} = sf{n};
        y = y+JJL_idwt2D(w{m,n},J,SF);
        
        %lo = w{J+1}{m}{n};
        %for j = J:-1:2
        %    lo = sfb2D(lo, w{j}{m}{n}, sf{m}, sf{n});
        %end
        %lo = sfb2D(lo, w{1}{m}{n}, Fsf{m}, Fsf{n});
        %yc = yc + lo;
    end
end

% normalization
y = y/2;

