function y = JJL_icplxdtdwt3D(w, J, Fsf, sf)

% Inverse 3D Complex Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   y = icplxdual3D(w, J, Fsf, sf)
% INPUT:
%   J - number of stages
%   Fsf - synthesis filter for last stage
%   sf - synthesis filters for preceeding stages
% OUTPUT:
%   y - output array
% See cplxdual3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

s0 = size(w{1,1,1})/2^(J);
msk = ones(size(w{1,1,1}));
msk(1:s0(1),1:s0(2),1:s0(3)) = 0;
msk = logical(msk);
[w{1,1,1}(msk) w{2,2,1}(msk) w{2,1,2}(msk) w{1,2,2}(msk)] = ...
    pm4inv(w{1,1,1}(msk), w{2,2,1}(msk), w{2,1,2}(msk), w{1,2,2}(msk));
[w{2,2,2}(msk) w{1,1,2}(msk) w{1,2,1}(msk) w{2,1,1}(msk)] = ...
    pm4inv(w{2,2,2}(msk), w{1,1,2}(msk), w{1,2,1}(msk), w{2,1,1}(msk));
    

%for j = 1:J
%    for m = 1:7
%        [w{j}{1}{1}{1}{m} w{j}{2}{2}{1}{m} w{j}{2}{1}{2}{m} w{j}{1}{2}{2}{m}] = ...
%            pm4inv(w{j}{1}{1}{1}{m}, w{j}{2}{2}{1}{m}, w{j}{2}{1}{2}{m}, w{j}{1}{2}{2}{m});
%         [w{j}{2}{2}{2}{m} w{j}{1}{1}{2}{m} w{j}{1}{2}{1}{m} w{j}{2}{1}{1}{m}] = ...
%            pm4inv(w{j}{2}{2}{2}{m}, w{j}{1}{1}{2}{m}, w{j}{1}{2}{1}{m}, w{j}{2}{1}{1}{m});
%    end
%end

y = zeros(size(w{1,1,1}));
SF = cell(2,3);
for m = 1:2
    for n = 1:2
        for p = 1:2
            SF{1,1} = Fsf{m};
            SF{1,2} = Fsf{n};
            SF{1,3} = Fsf{p};
            SF{2,1} = sf{m};
            SF{2,2} = sf{n};
            SF{2,3} = sf{p};
            
            y = y + JJL_idwt3D(w{m,n,p},J,SF);
            %{
            lo = w{J+1}{m}{n}{p};
            for j = J:-1:2
                lo = sfb3D(lo, w{j}{m}{n}{p}, sf{m}, sf{n}, sf{p});
            end
            lo = sfb3D(lo, w{1}{m}{n}{p}, Fsf{m}, Fsf{n}, Fsf{p});
            y = y + lo;
            %}
        end
    end
end

% normalization
y = y/sqrt(8);


