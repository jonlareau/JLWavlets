function [out] = JJL_dwt3D(x, J, af)
% 3-D Discrete Wavelet Transform
%
%
%Adapted From:
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
s = size(x);
out = x;
for k = 1:J
    if iscell(af)
        [sr, sc] = size(af);
        if (sr <= 1) || (sc <= 1)
            afk1 = af{min(k,length(af))};
            afk2 = af{min(k,length(af))};
            afk3 = af{min(k,length(af))};
        else
            afk1 = af{min(k,sr),1};
            afk2 = af{min(k,sr),2};
            afk3 = af{min(k,sr),3};
        end
    else
        afk1 = af; afk2 = af; afk3 = af;
    end
    [out(1:s(1),1:s(2),1:s(3))] = JJL_afb3D(out(1:s(1),1:s(2),1:s(3)), afk1, afk2, afk3);
    s = s/2;
end

