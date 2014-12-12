function [out] = poly2JJL(w,type)
%Conversion utility to convert from the Brooklyn Polytechnic University
%version of wavelet coefficients to JJL version...
%
%w_JJL = poly2JJL(w_poly);
%
%THIS FUNCTION HAS NOT BEEN IMPLEMENTED FOR ALL WAVELET TYPES YET
%
if nargin < 2 || isempty(type)
    [type,J,sz] = detectPolyWaveletType(w);
else
    J = length(w)-1;
    sz = size(w{J+1}{1})*2^J;
end

s0 = sz/2^(J);
switch type
    case 'dwt1D'
        out = zeros(sz);
        s = max(sz);
        for k = 1:J
            out((s/2+1):s) = w{k}(:);
            s = s/2;
        end
        out(1:s) = w{J+1}(:);
    case 'dwt2D'
        out = zeros(sz);
        s = sz;
        for k = 1:J
            out((s(1)/2+1):s(1), 1:s(2)/2)           = w{k}{1}; %LH;
            out(1:s(1)/2,        (s(2)/2+1):s(2))    = w{k}{2}; %HL;
            out((s(1)/2+1):s(1), (s(2)/2+1):s(2))    = w{k}{3}; %HH;
            s = s/2;
        end
        out(1:s(1),1:s(2)) = w{J+1};
    case 'dwt3D'
        out = zeros(sz);
        s = sz;
        for k=1:J
            out((s(1)/2+1):s(1), 1:s(2)/2,           1:s(3)/2) = w{k}{1}; %LLH;
            out(1:s(1)/2,        (s(2)/2+1):s(2),    1:s(3)/2) = w{k}{2}; %LHL;
            out((s(1)/2+1):s(1), (s(2)/2+1):s(2),    1:s(3)/2) = w{k}{3}; %LHH;
            out(1:s(1)/2,        1:s(2)/2,           (s(3)/2+1):s(3)) = w{k}{4}; %HLL;
            out((s(1)/2+1):s(1), 1:s(2)/2,           (s(3)/2+1):s(3)) = w{k}{5}; %HLH;
            out(1:s(1)/2,        (s(2)/2+1):s(2),    (s(3)/2+1):s(3)) = w{k}{6}; %HHL;
            out((s(1)/2+1):s(1), (s(2)/2+1):s(2),    (s(3)/2+1):s(3)) = w{k}{7}; %HHH;
            s = s/2;
        end
        out(1:s(1),1:s(2),1:s(3)) = w{J+1};
    case 'dtdwt1D'
        out = dtunwrap(w,J,type,sz);
    case 'dtdwt2D'
        out = dtunwrap(w,J,type,sz);
    case 'dtdwt3D'
        out = dtunwrap(w,J,type,sz);
    case 'cmplxdtdwt1D'
            error('function not implemented yet for 1-d transforms');
    case 'cmplxdtdwt2D'
        out{1,1} = zeros(sz); %wr1-Real part of section 1 (orientations 1,2,&3)
        out{1,2} = zeros(sz); %wr2-Real part of section 2 (orientations 1,2,&3)
        out{2,1} = zeros(sz); %wi1-Imag part of section 1 (orientations 1,2,&3)
        out{2,2} = zeros(sz); %wi2-Imag part of section 2 (orientations 1,2,&3)
        s = sz;
        for k = 1:J %scale
            for m = 1:2 %real or imaginary parts
                for n=1:2 %tree section
                    out{m,n}((s(1)/2+1):s(1), 1:s(2)/2)           = w{k}{m}{n}{1}; %LH;
                    out{m,n}(1:s(1)/2,        (s(2)/2+1):s(2))    = w{k}{m}{n}{2}; %HL;
                    out{m,n}((s(1)/2+1):s(1), (s(2)/2+1):s(2))    = w{k}{m}{n}{3}; %HH;
                end
            end
            s = s/2;
        end
        for m = 1:2
            for n = 1:2
                out{m,n}(1:s(1),1:s(2)) = w{J+1}{m}{n}; %LL
            end
        end

    case 'cmplxdtdwt3D'
        s = sz;
        out = cell(2,2,2);
        for m =1:2
            for n =1:2
                for p = 1:2
                    out{m,n,p} = zeros(sz);
                end
            end
        end
        for k = 1:J %scale
            for m = 1:2 %Orientation
                for n=1:2 %Orientation
                    for p = 1:2 %Orientation
                        out{m,n,p}( (s(1)/2+1):s(1), 1:s(2)/2, 1:s(3)/2 )               = w{k}{m}{n}{p}{1}; %LH;
                        out{m,n,p}( 1:s(1)/2, (s(2)/2+1):s(2), 1:s(3)/2 )               = w{k}{m}{n}{p}{2}; %HL;
                        out{m,n,p}( (s(1)/2+1):s(1), (s(2)/2+1):s(2), 1:s(3)/2 )        = w{k}{m}{n}{p}{3}; %HH;
                        out{m,n,p}( 1:s(1)/2, 1:s(2)/2, (s(3)/2+1):s(3) )               = w{k}{m}{n}{p}{4}; %HLL;
                        out{m,n,p}( (s(1)/2+1):s(1), 1:s(2)/2, (s(3)/2+1):s(3) )        = w{k}{m}{n}{p}{5}; %HLH;
                        out{m,n,p}( 1:s(1)/2, (s(2)/2+1):s(2), (s(3)/2+1):s(3) )        = w{k}{m}{n}{p}{6}; %HHL;
                        out{m,n,p}( (s(1)/2+1):s(1), (s(2)/2+1):s(2), (s(3)/2+1):s(3) ) = w{k}{m}{n}{p}{7}; %HHH;
                    end
                end
            end
            s = s/2;
        end
        for m = 1:2
            for n = 1:2
                for p = 1:2
                    out{m,n,p}(1:s(1),1:s(2),1:s(3)) = w{J+1}{m}{n}{p}; %LL
                end
            end
        end
    otherwise
        error('not a known type');
end

function out = dtunwrap(w,J,type,sz)
%w is from a dual tree or dual tree complex wavelet transform....
s = sz;
nd = ndims(w{J+1}{1});
l = length(w{J+1});
out = [];
for i = 1:l
    out{i} = zeros(s);
end

for k = 1:J
    for j=1:l
        if nd == 3
            out{j}((s(1)/2+1):s(1), 1:s(2)/2,           1:s(3)/2) = w{k}{j}{1}; %LLH;
            out{j}(1:s(1)/2,        (s(2)/2+1):s(2),    1:s(3)/2) = w{k}{j}{2}; %LHL;
            out{j}((s(1)/2+1):s(1), (s(2)/2+1):s(2),    1:s(3)/2) = w{k}{j}{3}; %LHH;
            out{j}(1:s(1)/2,        1:s(2)/2,           (s(3)/2+1):s(3)) = w{k}{j}{4}; %HLL;
            out{j}((s(1)/2+1):s(1), 1:s(2)/2,           (s(3)/2+1):s(3)) = w{k}{j}{5}; %HLH;
            out{j}(1:s(1)/2,        (s(2)/2+1):s(2),    (s(3)/2+1):s(3)) = w{k}{j}{6}; %HHL;
            out{j}((s(1)/2+1):s(1), (s(2)/2+1):s(2),    (s(3)/2+1):s(3)) = w{k}{j}{7}; %HHH;
        elseif (nd == 2 && all(s~=1))
            %it is 2-d
            out{j}((s(1)/2+1):s(1), 1:s(2)/2)           = w{k}{j}{1}; %LLH;
            out{j}(1:s(1)/2,        (s(2)/2+1):s(2))    = w{k}{j}{2}; %LHL;
            out{j}((s(1)/2+1):s(1), (s(2)/2+1):s(2))    = w{k}{j}{3}; %LHH;
        elseif (nd == 2 && xor(s(1)==1,s(2)==1))
            %it is 1-d
            error('function not implemented yet for 1-d transforms');
        end
    end
    s = s/2;
end

for j = 1:l
        if nd == 3
            out{j}(1:s(1),1:s(2),1:s(3)) = w{J+1}{j};
        elseif (nd == 2 && all(s~=1))
            out{j}(1:s(1),1:s(2)) = w{J+1}{j};
        elseif (nd == 2 && xor(s(1)==1,s(2)==1))
            %it is 1-d
            error('function not implemented yet for 1-d transforms');
        end
end