function [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J,Faf,af,Fsf,sf,ALTERNATE)
%Function for making JJL version of the dualtree filters.  Formats either
%the default or user defined filters into JJL_[xxx]dtdwt[N]D() format.
%
%The default filters to use are the Farras designed filters from the POLY
%implementation.
%
%USAGE:
%   [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J)
%   [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J,Faf,af,Fsf,sf)
%   [haf,gaf,hsf,gsf] = JJL_makeDTFilters(J,Faf,af,Fsf,sf,ALTERNATE)
%
%INPUTS:
%   J - number of scales to use
%   Faf - First stage analysis filters
%   af - subsequent stage analysis filters
%   Fsf - First Stage Synthesis filters
%   sf - Subsequent Stage Synthesis filters
%   ALTERNATE - [(1),0] Should the subsequent stage filters be alternated
%       for stability?
%
%OUTPUTS:
%   haf - High Frequency analysis filters
%   gaf - Low Frequency analysis filters
%   hsf - High Frequency synthisis filters
%   hsf - Low Frequency synthisis filters

if nargin <2 || isempty(Faf) || isempty(af) %|| isempty(Fsf) || isempty(sf)
    [Faf, Fsf1] = FSfarras;
    [af, sf1] = dualfilt1;
end
if nargin < 6
    ALTERNATE = 1;
end

haf = cell(1,J);
gaf = cell(1,J);

%Analysis or filters
haf{1} = Faf{1};
gaf{1} = Faf{2};
if ~ALTERNATE
    for k = 2:J
        haf{k} = af{1};
        gaf{k} = af{2};
    end
else
    i = 2:2:J;
    j = 3:2:J;
    for k =i;
        haf{k} = af{1};
        gaf{k} = af{2};
    end
    for k = j
        haf{k}  = af{2};
        gaf{k}  = af{1};
    end
end
%}
if nargout >= 3 
    if  nargin < 4 || isempty(Fsf) || isempty(sf)
        [Faf1, Fsf] = FSfarras;
        [af1, sf] = dualfilt1;
    end
    
    hsf = cell(1,J);
    gsf = cell(1,J);
    %synthesis filters
    hsf{1} = Fsf{1};
    gsf{1} = Fsf{2};
    if ~ALTERNATE
        for k = 2:J
            hsf{k} = sf{1};
            gsf{k} = sf{2};
        end
    else
        i = 2:2:J;
        j = 3:2:J;
        for k =i;
            hsf{k}  = sf{1};
            gsf{k}  = sf{2};
        end
        for k = j
            hsf{k}  = sf{2};
            gsf{k}  = sf{1};
        end
    end
else
    hsf = [];
    gsf = [];
end