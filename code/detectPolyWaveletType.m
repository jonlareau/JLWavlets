function [type,J,sz] = detectPolyWaveletType(w)
%This function was meant to help me write conversion utilities from the
%POLYT wavelet forms into standard Matrix formated wavelet. 
%
%The POLYT wavelets are not formatted in a way that can be speedily
%processed for various techniques, or that can handle larger matrices.
%Therfore I've been working to re-write their wavelet code so that the
%format of the wavelets can be used for better efficiency.
%
%   [type,J,sz] = detectPolyWaveletType(w)
J = length(w)-1;
if ~iscell(w{J+1})
    if iscell(w{1})
        if length(w{1}) == 3
            type = 'dwt2D';
            sz = size(w{J+1})*2^J;
        elseif length(w{1}) == 7
            type = 'dwt3D';
            sz = size(w{J+1})*2^J;
        end
    else
        type = 'dwt1D';
        sz = [length(w{J+1})*2^J,1];
    end
elseif ~iscell(w{J+1}{1})
    if iscell(w{1}{1})
        if length(w{1}{1}) == 3
            type = 'dtdwt2D';
            sz = size(w{J+1}{1})*2^J;
        elseif length(w{1}{1}) == 7
            type = 'dtdwt3D';
            sz = size(w{J+1}{1})*2^J;
        end
    else
        type = 'dtdwt1D';
        sz = length(w{J+1}{1})*2^J;
    end
elseif ~iscell(w{J+1}{1}{1})
    type = 'cmplxdtdwt2D';
    sz = size(w{J+1}{1}{1})*2^J;
elseif ~iscell(w{J+1}{1}{1}{1})
    type = 'cmplxdtdwt3D';
    sz = size(w{J+1}{1}{1}{1})*2^J;
else
    error('type unknown');
end
