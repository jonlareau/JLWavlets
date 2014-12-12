function [af, sf] = CDF97
% Cohen-Daubechies-Feauveau wavelet
%
% USAGE:
%    [af, sf] = CDF97
% OUTPUT:
%    af - analysis filters
%    sf - synthesis filters
% REFERENCE:
%   http://en.wikipedia.org/wiki/Cohen-Daubechies-Feauveau_wavelet
%
af =[0 0
    0.026748757411 0 
    -0.016864118443 0.091271763114 
    -0.078223266529 -0.057543526229 
    0.266864118443  -0.591271763114 
    0.602949018236  1.11508705 
    0.266864118443  -0.591271763114 
    -0.078223266529 -0.057543526229 
    -0.016864118443 0.091271763114 
    0.026748757411  0 ];

sf = [  0 0
        0               0.026748757411 
        -0.091271763114 0.016864118443 
        -0.057543526229 -0.078223266529 
        0.591271763114  -0.266864118443 
        1.11508705      0.602949018236 
        0.591271763114  -0.266864118443 
        -0.057543526229 -0.078223266529 
        -0.091271763114 0.016864118443 
        0               0.026748757411 ];

