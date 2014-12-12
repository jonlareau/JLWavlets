function y = softThreshold(x,T,mode)
%%
%y = softThreshold(x,T)
%
%Applies a soft threshold <T> to the values in <x>
%%
if nargin < 3
    mode = 'poly';
end

if strcmpi(mode,'poly')
    %Perform the soft thresholding as defined by the Brooklyn Polytechnic
    %University
    y = max(abs(x) - T, 0);
    y = y./(y+T).*x;
else
    y = max(abs(x)-T,0).*sign(x);
end