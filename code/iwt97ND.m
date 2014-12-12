function x = iwt97ND(x,J,D)
nd = 1:ndims(x);
if nargin < 2
    szx = size(x);
    J = min(floor(log2(szx(D))));
end
if nargin < 3
    D = ndims(x):-1:1;
end
for i = D
    n  = nd;
    n(i) = 1;
    n(1) = i;
    x = permute(x,n);
    x = iwt97(x,J);
    x = permute(x,n);
end