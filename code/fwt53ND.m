function x = fwt53ND(x,J,D)
nd = 1:ndims(x);
if nargin < 2
    J = 3;
end
if nargin < 3
    D = 1:ndims(x);
end
for i = D
    n  = nd;
    n(i) = 1;
    n(1) = i;
	x = permute(fwt53(permute(x,n),J),n);
end