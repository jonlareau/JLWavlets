function x = fwt972D(x)
x = fwt97(x);
x = permute(fwt97(permute(x,[2,1])),[2,1]);