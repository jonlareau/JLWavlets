function [y,w,T] = JJL_den2D(x,T,J)
%%
%Denoising of a 2D data set using the complex dual-tree wavelet
%transform.
%[y,w] = JJL_den2D(x,T,J)
%%
if nargin == 0
    K = 64;
    D = 20;
    in = generateTestTarget();
    in = in+rand(size(in));
    imshow(in,[]); title('input');
    
    x = in;
    T = 1;
end

%%
%Setup parameters...
szx = size(x);
npow = nextpow2(max(szx));
x = padarray(x,max(0,2^(npow)-szx),'symmetric','post');
[Faf,Fsf] = FSfarras;
[af,sf] = dualfilt1;
if nargin < 3 || isempty(J)
    %According to my conversations with Dr. Farras & Dr. Selesnick, we
    %should try to ensure that the size of the highest scale Low-Frequency
    %component is >= Nx / 2^J.  Otherwise the filters will be longer than
    %the data vector at that level and errors will occur...
    J = floor(min(npow,max(log2(size(x)./size(af{1},1)))))
end
s0 = size(x)/2^(J);

%Compute the wavelet coefficients
w = JJL_cplxdtdwt2D(x,J,Faf,af);

%Set the threshold value
if nargin < 2 || isempty(T)
    T = std(w{1,1}(:))
end

%Filter the wavelet coefficients...
msk = true(size(w{1,1}));
msk(1:s0(1)-1,1:s0(2)-1) = 0;
for m = 1:2
    for n = 1:2
       w{m,n}(msk) = softThreshold(w{m,n}(msk),T);
    end
end
y = JJL_icplxdtdwt2D(w,J,Fsf,sf);

%Remove the padding.....
y = y(1:szx(1),1:szx(2));

if nargout == 0
    figure, imshow(y,[]); title('output');
    y = [];
    w = [];
    
end
