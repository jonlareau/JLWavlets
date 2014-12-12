function [y,w,T] = JJL_den3D(x,T,J)
%%
%Denoising of a 3D Volumetric data set using the complex dual-tree wavelet
%transform.
%[y,w] = JJL_den3D(x,T,J)
%%

%For Debugging the auto-thresholding setting
%global mnD;
%global mxD;
%global stD;
%global nnzD;


if nargin == 0
    K = 64;
    D = 20;
    in = .1*rand(K*[1,.6,2]);
    in(K/2-D:K/2+D,K/2-D:K/2+D,K/2-D:K/2+D) = rand((2*D+1)*[1,1,1]);
    
    in = in - mean(in(:));
    in = .5* in / max(abs(in(:)));
    
    mean(in(:))
    
    x = in;
    T = 1;
end

%%
%Setup parameters...
szx = size(x);
npow = nextpow2(min(szx));
x = padarray(x,max(0,2^(npow)-szx),'post');
[Faf,Fsf] = FSfarras;
[af,sf] = dualfilt1;
if nargin < 3 || isempty(J)
    %According to my conversations with Dr. Farras & Dr. Selesnick, we
    %should try to ensure that the size of the highest scale Low-Frequency
    %component is >= Nx / 2^J.  Otherwise the filters will be longer than
    %the data vector at that level and errors will occur...
    J = floor(min(npow,max(log2(size(x)./size(af{1},1)))));
end
s0 = size(x)/2^(J);

%Compute the wavelet coefficients
w = JJL_cplxdtdwt3D(x,J,Faf,af);

%Set the threshold value
if nargin < 2 || isempty(T)
    
    nnz = length(find(x~=0));
    T = std(w{1,1,1}(:));
    
    %For debugging the auto threshold generation
    %mnD(end+1) = min(w{1,1,1}(:));
    %mxD(end+1) = max(w{1,1,1}(:));
    %stD(end+1) = T;
    %nnzD(end+1) = length(find(x~=0));
end

%Filter the wavelet coefficients...
msk = true(size(w{1,1,1}));
msk(1:s0(1)-1,1:s0(2)-1,1:s0(3)-1) = 0;
for m = 1:2
    for n = 1:2
        for p = 1:2
            w{m,n,p}(msk) = softThreshold(w{m,n,p}(msk),T);
        end
    end
end
y = JJL_icplxdtdwt3D(w,J,Fsf,sf);

%Remove the padding.....
y = y(1:szx(1),1:szx(2),1:szx(3));

if nargout == 0
    y = [];
    w = [];
end
