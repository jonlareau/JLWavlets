% cplxdual2D_plots
% DISPLAY 2D WAVELETS OF cplxdual2D.M

J = 4;
L = 3*2^(J+1);
N = L/2^J;
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
x = zeros(2*L,6*L);
szx = size(x)/2^J;

if 1
    w = JJL_cplxdtdwt2D(x,J,Faf,af);
    
    w{1,1}(3855) = 1;
    w{1,1}(8463) = 1;
    w{1,1}(13059) = 1;
    
    w{1,2}(2703) = 1;
    w{1,2}(7299) = 1;
    w{1,2}(11919) = 1;

    w{2,1}(3861) = 1;
    w{2,1}(8469) = 1;
    w{2,1}(13065) = 1;
    
    w{2,2}(2709) = 1;
    w{2,2}(7305) = 1;
    w{2,2}(11925) = 1;
    
    y = JJL_icplxdtdwt2D(w,J,Fsf,sf);
else
    %Use the POLY code to validate...
    w = cplxdual2D(x, J, Faf, af);
    w{J}{1}{1}{1}(N/2,      N/2+3*N) = 1;
    w{J}{1}{1}{2}(N/2,      N/2+5*N) = 1;
    w{J}{1}{1}{3}(N/2,      N/2+1*N) = 1;
    w{J}{1}{2}{1}(N/2,      N/2+2*N) = 1;
    w{J}{1}{2}{2}(N/2,      N/2+0*N) = 1;
    w{J}{1}{2}{3}(N/2,      N/2+4*N) = 1;
    w{J}{2}{1}{1}(N/2+N,    N/2+3*N) = 1;
    w{J}{2}{1}{2}(N/2+N,    N/2+5*N) = 1;
    w{J}{2}{1}{3}(N/2+N,    N/2+1*N) = 1;
    w{J}{2}{2}{1}(N/2+N,    N/2+2*N) = 1;
    w{J}{2}{2}{2}(N/2+N,    N/2+0*N) = 1;
    w{J}{2}{2}{3}(N/2+N,    N/2+4*N) = 1;
    y = icplxdual2D(w, J, Fsf, sf);
end

y = [y; sqrt(y(1:L,:).^2+y(L+[1:L],:).^2)];
figure(1)
clf
imagesc(y);
title('2D Dual-Tree Complex Wavelets')
axis image
axis off
colormap(gray(128))
print -djpeg95 cplxdual2D_plots