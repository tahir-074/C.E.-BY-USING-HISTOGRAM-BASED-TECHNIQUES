function y=periodogram2(x)
% Author: Wang Xianju
% April 1th, 2002
% This function computes the 2-D periodogram based on Fourier Method.
   
x=fft2(x);
x=x.*conj(x);
x=x/size(x,1)/size(x,2);
x=fftshift(x);
y=x;