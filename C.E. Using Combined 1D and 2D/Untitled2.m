clc;
clear all;
close all;
a1 = imread('06.jpg');
a = im2double(a1);
[m,n,k] = size(a);
ar = a(:,:,1);
ag = a(:,:,2);
ab = a(:,:,3);
br = 0.2*ar;
bg = 0.2*ag;
bb = 0.2*ab;
c = zeros(m,n,k);
c(:,:,1) = br;
c(:,:,2) = bg;
c(:,:,3) = bb;
figure,
imshow(c);

