clear all
close all

addpath(genpath('shearlet'));
[namefile,pathname]=uigetfile({'*.png;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.png)'},'Choose GrayScale Image');
[InImg_,map]=imread(namefile);
Si=ndims(InImg_);
if Si>2
    InImg=rgb2gray(InImg_);         
else
    InImg=InImg_;
end
img = InImg_;
figure,
imshow(img);
title('Original image');
cform = makecform('srgb2lab');
I = applycform(img,cform);
 HSV=rgb2hsv(img);
 figure;imshow(HSV);title('Histogram equalization (HE)');
 figure,
 rgb = hsv2rgb(HSV);
 surf(peaks);
 title('Mapping');
colorbar
figure,
imshow(I);
title('Histogram modification framework (HMF)');
L = I(:,:,1); a = I(:,:,2); b = I(:,:,3);
% % 
PyramidL = BuildPyramid(L);
Pyramida = BuildPyramid(a);
Pyramidb = BuildPyramid(b);
[N1,N2,N3] = size(PyramidL{1});
pl1 = PyramidL{1};
pl2 = imresize(PyramidL{2},[N1,N2]);
pl3 = imresize(PyramidL{3},[N1,N2]);

lmap = (0.4*pl2+0.6*pl3).*(pl1);
Thl = Threshold(lmap);
map1 = heaviside(lmap - Thl);

pa1 = Pyramida{1};
pa2 = imresize(Pyramida{2},[N1,N2]);
pa3 = imresize(Pyramida{3},[N1,N2]);

amap = (0.4*pa2+0.6*pa3).*(pa1);
Tha = Thl/1.5;
map2 = heaviside(amap - Tha);


pb1 = Pyramidb{1};
pb2 = imresize(Pyramidb{2},[N1,N2]);
pb3 = imresize(Pyramidb{3},[N1,N2]);

bmap = (0.4*pb2+0.6*pb3).*(pb1);
Thb = Tha/1.5;
map3 = heaviside(bmap - Thb);

Background = map1|map2|map3;

%%obtain the reflection component
se = strel('line',10,10);

w = fspecial('sobel');
for ch = 1:size(img,7)
    img = im2double(img);
    gx(:,:,ch) = imfilter(img(:,:,ch),w);
    gy(:,:,ch) = imfilter(img(:,:,ch),w');
    figure(5);
    imshow(img);
    colormap hot;
    title('Contextual and variational CE (CVC)');
    
end

grad = sqrt(gx.^2+gy.^2);
grad = max(grad,[],3);
reflectionPoints = find(grad<0.3&grad>0.05);
  map4 = zeros(N1,N2);
map4(reflectionPoints) = 1;

bw2 = map4;
reflectionpoints = find(bw2 == 1);
for i = 1:length(reflectionpoints)
    p = reflectionpoints(i);
    [Hp,rows,cols] = getpatch([N1,N2],p);
    o = Background(Hp);
    flag = find(o==1);
    if(Background(p)==1|length(flag)~=0)
        bw2(p) = 0;
    end
end

  Reflection = bw2;


indF = find(Reflection == 1);
indB = find(Background == 1);

[h,w,d] = size(grad);
G = struct;
[G.Gx,G.Gy,G.Gxx,G.Gxy,G.Gyy]=getGMat(w,h);

for ch = 1:3
    [I1(:,:,ch),I2(:,:,ch)]=layerSepIRLS (img(:,:,ch),G,indF,indB);
end

figure
imshow(I2);
title('Histogram-based Locality-Preserving CE (HBLPCE)');

L=imadjust(a);
 
%% Enhance contrast using histogram equalization 
H = histeq(a);

%% Brighten most of the details in the image except the nucleus
R1=imadd(L,H);
 
%% Highlight all the objects and its borders in the image including the cell nucleus
R2 = imsubtract(R1,H);

%% Remove almost all the other blood components while retaining the nucleus with minimum affect of distortion on the nucleus part of the white blood cell.
R3=imadd(R1,R2);
 
filterStart=tic;

for i=1:1
    R3=ordfilt2(R3,1,ones(3,3));
end


%=====================================================
level=graythresh(R3);
bw=im2bw(R3,level);
%Complement image
bw = imcomplement(bw);


out_mask = bw;

    %Re-count objects in the image
    cells = bwconncomp(bw);
    no_of_WBCs=cells.NumObjects;
   a1=InImg_;
    a2 = a1;
    for j=1:no_of_WBCs
        %return the coordinates for the pixels in object j
        [r, c] = find(bwlabel(bw)==j);
        rc = [r c];
        for i=1:max(size(rc))
            a2(rc(i,1),rc(i,2),:)=uint8(a1(rc(i,1),rc(i,2),:)*1.5);
        end
    end 
rgb = I2;
gray_image = rgb2gray(rgb);

[centers, radii] = imfindcircles(rgb,[35 60],'ObjectPolarity','dark','Sensitivity',0.9)

I = InImg_;
a = im2double(img);
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
 title('Spatial Entropy-based CE (SECE)');

    figure,imshow(InImg_);
    title('Residual Spatial Entropy-based CE (RSECE) ');
    figure
    RGB = hsv2rgb(HSV)
    imshow(RGB);
    title('SMIR');
    %show the final segmented image
    bw1 = imcomplement(bw);
    amask = a;
    amask(bw1) = 255; 

J = BIMEF(I); 
figure,
subplot 121; imshow(I); title('Original Image');
subplot 122; imshow(J); title('OUTPUT'); 



