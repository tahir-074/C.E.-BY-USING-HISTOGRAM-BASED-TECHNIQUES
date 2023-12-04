clc
h=ones(1,15)/15;
h=imrotate(h,46);
y=conv2(x,h);
imshow(uint8(y))
figure
Y=fft2(y);
[a b]=size(Y);
[a1 b1]=size(x);
X=fft2(x,a,b);
H1=Y./X;
H=fft2(h,a,b);
imshow(H)
figure
imshow(H1)
h1=ifft2(H);
h1=abs(real(h1(1:a-a1+1,1:b-b1+1)));
[a b]=size(h1);
%estimacija n i teta
for teta=0:180
    h2=imrotate(h1,-teta);
    a=sum(h2');
    a(a<0.0001)=0;
    if(sum(a>0)==1)
        'uspjelo'
        teta      %ima vise rjrsenja teta, zbog rotacije
        N=sum(sum(h2>0.001))
    end
end