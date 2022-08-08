function[u1]=propTF_inverse(u2,L,lambda,z);
% propagation - transfer function approach
% assumes same x and y side lengths and      xy长度相同
% uniform sampling       均匀采样
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

[M,N] = size(u2);  %get input field array size
dx = L/M;          %sample interval
k=2*pi/lambda;     %wavenumber
fx=-1/(2*dx):1/L:1/(2*dx)-1/L;   %freq coords
[FX,FY]=meshgrid(fx,fx);    %FX=FY

H = exp(-j*pi*lambda*z.*(FX.^2+FY.^2)); %trans func 其中exp(jkz)项被忽略，因为本项不影响观测平面结果的横向空间结构
H=fftshift(H);              %shift trans func
U2=fft2(fftshift(u2));      %shift,fft src field
U1=U2./H;                   %multiply
u1=ifftshift(ifft2(U1));    %inv fft,center obs field
end



