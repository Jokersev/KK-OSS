% 衍射图中减去透镜场再逆变换进行实空间约束
function [initial_object]= realspace_Fraunhoffer(initial_diffraction,paddingx,paddingy,support_pad,ulens,iterationNum) 
%衍射图中减去透镜场
ulens_diffraction = fft2(ulens);
initial_diffraction = initial_diffraction-ulens_diffraction;
%衍射图逆傅里叶变换到实空间
initial_object = ifft2(initial_diffraction);
%实空间约束
% initial_object(initial_object<0) = 0;   %小于0 的部分为0，正向约束
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %支撑域约束
% figure
% imagesc(abs(initial_object));axis square;colormap('gray');
% title(['实空间 iterationNum = ',num2str(iterationNum)])
end




