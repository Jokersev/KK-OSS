%实空间取消了支持域约束
function [initial_object]= lens_Fresnelrealspace_TF2(initial_diffraction,support_pad,paddingx,paddingy,L,lambda,z,iterationNum) 
%衍射图逆傅里叶变换到实空间
initial_object = propTF_inverse(initial_diffraction,L,lambda,z);
%实空间约束
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %支撑域约束
% figure
% imagesc(abs(initial_object));axis square;colormap('gray');
% title(['实空间 iterationNum = ',num2str(iterationNum)])
end




