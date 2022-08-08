%实空间取消了支持域约束
function [initial_object]= lens_realspace_Fraunhoffer2(initial_diffraction,support_pad,L,lambda,z,paddingx,paddingy,iterationNum) 
%衍射图逆傅里叶变换到实空间
initial_object = ifft2(initial_diffraction);
initial_object = propTF_inverse(initial_object,L,lambda,z);
%实空间约束
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %支撑域约束
% figure
% imagesc(abs(initial_object(vecX + ncX_big, vecY + ncY_big)));axis square;colormap('gray');
% title(['实空间 iterationNum = ',num2str(iterationNum)])
end




