%ʵ�ռ�ȡ����֧����Լ��
function [initial_object]= lens_Fresnelrealspace_TF2(initial_diffraction,support_pad,paddingx,paddingy,L,lambda,z,iterationNum) 
%����ͼ�渵��Ҷ�任��ʵ�ռ�
initial_object = propTF_inverse(initial_diffraction,L,lambda,z);
%ʵ�ռ�Լ��
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %֧����Լ��
% figure
% imagesc(abs(initial_object));axis square;colormap('gray');
% title(['ʵ�ռ� iterationNum = ',num2str(iterationNum)])
end




