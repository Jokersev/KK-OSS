% ����ͼ�м�ȥ͸��������任����ʵ�ռ�Լ��
function [initial_object]= Fresnelrealspace_TF(initial_diffraction,paddingx,paddingy,L,lambda,z,support_pad,ulens,iterationNum) 
%����ͼ�м�ȥ͸����
ulens_diffraction = propTF(ulens,L,lambda,z);
initial_diffraction = initial_diffraction-ulens_diffraction;
%����ͼ�渵��Ҷ�任��ʵ�ռ�
initial_object = propTF_inverse(initial_diffraction,L,lambda,z);
%ʵ�ռ�Լ��
% initial_object(initial_object<0) = 0;   %С��0 �Ĳ���Ϊ0������Լ��
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %֧����Լ��
% figure
% imagesc(abs(initial_object));axis square;colormap('gray');
% title(['ʵ�ռ� iterationNum = ',num2str(iterationNum)])
end




