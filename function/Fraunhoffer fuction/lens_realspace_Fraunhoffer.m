%ʵ�ռ�ȡ����֧����Լ��
function [initial_object]= lens_realspace_Fraunhoffer(initial_diffraction,support_pad,paddingx,paddingy,iterationNum) 
%����ͼ�渵��Ҷ�任��ʵ�ռ�
initial_object = ifft2(initial_diffraction);
%ʵ�ռ�Լ��
initial_object = padarray(initial_object,[paddingx paddingy]);
% initial_object = initial_object.*support_pad;  %֧����Լ��
% figure
% imagesc(abs(initial_object(vecX + ncX_big, vecY + ncY_big)));axis square;colormap('gray');
% title(['ʵ�ռ� iterationNum = ',num2str(iterationNum)])
end




