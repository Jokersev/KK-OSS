function [initial_object]= ER_realspace(diffraction_update,support_pad,vecX,ncX_big,vecY,ncY_big,paddingx,paddingy) 
%����ͼ�渵��Ҷ�任��ʵ�ռ�
initial_object = ifft2(diffraction_update);
%ʵ�ռ�Լ��
initial_object = initial_object(vecX + ncX_big, vecY + ncY_big);
initial_object(initial_object<0) = 0;   %С��0 �Ĳ���Ϊ0������Լ��
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %֧����Լ��
end




