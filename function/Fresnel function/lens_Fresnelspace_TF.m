function [diffraction_update]= lens_Fresnelspace_TF(initial_object,L1,lambda,z,am_p,vecX,ncX_big,vecY,ncY_big,iterationNum) 
%ȥ��ʵ�ռ������
initial_object = initial_object(vecX + ncX_big, vecY + ncY_big);
initial_object(isnan(initial_object)) = 0;
%�������任��Ƶ��ռ�
diffraction=propTF(initial_object,L1,lambda,z);
%���׿ռ�Լ�������µ���λ��ʵ��ķ�ֵ���Ϊ����ͼ��������ͼ���£��ⳡ��ʽ
diffraction_update = am_p.*diffraction./abs(diffraction);
% figure
% % imshow(log(real(fftshift(diffraction_update))+1))
% imagesc(abs(diffraction_update));axis square;colormap('gray');
% title(['����ռ�����ͼ iterationNum = ',num2str(iterationNum)])
% %%monitor error
% errInd = find(pic_fft_pad~=0);
% errK(iterationNum) = sum(abs(pic_fft_pad(errInd)-diffraction_update(errInd)))./sum(abs(pic_fft_pad(errInd)));   
end




