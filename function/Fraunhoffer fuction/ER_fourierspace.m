function [diffraction_update]= ER_fourierspace(initial_object,am_p_pad) 
%����Ҷ�任�����׿ռ�
diffraction = fft2(initial_object);
%���׿ռ�Լ�������µ���λ��ʵ��ķ�ֵ���Ϊ����ͼ��������ͼ����
diffraction_update = am_p_pad.*diffraction./abs(diffraction);
%figure
%imshow(log(real(fftshift(diffraction_update))+1))
%title(['����ͼ iterationNum = ',num2str(iterationNum)])
%%monitor error
%errInd = find(pic_fft_pad~=0);
%errK(iterationNum) = sum(abs(pic_fft_pad(errInd)-diffraction_update(errInd)))./sum(abs(pic_fft_pad(errInd)));   
end




