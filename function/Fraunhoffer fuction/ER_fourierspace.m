function [diffraction_update]= ER_fourierspace(initial_object,am_p_pad) 
%傅里叶变换到倒易空间
diffraction = fft2(initial_object);
%倒易空间约束：更新的相位与实验的幅值组合为衍射图，即衍射图更新
diffraction_update = am_p_pad.*diffraction./abs(diffraction);
%figure
%imshow(log(real(fftshift(diffraction_update))+1))
%title(['衍射图 iterationNum = ',num2str(iterationNum)])
%%monitor error
%errInd = find(pic_fft_pad~=0);
%errK(iterationNum) = sum(abs(pic_fft_pad(errInd)-diffraction_update(errInd)))./sum(abs(pic_fft_pad(errInd)));   
end




