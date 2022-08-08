function [diffraction_update]= lens_fourierspace2(initial_object,am_p,vecX,ncX_big,vecY,ncY_big,iterationNum,L,lambda,z) 
%去掉实空间零填充
initial_object = initial_object(vecX + ncX_big, vecY + ncY_big);
initial_object(isnan(initial_object)) = 0;
%傅里叶变换到频域空间
initial_object = propTF(initial_object,L,lambda,z);
diffraction=fft2(initial_object);
%倒易空间约束：更新的相位与实验的幅值组合为衍射图，即衍射图更新，光场形式
diffraction_update = am_p.*diffraction./abs(diffraction);
% figure
% % imshow(log(real(fftshift(diffraction_update))+1))
% imagesc(abs(diffraction_update));axis square;colormap('gray');
% title(['衍射空间衍射图 iterationNum = ',num2str(iterationNum)])
% %%monitor error
% errInd = find(pic_fft_pad~=0);
% errK(iterationNum) = sum(abs(pic_fft_pad(errInd)-diffraction_update(errInd)))./sum(abs(pic_fft_pad(errInd)));   
end




