%样品重建完毕后最后再去掉透镜作用，针对main2 添加收敛条件
function [rec,rec_lens] = ER_reconstruction_Fraunhoffer2(x1,y1,Iterations,initial_field,initial_diffraction,initial_lens_field,initial_lens_diffraction,support_pad,supportlens_pad,paddingx,paddingy,am_p,am_p_lens,vecX,ncX_big,vecY,ncY_big,ph_p,ph_p_lens,L,lambda,z,lens_modelimg,modelimg)
diffraction_update = initial_diffraction;
diffraction_lens_update = initial_lens_diffraction;
bestErr = 1e30;%reset best error
errK = zeros(1,Iterations,'single');
errK_lens = zeros(1,Iterations,'single');
% Preallocate error arrays 预分配错误数组
RfacF_lens = zeros(ceil(Iterations/2),1,'single');  counter1_lens=0; errorF_lens=1;     %wy chagned
RfacR_lens = zeros(ceil(Iterations/2),1,'single');  counter2_lens=0; errorR_lens=1;     %wy chagned
% Preallocate error arrays 预分配错误数组
RfacF = zeros(ceil(Iterations/2),1,'single');  counter1=0; errorF=1;     %wy chagned
RfacR = zeros(ceil(Iterations/2),1,'single');  counter2=0; errorR=1;     %wy chagned
% figure
for iterationNum = 1:Iterations 
    %透镜衍射重建
    [initial_lens_field]= lens_realspace_Fraunhoffer2(diffraction_lens_update,supportlens_pad,L,lambda,z,paddingx,paddingy,iterationNum);
    [diffraction_lens_update]= lens_fourierspace2(initial_lens_field,am_p_lens,vecX,ncX_big,vecY,ncY_big,iterationNum,L,lambda,z);
    rec_lens = ifft2(diffraction_lens_update);
%     figure;imagesc(x1,y1,abs(rec_lens));title(['lens reconstruction iterationNum = ',num2str(iterationNum)],'fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');

%   透镜误差reciprocal space
    k_lens = fft2(rec_lens);
    errorF_lens = sum(sum(abs(k_lens(am_p_lens~=-1)-am_p_lens(am_p_lens~=-1)))) / sum(sum(am_p_lens(am_p_lens~=-1)));
    counter1_lens=counter1_lens+1; RfacF_lens(counter1_lens) = errorF_lens;
%   透镜误差real space
    errorR_lens= sum(sum(abs(lens_modelimg-rec_lens))) / sum(sum(abs(lens_modelimg+rec_lens)));  
    counter2_lens=counter2_lens+1; RfacR_lens(counter2_lens) = errorR_lens; 

   %     物体重建
    [initial_field]= realspace_Fraunhoffer(diffraction_update,paddingx,paddingy,support_pad,rec_lens,iterationNum);
     rec_initial_field = initial_field(vecX + ncX_big, vecY + ncY_big);
    % 画迭代过程中的重建图像
%     if mod(iterationNum,10)==0 
% %         figure;imagesc(x1,y1,abs(rec_lens));axis square;colormap('gray');
% %         title(['lens iterationNum = ',num2str(iterationNum)],'fontsize',22)
%         figure;imagesc(x1,y1,abs(rec_initial_field));axis square;colormap('gray');
%         title(['iterationNum = ',num2str(iterationNum)],'fontsize',22)    
%     else if iterationNum <= 10
% %              figure;imagesc(x1,y1,abs(rec_lens));axis square;colormap('gray');
% %              title(['lens iterationNum = ',num2str(iterationNum)],'fontsize',22)
%              figure;imagesc(x1,y1,abs(rec_initial_field));axis square;colormap('gray');
%              title(['iterationNum = ',num2str(iterationNum)],'fontsize',22)
%         end
%     end
    [diffraction_update]= fourierspace(initial_field,am_p,vecX,ncX_big,vecY,ncY_big,rec_lens,iterationNum);

%   物体误差reciprocal space
    k_sample = fft2(rec_initial_field);
    k = k_lens + k_sample;
    errorF = sum(sum(abs(k(am_p~=-1)-am_p(am_p~=-1)))) / sum(sum(am_p(am_p~=-1)));
    counter1=counter1+1; RfacF(counter1) = errorF;
%   物体误差real space
    errorR = sum(sum(abs(modelimg-rec_initial_field))) / sum(sum(abs(modelimg+rec_initial_field)));  
    counter2=counter2+1; RfacR(counter2) = errorR;   
    figure(63), plot(RfacF),xlabel('iterationNum/2');ylabel('rec errorF') , title(['errorF = ',int2str(errorF*100),'%']); %倒易空间误差%wy changed
    figure(64), plot(RfacR),xlabel('iterationNum');ylabel('rec errorR') , title(['errorR = ',mat2str(errorR,3)]); %实空间误差%wy changed
    
    if errK(iterationNum)<bestErr %if current reconstruction has better error, update best error and best reconstruction
%     fprintf('GENFIRE: new best object, iteration %d\n',iterationNum)
        bestErr = errK(iterationNum);
        rec = rec_initial_field;
    end
%     fprintf('KCDI_ER: Iteration %d: Error = %d\n',iterationNum, errK(iterationNum));
%     fprintf('KCDI_ER: Iteration %d: lens Error = %d\n',iterationNum, errK_lens(iterationNum)); 
end
%显示重建结果
figure;imagesc(x1,y1,abs(rec));title('reconstruction result','fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% figure;plot(errK_lens);xlabel('iterationNum');ylabel('error');title('lens errK')
% figure;plot(errK);xlabel('iterationNum');ylabel('error');title('sample errK')
end