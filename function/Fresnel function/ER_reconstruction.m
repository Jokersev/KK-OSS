%在透镜面进行支持域约束，样品重建完毕后最后再去掉透镜作用，针对main2 添加收敛条件
function [rec,rec_lens] = ER_reconstruction(x1,y1,Iterations,initial_field,initial_diffraction,initial_lens_field,initial_lens_diffraction,support_pad,supportlens_pad,paddingx,paddingy,L1,lambda,z2,z,am_p,am_p_lens,vecX,ncX_big,vecY,ncY_big,ph_p,ph_p_lens)
diffraction_update = initial_diffraction;
diffraction_lens_update = initial_lens_diffraction;
bestErr = 1e30;%reset best error
errK = zeros(1,Iterations,'single');
errK_lens = zeros(1,Iterations,'single');
% figure
for iterationNum = 1:Iterations 
    %透镜衍射重建
    [initial_lens_field]= lens_Fresnelrealspace_TF2(diffraction_lens_update,supportlens_pad,paddingx,paddingy,L1,lambda,z,iterationNum);
    [diffraction_lens_update]= lens_Fresnelspace_TF(initial_lens_field,L1,lambda,z,am_p_lens,vecX,ncX_big,vecY,ncY_big,iterationNum);
    rec_lens = propTF_inverse(diffraction_lens_update,L1,lambda,z2);
%     figure;imagesc(x1,y1,abs(rec_lens));title(['lens reconstruction iterationNum = ',num2str(iterationNum)],'fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
    %monitor lens error
    k_lens = propTF(rec_lens,L1,lambda,z2);
%     %GENFIRE透镜误差
    measuredK_lens = am_p_lens.*exp(1i*ph_p_lens);
    errInd_lens = find(measuredK_lens~=0);
    errK_lens(iterationNum) = sum(abs(k_lens(errInd_lens)-measuredK_lens(errInd_lens)))./sum(abs(measuredK_lens(errInd_lens))); 
%     物体重建
    [initial_field]= Fresnelrealspace_TF(diffraction_update,paddingx,paddingy,L1,lambda,z2,support_pad,rec_lens,iterationNum);
     rec_initial_field = initial_field(vecX + ncX_big, vecY + ncY_big);
    % 画迭代过程中的重建图像
    if mod(iterationNum,10)==0 
%         figure;imagesc(x1,y1,abs(rec_lens));axis square;colormap('gray');
%         title(['lens iterationNum = ',num2str(iterationNum)],'fontsize',22)
        figure;imagesc(x1,y1,abs(rec_initial_field));axis square;colormap('gray');
        title(['iterationNum = ',num2str(iterationNum)],'fontsize',22)    
    else if iterationNum <= 10
%              figure;imagesc(x1,y1,abs(rec_lens));axis square;colormap('gray');
%              title(['lens iterationNum = ',num2str(iterationNum)],'fontsize',22)
             figure;imagesc(x1,y1,abs(rec_initial_field));axis square;colormap('gray');
             title(['iterationNum = ',num2str(iterationNum)],'fontsize',22)
        end
    end
    [diffraction_update]= Fresnelspace_TF(initial_field,L1,lambda,z2,am_p,vecX,ncX_big,vecY,ncY_big,rec_lens,iterationNum);
    %monitor error
    k_sample = propTF(rec_initial_field,L1,lambda,z2);
    k = k_lens + k_sample;
%     %GENFIRE误差
    measuredK = am_p.*exp(1i*ph_p);
    errInd = find(measuredK~=0);
    errK(iterationNum) = sum(abs(k(errInd)-measuredK(errInd)))./sum(abs(measuredK(errInd))); 
%     %KCDI文献里的误差
%     errK(iterationNum) = sum((am_p(errInd)-k_am(errInd)).^2)./sum(am_p(errInd).^2); 
    if errK(iterationNum)<bestErr %if current reconstruction has better error, update best error and best reconstruction
%     fprintf('GENFIRE: new best object, iteration %d\n',iterationNum)
        bestErr = errK(iterationNum);
        rec = rec_initial_field;
    end
    fprintf('KCDI_ER: Iteration %d: Error = %d\n',iterationNum, errK(iterationNum));
    fprintf('KCDI_ER: Iteration %d: lens Error = %d\n',iterationNum, errK_lens(iterationNum)); 
end
%显示重建结果
figure;imagesc(x1,y1,abs(rec));title('reconstruction result','fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% figure;plot(errK_lens);xlabel('iterationNum');ylabel('error');title('lens errK')
% figure;plot(errK);xlabel('iterationNum');ylabel('error');title('sample errK')
end