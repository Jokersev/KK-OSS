%��Ʒ�ؽ���Ϻ������ȥ��͸�����ã����main2 �����������
function [rec,rec_lens] = reconstruction_obj6_convergence2_Fraunhoffer(x1,y1,Iterations,initial_field,initial_diffraction,initial_lens_field,initial_lens_diffraction,support_pad,paddingx,paddingy,am_p,am_p_lens,vecX,ncX_big,vecY,ncY_big,ph_p,ph_p_lens)
diffraction_update = initial_diffraction;
diffraction_lens_update = initial_lens_diffraction;
bestErr = 1e30;%reset best error
errK = zeros(1,Iterations,'single');
errK_lens = zeros(1,Iterations,'single');
% figure
for iterationNum = 1:Iterations 
    %͸�������ؽ�
    [initial_lens_field]= lens_realspace_Fraunhoffer(diffraction_lens_update,support_pad,paddingx,paddingy,iterationNum);
    [diffraction_lens_update]= lens_fourierspace(initial_lens_field,am_p_lens,vecX,ncX_big,vecY,ncY_big,iterationNum);
    rec_lens = ifft2(diffraction_lens_update);
%     figure;imagesc(x1,y1,abs(rec_lens));title(['lens reconstruction iterationNum = ',num2str(iterationNum)],'fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
    %monitor lens error
    k_lens = fft2(rec_lens);
%     %GENFIRE͸�����
    measuredK_lens = am_p_lens.*exp(1i*ph_p_lens);
    errInd_lens = find(measuredK_lens~=0);
    errK_lens(iterationNum) = sum(abs(k_lens(errInd_lens)-measuredK_lens(errInd_lens)))./sum(abs(measuredK_lens(errInd_lens))); 
%     �����ؽ�
    [initial_field]= realspace_Fraunhoffer(diffraction_update,paddingx,paddingy,support_pad,rec_lens,iterationNum);
     rec_initial_field = initial_field(vecX + ncX_big, vecY + ncY_big);
    % �����������е��ؽ�ͼ��
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
    %monitor error
    k_sample = fft2(rec_initial_field);
    k = k_lens + k_sample;
%     %GENFIRE���
    measuredK = am_p.*exp(1i*ph_p);
    errInd = find(measuredK~=0);
    errK(iterationNum) = sum(abs(k(errInd)-measuredK(errInd)))./sum(abs(measuredK(errInd))); 
%     %KCDI����������
%     errK(iterationNum) = sum((am_p(errInd)-k_am(errInd)).^2)./sum(am_p(errInd).^2); 
    if errK(iterationNum)<bestErr %if current reconstruction has better error, update best error and best reconstruction
%     fprintf('GENFIRE: new best object, iteration %d\n',iterationNum)
        bestErr = errK(iterationNum);
        rec = rec_initial_field;
    end
    fprintf('KCDI_ER: Iteration %d: Error = %d\n',iterationNum, errK(iterationNum));
    fprintf('KCDI_ER: Iteration %d: lens Error = %d\n',iterationNum, errK_lens(iterationNum)); 
end
%��ʾ�ؽ����
figure;imagesc(x1,y1,abs(rec));title('reconstruction result','fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% figure;plot(errK_lens);xlabel('iterationNum');ylabel('error');title('lens errK')
% figure;plot(errK);xlabel('iterationNum');ylabel('error');title('sample errK')
end