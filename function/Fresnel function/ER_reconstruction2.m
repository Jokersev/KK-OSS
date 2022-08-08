%�����������͸�������֧����Լ������Ʒ�ؽ���Ϻ������ȥ��͸�����ã����main2 �����������
function [rec,rec_lens] = ER_reconstruction2(x1,y1,Iterations,initial_field,initial_diffraction,initial_lens_field,initial_lens_diffraction,support_pad,supportlens_pad,paddingx,paddingy,L1,lambda,z2,z,am_p,am_p_lens,vecX,ncX_big,vecY,ncY_big,ph_p,ph_p_lens,lens_modelimg,modelimg)
diffraction_update = initial_diffraction;
diffraction_lens_update = initial_lens_diffraction;
bestErr = 1e30;%reset best error
errK = zeros(1,Iterations,'single');
errK_lens = zeros(1,Iterations,'single');
% Preallocate error arrays Ԥ�����������
RfacF_lens = zeros(ceil(Iterations/2),1,'single');  counter1_lens=0; errorF_lens=1;     %wy chagned
RfacR_lens = zeros(ceil(Iterations/2),1,'single');  counter2_lens=0; errorR_lens=1;     %wy chagned
% Preallocate error arrays Ԥ�����������
RfacF = zeros(ceil(Iterations/2),1,'single');  counter1=0; errorF=1;     %wy chagned
RfacR = zeros(ceil(Iterations/2),1,'single');  counter2=0; errorR=1;     %wy chagned
% figure
for iterationNum = 1:Iterations 
    %͸�������ؽ�
    [initial_lens_field]= lens_Fresnelrealspace_TF2(diffraction_lens_update,supportlens_pad,paddingx,paddingy,L1,lambda,z,iterationNum);
    [diffraction_lens_update]= lens_Fresnelspace_TF(initial_lens_field,L1,lambda,z,am_p_lens,vecX,ncX_big,vecY,ncY_big,iterationNum);
    rec_lens = propTF_inverse(diffraction_lens_update,L1,lambda,z2);
%     figure;imagesc(x1,y1,abs(rec_lens));title(['lens reconstruction iterationNum = ',num2str(iterationNum)],'fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');

%   ͸�����reciprocal space
    k_lens = abs(propTF(rec_lens,L1,lambda,z2));
    errorF_lens = sum(sum(abs(k_lens(am_p_lens~=-1)-am_p_lens(am_p_lens~=-1)))) / sum(sum(am_p_lens(am_p_lens~=-1)));
    counter1_lens=counter1_lens+1; RfacF_lens(counter1_lens) = errorF_lens;
%   ͸�����real space
    errorR_lens= sum(sum(abs(lens_modelimg-rec_lens))) / sum(sum(abs(lens_modelimg+rec_lens)));  
    counter2_lens=counter2_lens+1; RfacR_lens(counter2_lens) = errorR_lens;      
    
    %     �����ؽ�
    [initial_field]= Fresnelrealspace_TF(diffraction_update,paddingx,paddingy,L1,lambda,z2,support_pad,rec_lens,iterationNum);
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
    [diffraction_update]= Fresnelspace_TF(initial_field,L1,lambda,z2,am_p,vecX,ncX_big,vecY,ncY_big,rec_lens,iterationNum);
%   �������reciprocal space
    k_sample = propTF(rec_initial_field,L1,lambda,z2);
    k = k_lens + k_sample;
    errorF = sum(sum(abs(k(am_p~=-1)-am_p(am_p~=-1)))) / sum(sum(am_p(am_p~=-1)));
    counter1=counter1+1; RfacF(counter1) = errorF;
    
%   �������real space
    errorR = sum(sum(abs(modelimg-rec_initial_field))) / sum(sum(abs(modelimg+rec_initial_field)));  
    counter2=counter2+1; RfacR(counter2) = errorR;   
    
    figure(63), plot(RfacF),xlabel('iterationNum/2');ylabel('rec errorF') , title(['errorF = ',int2str(errorF*100),'%']); %���׿ռ����%wy changed
    figure(64), plot(RfacR),xlabel('iterationNum');ylabel('rec errorR') , title(['errorR = ',mat2str(errorR,3)]); %ʵ�ռ����%wy changed
    %if current reconstruction has better error, update best error and best reconstruction
    if errK(iterationNum)<bestErr 
%     fprintf('GENFIRE: new best object, iteration %d\n',iterationNum)
        bestErr = errK(iterationNum);
        rec = rec_initial_field;
    end
%     fprintf('KCDI_ER: Iteration %d: Error = %d\n',iterationNum, errorR(iterationNum));
%     fprintf('KCDI_ER: Iteration %d: lens Error = %d\n',iterationNum, errorR_lens(iterationNum)); 
end
%��ʾ�ؽ����
figure;imagesc(x1,y1,abs(rec));title('reconstruction result','fontsize',22);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% figure;plot(errK_lens);xlabel('iterationNum');ylabel('error');title('lens errK')
% figure;plot(errK);xlabel('iterationNum');ylabel('error');title('sample errK')
end