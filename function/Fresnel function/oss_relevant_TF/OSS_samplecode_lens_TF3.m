% OSS Source code: �Ƿ��ڹ�����ʹ�����˲����ͽ���ѡ��

% function [mask RfacR RfacF RFD]= OSS_samplecode (F2D,supp,iter,beta,showim,modelimg,hiofirst) 
function [mask,RfacR,RfacF,sample]= OSS_samplecode_lens_TF3(F2D,supp,iter,beta,showim,modelimg,hiofirst,L,lambda,z,k_space) 
 
    %% General assignment of variables ������һ�㸳ֵ
    [Rsize,Csize] = size(F2D);
    R2D=zeros(Rsize,Csize,10,'single');
    toperrs=single(1:10:100);
    kfilter=zeros(Rsize,Csize,'single');
    realgroup=zeros(Rsize,Csize,2,'single');
    realgroup(:,:,1)=modelimg;

    %% Assign variables �������
    stopper = find(F2D==-1);
    filtercount=10;
%     filtnum=0;    %wy changed
    filtnum=1;
    store=0;

    %% Define support
    Rcenter = ceil(Rsize/2);
    Ccenter = ceil(Csize/2);
    Rsupport = supp(1);
    Csupport = supp(2);
    half_Rsupport = ceil(Rsupport/2);
    half_Csupport = ceil(Csupport/2);
    support = zeros(Rsize,Csize,'single');
    support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport-1,Ccenter-half_Csupport+1:Ccenter+half_Csupport-1) = 1;
    mask=support;

    %% Compute filter parameter alpha  ������������� alpha
    X=1:iter; 
    FX=(filtercount+1-ceil(X*filtercount/iter))*ceil(iter/(1*filtercount)); 
    FX=((FX-ceil(iter/filtercount))*(2*Rsize)/max(FX))+(2*Rsize/10);
    figure(98), plot(X,FX), axis tight; title('OSS Filter Size v. Iterations');

    %% Generate initial filter ���ó�ʼ������
    for kk=1:Rsize
        for jj=1:Csize
            kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(1)^2) ));
        end
    end
    kfilter=kfilter/max(max(kfilter));

    %% Assign random phases ���������λ
    rng('shuffle','twister');
    phase_angle = rand(Rsize,Csize,'single'); 

    %% Define initial k, r space ���ó�ʼ������ռ��ʵ�ռ�
%     initial_k=F2D; initial_k(initial_k==-1)=0;
%     k_space = initial_k.*exp(1i*phase_angle); %��ʼ����ͼ�������λ  
%     buffer_r_space = single(real(ifftn(ifftshift(k_space)))); %��ʼ���壬buffer����ʵ�ռ����� 
    buffer_r_space = propTF_inverse(k_space,L,lambda,z); %wy changed
    %% Preallocate error arrays Ԥ�����������
    RfacF = zeros(ceil(iter/2),1,'single');  counter1=0; errorF=1;     %wy chagned
    RfacR = zeros(ceil(iter/2),1,'single');  counter2=0; errorR=1;     %wy chagned

    %% Image display argument ͼ����ʾ����
    if showim==1     
        figure(1),
    end
    if nargin<7
        hiofirst=0;
    end

    %% OSS iterations
    for iteration = 1:iter

            %% OSS with Support & Positivity constraint
%             r_space = real(ifftn(ifftshift(k_space)));    %ʵ�ռ�
            r_space = propTF_inverse(k_space,L,lambda,z);
            sample = r_space.*mask;
            r_space = buffer_r_space-beta*r_space;
%             sample(sample<0)=r_space(sample<0);  %����Լ��

            %% Apply frequency filter (OSS)
            if hiofirst==0 || iteration>ceil(iter/filtercount)
                for kk=1:Rsize
                    for jj=1:Csize
                        kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(iteration)^2) ));
                    end
                end
                kfilter=kfilter/max(max(kfilter));
%                 ktemp=fftshift(fftn(r_space));
                ktemp=propTF(r_space,L,lambda,z);
                ktemp=ktemp.*kfilter;
%                 r_space=single(real(ifftn(ifftshift(ktemp))));  wy changed
                r_space=propTF_inverse(ktemp,L,lambda,z);
            end

            %% Use best result from last filter  ʹ�����һ������������ѽ�� wy changed
%             if mod(iteration,ceil(iter/filtercount))==0
%                 r_space=R2D(:,:,filtnum);  
%             else
%                 r_space(mask==1)=sample(mask==1);
%             end
%             
            % GENFIREE_OSS
            if mod(iteration,ceil(iter/filtercount))==0
                %start each new filter with the current best reconstruction
                %�Ե�ǰ��õ��ؽ���ʼ�µ��˲�
                r_space=bestRec; 
            else
                %otherwise put back the sample that was just isolated and  �������sample�Ż�ԭ�����Ҹ���
                %updated
                r_space(mask==1)=sample(mask==1);
            end

            %% Update reconstruction
            buffer_r_space = r_space;
%             k_space = fftshift(fftn(r_space));  %����ռ�
            k_space = propTF(r_space,L,lambda,z);
            phase_angle = angle(k_space);

            stopper_k_space = k_space(stopper);    
            k_space = F2D.*exp(1i*phase_angle);   
            k_space(stopper) = stopper_k_space;   

            %% Calculate errors    
            if rem(iteration,2)==0

                %% Calculate error in reciprocal space
                Ktemp = sample;
%                 Ktemp = abs(fftshift(fftn(Ktemp)));
                Ktemp = abs(propTF(Ktemp,L,lambda,z));
                
                errorF = sum(sum(abs(Ktemp(F2D~=-1)-F2D(F2D~=-1)))) / sum(sum(F2D(F2D~=-1)));
                
                counter1=counter1+1; RfacF(counter1) = errorF;

                %% Determine interations with best error  ȷ������������Ľ���
%                 filtnum=ceil(iteration*filtercount/iter);
%                 if errorF<= toperrs(filtnum) && iteration>store+2
%                     toperrs(filtnum)=errorF;
%                     R2D(:,:,filtnum)=r_space;
%                     store=iteration;
%                 end

                % GENFIREE_OSS
                %if current reconstruction has better error, update best error and best reconstruction
                %�����ǰ�ع��и��õ������������������ع�
                filtnum=ceil(iteration*filtercount/iter);
                if errorF<= toperrs(filtnum) && iteration>store+2 
                    bestRec = r_space;
                    store=iteration;
                end
                %% Calculate error in real space    
                realgroup(:,:,2)=sample;
                realgroup2=realgroup(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
%                 [realgroup2]=align2(realgroup2,0,iteration);
%                 errorR = sum(sum(abs(realgroup2(:,:,1)-realgroup2(:,:,2)))) / sum(sum(realgroup2(:,:,1)));  %wy changed
%                 counter2=counter2+1; RfacR(counter2) = errorR;     %wy changed
                errorR = sum(sum(abs(realgroup2(:,:,1)-realgroup2(:,:,2)))) / sum(sum(abs(realgroup2(:,:,1)+realgroup2(:,:,2))));  
                counter2=counter2+1; RfacR(counter2) = errorR;    

                %% Figure shows progress
                if showim==1
%                     figure(777),
%                     subplot(2,2,1), imagesc(abs(squeeze(realgroup2(:,:,1)))), axis square,colormap('gray'), title(strcat(int2str(FX(iteration)),'--OSS'));
%                     subplot(2,2,2), imagesc(abs(squeeze(realgroup2(:,:,2)))), axis square,colormap('gray'), title(int2str(iteration));
%                     subplot(2,2,3), plot(RfacF), axis([0 ceil(iteration/2) 0 0.8]), title(int2str(errorF*100)); %���׿ռ����
%                     subplot(2,2,4), plot(RfacR), axis([0 ceil(iteration/2) 0 0.8]), title(int2str(errorR*100)); %ʵ�ռ����
%                     figure(71), imagesc(abs(realgroup2(:,:,1))), axis square, colormap('gray');title(strcat(int2str(FX(iteration)),'--OSS')); %wy changed
%                     figure(72), imagesc(abs(realgroup2(:,:,2))), axis square, colormap('gray');title(['oss�ؽ� iteration = ',int2str(iteration)]);%wy changed
                    figure(73), plot(RfacF),xlabel('iterationNum/2');ylabel('lens errorF') , title(['errorF = ',int2str(errorF*100),'%']); %���׿ռ����%wy changed
                    figure(74), plot(RfacR),xlabel('iterationNum');ylabel('lens errorR') , title(['errorR = ',mat2str(errorR,3)]); %ʵ�ռ����%wy changed
                    
                    
                    
                    drawnow
                end
            end
    end

    %% Save results
%     if rem(iteration,iter)==0
%         save ('R2D.mat','R2D','RfacF','RfacR','mask','toperrs','r_space');
%     end

%     %% Show image: sum of best 4 steps
%     s=find(toperrs==min(toperrs));
%     RFD= squeeze(R2D(:,:,s)); 
% 
%     if showim==1
%     figure(2),
%     subplot(2,2,1), imagesc(squeeze(RFD)), axis image; 
%     subplot(2,2,2), imagesc(squeeze(mask)), axis image;
%     subplot(2,2,3), plot(RfacF), axis tight; 
%     subplot(2,2,4), plot(RfacR), axis tight;
%     end
end