% OSS Source code:

function [mask,RfacR,RfacF,sample]= OSS_samplecode3(F2D,supp,iter,beta,showim,modelimg,hiofirst,k_space) 
 
%% General assignment of variables
[Rsize,Csize] = size(F2D);
R2D=zeros(Rsize,Csize,10,'single');
toperrs=single(1:10:100);
kfilter=zeros(Rsize,Csize,'single');
realgroup=zeros(Rsize,Csize,2,'single');
realgroup(:,:,1)=modelimg;
 
%% Assign variables
stopper = find(F2D==-1);
filtercount=10;   %十种滤波器
% filtnum=0;   %wy changed
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
 
%% Compute filter parameter alpha 
X=1:iter; 
FX=(filtercount+1-ceil(X*filtercount/iter))*ceil(iter/(1*filtercount)); 
FX=((FX-ceil(iter/filtercount))*(2*Rsize)/max(FX))+(2*Rsize/10);
figure(98), plot(X,FX), axis tight; title('OSS Filter Size v. Iterations');
 
%% Generate initial filter
for kk=1:Rsize
    for jj=1:Csize
        kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(1)^2) ));
    end
end
kfilter=kfilter/max(max(kfilter));
 
%% Assign random phases
rng('shuffle','twister');
phase_angle = rand(Rsize,Csize,'single'); 
 
%% Define initial k, r space
% initial_k=F2D; initial_k(initial_k==-1)=0;     %wy changed
% k_space = initial_k.*exp(1i*phase_angle);      %wy changed
% buffer_r_space = single(real(ifftn(ifftshift(k_space))));    %wy changed
buffer_r_space = single(ifft2(k_space));   %返回实空间此时的图像应该是复值，所以不能只取实部real

%% Preallocate error arrays
RfacF = zeros(ceil(iter/2),1,'single');  counter1=0; errorF=1;
RfacR = zeros(ceil(iter/2),1,'single');  counter2=0; errorR=1;
 
%% Image display argument
if showim==1    
    figure(1),
end
if nargin<7
    hiofirst=0;
end
 
%% OSS iterations
for iteration = 1:iter
    
        %% OSS with Support & Positivity constraint
%         r_space = real(ifftn(ifftshift(k_space)));  %wy changed
        r_space = ifft2(k_space);
        sample = r_space.*mask;
        r_space = buffer_r_space-beta*r_space;
%         sample(sample<0)=r_space(sample<0);       %正则约束，先保留
        
        %% Apply frequency filter (OSS)
        if hiofirst==0 || iteration>ceil(iter/filtercount)
            for kk=1:Rsize
                for jj=1:Csize
                    kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(iteration)^2) ));
                end
            end
            kfilter=kfilter/max(max(kfilter));  %按照迭代次数的情况选择相应的滤波器
%             ktemp=fftshift(fftn(r_space)); %wy changed
            ktemp=fft2(r_space);  %支持域外的部分的傅里叶变换
            ktemp=ktemp.*kfilter;  %支持域外的部分应用滤波器
%             r_space=single(real(ifftn(ifftshift(ktemp))));    %wy chenged
            r_space=single(ifft2(ktemp));   %滤波之后支持域外的部分返回实空间
        end
 
        %% Use best result from last filter 从最后的滤波器中选择支持域外部的最佳结果
%         if mod(iteration,ceil(iter/filtercount))==0
%             r_space=R2D(:,:,filtnum);  
%         else
%             r_space(mask==1)=sample(mask==1);
%         end

        %% Update reconstruction
        buffer_r_space = r_space;
%         k_space = fftshift(fftn(r_space));      %wy changed
        k_space = fft2(r_space);  
        phase_angle = angle(k_space);
  
        stopper_k_space = k_space(stopper);    
        k_space = F2D.*exp(1i*phase_angle);   
        k_space(stopper) = stopper_k_space;   
                   
        %% Calculate errors    
        if rem(iteration,2)==0   %偶数次迭代计算误差
            
            %% Calculate error in reciprocal space
            Ktemp = sample;
%             Ktemp = abs(fftshift(fftn(Ktemp)));        %wy changed
            Ktemp = abs(fft2(Ktemp));
            errorF = sum(sum(abs(Ktemp(F2D~=-1)-F2D(F2D~=-1)))) / sum(sum(F2D(F2D~=-1)));
            counter1=counter1+1; RfacF(counter1) = errorF;
            
            %% Determine interations with best error 确定具有最佳误差的交互作用
%             filtnum=ceil(iteration*filtercount/iter);
%             if errorF<= toperrs(filtnum) && iteration>store+2
%                 toperrs(filtnum)=errorF;
%                 R2D(:,:,filtnum)=r_space;
%                 store=iteration;
%             end
            
            %% Calculate error in real space    
            realgroup(:,:,2)=sample;
            realgroup2=realgroup(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
            [realgroup2]=align2(realgroup2,0,iteration);

            %wy changed
%             errorR = sum(sum(abs(realgroup2(:,:,1)-realgroup2(:,:,2)))) / sum(sum(realgroup2(:,:,1)));
%             counter2=counter2+1; RfacR(counter2) = errorR;
             errorR = sum(sum(abs(realgroup2(:,:,1)-realgroup2(:,:,2)))) / sum(sum(abs(realgroup2(:,:,1)+realgroup2(:,:,2))));  
             counter2=counter2+1; RfacR(counter2) = errorR;     
              
            %% Figure shows progress
            if showim==1
                %wy changed
                figure (65),     
                subplot(2,2,1), imagesc(squeeze(realgroup2(:,:,1))), axis image, title(strcat(int2str(FX(iteration)),'--OSS'));
                subplot(2,2,2), imagesc(squeeze(realgroup2(:,:,2))), axis image, title(int2str(iteration));
                subplot(2,2,3), plot(RfacF), axis([0 ceil(iteration/2) 0 0.8]), title(int2str(errorF*100)); 
                subplot(2,2,4), plot(RfacR), axis([0 ceil(iteration/2) 0 0.8]), title(int2str(errorR*100)); 
%                 figure(61), imagesc(abs(realgroup2(:,:,1))), axis square, colormap('gray');title(strcat(int2str(FX(iteration)),'--OSS')); %wy changed
%                 figure(62), imagesc(abs(realgroup2(:,:,2))), axis square, colormap('gray');title(['oss重建 iteration = ',int2str(iteration)]);%wy changed
%                 figure(63), plot(RfacF),xlabel('iterationNum/2');ylabel('rec errorF') , title(['errorF = ',int2str(errorF*100),'%']); %倒易空间误差%wy changed
%                 figure(64), plot(RfacR),xlabel('iterationNum');ylabel('rec errorR') , title(['errorR = ',mat2str(errorR,3)]); %实空间误差%wy changed

                drawnow
            end
        end
end
 
%% Save results     %wy changed
% if rem(iteration,iter)==0
%     save ('R2D.mat','R2D','RfacF','RfacR','mask','toperrs','r_space');
% end
 
%% Show image: sum of best 4 steps   %wy changed
% s=find(toperrs==min(toperrs));
% RFD= squeeze(R2D(:,:,s)); 
%  
% if showim==1
% figure(2),
% subplot(2,2,1), imagesc(squeeze(RFD)), axis image; 
% subplot(2,2,2), imagesc(squeeze(mask)), axis image;
% subplot(2,2,3), plot(RfacF), axis tight; 
% subplot(2,2,4), plot(RfacR), axis tight;
end
