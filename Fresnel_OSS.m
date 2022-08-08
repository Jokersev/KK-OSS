%% KCDI_oss ʹ�÷����������µ�oss
%% 
close all
clear
%% ����KCDIʵ������
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDsample_intensity_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDsample_angle_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDlens_intensity_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDlens_angle_IR.mat')
% load('F:\users\wy\KCDI_ER\oss data\KCDI_sample_supportfield.mat')
% load('F:\users\wy\KCDI_ER\oss data\KCDI_lens_supportfield.mat')
load('KCDI_CCDsample_intensity_TF.mat')
load('KCDI_CCDsample_angle_TF.mat')
load('KCDI_CCDlens_intensity_TF.mat')
load('KCDI_CCDlens_angle_TF.mat')
load('KCDI_lens_supportfield_TF.mat')
load('KCDI_sample_supportfield_TF.mat')
I=Intensity;
ph_p_orig = Angle;
I_lens = lens_Intensity;
ph_p_orig_lens = lens_Angle;

%% ��ʾԭ����
% figure
% subplot(2,2,1);imagesc(I);axis square;colormap('gray');
% subplot(2,2,2);imagesc(ph_p_orig);axis square;
% subplot(2,2,3);imagesc(I_lens);axis square;colormap('gray');
% subplot(2,2,4);imagesc(ph_p_orig_lens);axis square;
%% ��������
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% % export_fig(gcf,'-eps','-r300','-painters','./δ��������.eps');
% % figure;imagesc(I_lens);axis square;colormap('gray');title('��������ǰ')
[M,N] = size(I);
lam = 0.1; %����5%=0.05����
r = poissrnd(lam,[M,N]); %���Ӳ�������
I = I+I.*r;
I_lens = I_lens+I_lens.*r;
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./����10%������.eps');
figure;imagesc(I_lens);axis square;colormap('gray');title('͸������������')
%% ���ò���
L1=28672e-6; 
M=2048;       %number of samples
lambda=632e-9;     %wavelength(m)
z1=0.601;             %ͨ��͸�������������601mm������Ʒ
z2=0.463;             %��ccd�ķ������������463mm
zf=0.474;           %������۵ľ��루���ࣩ474mm
wl=400;    %����ֵ͸���ھ�(����Ϊ��λ)
%%
dx1=L1/M;    %src sample interval
x1=-L1/2:dx1:L1/2-dx1;    %src coords
y1=x1;
[X1,Y1]=meshgrid(x1,y1);
k=2*pi/lambda;      %wavenumber
[x_array,y_array] = meshgrid(1:M,1:M); 
x_array = x_array - floor(max(x_array(:))/2+1); % center of image to be zero 
y_array = y_array - floor(max(y_array(:))/2+1); % center of image to be zero 
%����͸������z1��������䳡
u1=(x_array./wl).^2+(y_array./wl).^2 <= 1; 
uout=u1.*exp(-1i*k/(2*zf)*(X1.^2+Y1.^2));
ulens_1=propTF(uout,L1,lambda,z1); 
%����͸�����ʺ̷Ѵ�����CCD�����䳡
ulens_2=fftshift(fft2(ulens_1)); 
% figure
% imagesc(x1,y1,I)
% title('CCD intensity')
% axis square;xlabel('x/m');ylabel('y/m');colormap('gray')
%% ER�㷨
%% ��ʼ�������ֵ
siza = size(I);                                     
am_p = sqrt(I);
am_p_lens = sqrt(I_lens);
%͸���ĳ�ʼ��λ
ph_p_lens = rand(siza)*2*pi;% rand(N(1),N(1))*pi; %�����λ
% ph_p_lens = ph_p_orig_lens; %͸����ԭͼ��λ
% ph_p_lens =imnoise(ph_p_orig_lens,'salt & Pepper');  %ԭͼ��λ�ӽ�������
% ph_p_lens =imnoise(ph_p_orig_lens,'gaussian',0.4);  %ԭͼ��λ�Ӹ�˹���� 
% ph_p_lens =imnoise(ph_p_orig_lens,'poisson');  %ԭͼ��λ�Ӳ������� 
%����ĳ�ʼ��λ
ph_p = rand(siza)*2*pi;% rand(N(1),N(1))*pi; %�����λ
% ph_p = ph_p_orig; %�����ԭͼ��λ
% ph_p = imnoise(ph_p_orig,'salt & Pepper');  %ԭͼ��λ�ӽ�������
% ph_p = imnoise(ph_p_orig,'gaussian',0.4);  %ԭͼ��λ�Ӹ�˹���� 
% ph_p = imnoise(ph_p_orig,'poisson');  %ԭͼ��λ�Ӳ������� 
if_kk_lens = 0;
if_kk = 0;
tic;
[initial_lens_field,initial_lens_diffraction,ws] = initial_guess(x1,y1,if_kk_lens ,ulens_2,I_lens,L1,lambda,M,ph_p_lens,am_p_lens,siza,z2);
[initial_field,initial_diffraction,ws] = initial_guess(x1,y1,if_kk ,ulens_2,I,L1,lambda,M,ph_p,am_p,siza,z2);
% figure
% subplot(2,2,1)
% imagesc(x1,y1,abs(initial_lens_field))
% title('initial guess lens field','fontsize',18);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% subplot(2,2,2)
% imagesc(x1,y1,abs(initial_lens_diffraction))
% title('initial guess lens diffraction','fontsize',18);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% subplot(2,2,3)
% imagesc(x1,y1,abs(initial_field))
% title('initial guess object field','fontsize',18);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
% subplot(2,2,4)
% imagesc(x1,y1,abs(initial_diffraction))
% title('initial guess object diffraction','fontsize',18);axis square;colormap('gray');xlabel('x/m');ylabel('y/m');
%% support֧�����С�;���͸����������Ʒ�ϵĴ�Сһ�¡̣��ٰ���ԭ���ķ�����������䡣  
% support = (x_array./ws).^2+(y_array./ws).^2 <= 1; 
%% �������ƶ�
% image1 = ph_p;
% image2 = ph_p_orig;
% glcms1 = graycomatrix(image1);
% glcms2 = graycomatrix(image2);
% similarity = cosin_similarity(glcms1,glcms2);
% cos1 = acos(similarity);
% angle = cos1*180/pi; %�Ƕ�
% fprintf('����ֵ = %f\n',similarity)
% fprintf('���ҽǶ� = %f��\n',angle) 
%% OSS�ؽ���֧����Ϊ200*200�������Σ�rec��2048*2048��Ϊ��Ʒ��������Ľ��
rec_modelimg = kcdi_real; %�������������ģ��
lens_modelimg = kcdi_lens_real; %����������Ʒ���ϵ�͸��ģ��
F2D_rec = am_p; %���������������������
F2D_lens = am_p_lens;%���������͸�����������
supp = [200 200]; %�����֧�����С
[Rsize,Csize] = size(I);
Rcenter = ceil(Rsize/2);
Ccenter = ceil(Csize/2);
Rsupport = supp(1);
Csupport = supp(2);
half_Rsupport = ceil(Rsupport/2);
half_Csupport = ceil(Csupport/2);
iter = 100;
beta = 0.9;
showim = 1;
hiofirst = 0;
%kk��λ������Ҫ�˲���ѡ��
% [mask_sample,RfacR_sample,RfacF_sample,rec]= OSS_samplecode_TF2(F2D_rec,supp,iter,beta,showim,rec_modelimg,hiofirst,L1,lambda,z2,initial_diffraction);  
% [mask_lens,RfacR_lens,RfacF_lens,rec_lens]= OSS_samplecode_lens_TF2(F2D_lens,supp,iter,beta,showim,lens_modelimg,hiofirst,L1,lambda,z2,initial_lens_diffraction);  
%�����λ����Ҫ�˲���ѡ��
[mask_sample,RfacR_sample,RfacF_sample,rec]= OSS_samplecode_TF3(F2D_rec,supp,iter,beta,showim,rec_modelimg,hiofirst,L1,lambda,z2,initial_diffraction);  
[mask_lens,RfacR_lens,RfacF_lens,rec_lens]= OSS_samplecode_lens_TF3(F2D_lens,supp,iter,beta,showim,lens_modelimg,hiofirst,L1,lambda,z2,initial_lens_diffraction);  

%% ��ʾ�ؽ�����Լ��ؽ�ʱ��
reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('OSS: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);
rec_sample_dif = propTF(rec,L1,lambda,z2);
rec_lens_dif = propTF(rec_lens,L1,lambda,z2);
rec_dif = rec_sample_dif - rec_lens_dif;
rec = propTF_inverse(rec_dif,L1,lambda,z2);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./�ؽ����.eps');
rec = rec(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./�Ŵ���ؽ����.eps');