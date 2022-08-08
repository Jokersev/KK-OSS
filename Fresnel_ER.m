%% KCDI_ER ���͸��֧���򣬵����ؽ���Ϻ���ȥ��͸������,���������TFģ�⣬oss��ʵ��ֵ��ER�����ȶ���Ʒ���й��������KCDI
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

%% ������ʾ
% figure
% subplot(2,2,1);imagesc(I);axis square;colormap('gray');
% subplot(2,2,2);imagesc(ph_p_orig);axis square;
% subplot(2,2,3);imagesc(I_lens);axis square;colormap('gray');
% subplot(2,2,4);imagesc(ph_p_orig_lens);axis square;
%% �������
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% % export_fig(gcf,'-eps','-r300','-painters','./δ�������.eps');
% % figure;imagesc(I_lens);axis square;colormap('gray');title('�������ǰ')
[M,N] = size(I);
lam = 0.1; %���5%=0.05����
r = poissrnd(lam,[M,N]); %��Ӳ�������
I = I+I.*r;
I_lens = I_lens+I_lens.*r;
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./���10%������.eps');
figure;imagesc(I_lens);axis square;colormap('gray');title('͸�����������')

%% ���ò���
L1=28672e-6; 
M=2048;       %number of samples
lambda=632e-9;     %wavelength(m)
z1=0.601;             %ͨ��͸�������������601mm������Ʒ
z2=0.463;             %��ccd�ķ������������463mm
zf=0.474;           %������۵ľ��루���ࣩ474mm
z=0.474+0.59;
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
%����͸��������CCD�����䳡
ulens_2=propTF(ulens_1,L1,lambda,z2); 
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
% [initial_lens_field,initial_lens_diffraction,ws] = initial_guess_lens(L1,lambda,ph_p_lens,am_p_lens,z2);
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
%% ��ʼ����������֧����
size_pixel = size(I);
N = size_pixel(1);
oversamplingRatio = 3;
newDimx = size(I,1)*oversamplingRatio;%size of oversampled array
newDimy = size(I,2)*oversamplingRatio;%size of oversampled array
paddingx = round((newDimx-size(I,1))/2);%how many zeros to add
paddingy = round((newDimy-size(I,2))/2);%how many zeros to add
initial_lens_field_pad = padarray(initial_lens_field,[paddingx paddingy]);
initial_field_pad = padarray(initial_field,[paddingx paddingy]);
%% support֧�����С�;���͸����������Ʒ�ϵĴ�Сһ�¡̣��ٰ���ԭ���ķ�����������䡣  
support = (x_array./ws).^2+(y_array./ws).^2 <= 1;  
support_pad = padarray(support,[paddingx paddingy]);
initial_lens_field = initial_lens_field_pad.*support_pad;
initial_field = initial_field_pad.*support_pad;
%͸����֧����
supportlens_pad = padarray(u1,[paddingx paddingy]);  
%% ��ʼֵ��ʾ
ncX_big = round((size(initial_field,1)+1)/2);
ncY_big = round((size(initial_field,2)+1)/2);
vecX = 1:size(I,1); ncX = round((size(I,1)+1)/2); vecX = vecX - ncX;
vecY = 1:size(I,2); ncY = round((size(I,2)+1)/2); vecY = vecY - ncY;
initial_guess_field = initial_field(vecX + ncX_big, vecY + ncY_big);
initial_guess_lens_field = initial_lens_field(vecX + ncX_big, vecY + ncY_big);
%% ER�����ؽ�
Iterations = 100;
modelimg = kcdi_real; %�������������ģ��
lens_modelimg = kcdi_lens_real; %����������Ʒ���ϵ�͸��ģ��

[rec,rec_lens] = ER_reconstruction2(x1,y1,Iterations,initial_field,initial_diffraction,initial_lens_field,initial_lens_diffraction,support_pad,supportlens_pad,paddingx,paddingy,L1,lambda,z2,z,am_p,am_p_lens,vecX,ncX_big,vecY,ncY_big,ph_p_orig,ph_p_orig_lens,lens_modelimg,modelimg);

%% er�ؽ����
supp=[200 200];
[Rsize,Csize] = size(I);
Rcenter = ceil(Rsize/2);
Ccenter = ceil(Csize/2);
Rsupport = supp(1);
Csupport = supp(2);
half_Rsupport = ceil(Rsupport/2);
half_Csupport = ceil(Csupport/2);
reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('ER: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./�ؽ����.eps');
rec = rec(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./�Ŵ���ؽ����.eps');
%% �ؽ�����λ��ԭͼ��λ�ıȽ�,������ʱ�ɲ���
% figure
% subplot(1,2,1);plot(x1,unwrap(ph_p_orig(M/2+10,:)));xlabel('x(m)');ylabel('y(rad)');axis square;title('ccd free noise E(r) angle M/2+10')
% subplot(1,2,2);plot(x1,unwrap(phase(rec_sample_dif(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;title('rec E(r) angle M/2+10')