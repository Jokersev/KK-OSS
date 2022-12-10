%% KCDI_ER迭代重建完毕后再去掉透镜作用,夫朗禾费衍射傅里叶变换，oss的实验值和ER都是先对样品进行过采样后的KCDI
%% 
close all
clear
%% 导入KCDI实验数据
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDsample_intensity_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDsample_angle_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDlens_intensity_IR.mat')
% load('F:\users\wy\KCDI_ER\IR data\KCDI_CCDlens_angle_IR.mat')
% load('F:\users\wy\KCDI_ER\oss data\KCDI_sample_supportfield.mat')
% load('F:\users\wy\KCDI_ER\oss data\KCDI_lens_supportfield.mat')
load('KCDI_CCD_angle_Fraunhoffer.mat')
load('KCDI_CCD_intensity_Fraunhoffer.mat')
load('KCDI_CCDlens_angle_Fraunhoffer.mat')
load('KCDI_CCDlens_intensity_Fraunhoffer.mat')
load('KCDI_lens_supportfield_Fraunhoffer.mat')
load('KCDI_sample_supportfield_Fraunhoffer.mat')
I=Intensity;
ph_p_orig = Angle;
I_lens = lens_Intensity;
ph_p_orig_lens = lens_Angle;

%% 显示原数据
% figure
% subplot(2,2,1);imagesc(I);axis square;colormap('gray');
% subplot(2,2,2);imagesc(ph_p_orig);axis square;
% subplot(2,2,3);imagesc(I_lens);axis square;colormap('gray');
% subplot(2,2,4);imagesc(ph_p_orig_lens);axis square;
%% 添加噪声
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./未添加噪声.eps');
% figure;imagesc(I_lens);axis square;colormap('gray');title('添加噪声前')
[M,N] = size(I);
lam = 0.1; %添加5%=0.05噪声
r = poissrnd(lam,[M,N]); %添加泊松噪声
I = I+I.*r;
I_lens = I_lens+I_lens.*r;
figure('color',[1 1 1]);imagesc(I);axis square;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./添加10%噪声后.eps');
% figure;imagesc(I_lens);axis square;colormap('gray');title('透镜添加噪声后')


%% 设置参数
L1=28672e-6; 
M=2048;       %number of samples
lambda=632e-9;     %wavelength(m)
z1=0.601;             %通过透镜后菲涅尔传播601mm到达样品
z2=0.463;             %到ccd的菲涅尔传输距离463mm
zf=0.474;           %光束汇聚的距离（焦距）474mm
wl=400;    %常用值透镜口径(像素为单位)
%%
dx1=L1/M;    %src sample interval
x1=-L1/2:dx1:L1/2-dx1;    %src coords
y1=x1;
[X1,Y1]=meshgrid(x1,y1);
k=2*pi/lambda;      %wavenumber
[x_array,y_array] = meshgrid(1:M,1:M); 
x_array = x_array - floor(max(x_array(:))/2+1); % center of image to be zero 
y_array = y_array - floor(max(y_array(:))/2+1); % center of image to be zero 
%计算透镜传播z1距离的衍射场
u1=(x_array./wl).^2+(y_array./wl).^2 <= 1; 
uout=u1.*exp(-1i*k/(2*zf)*(X1.^2+Y1.^2));
ulens_1=propTF(uout,L1,lambda,z1); 
%计算透镜夫朗禾费传播到CCD的衍射场
ulens_2=fftshift(fft2(ulens_1)); 
% figure
% imagesc(x1,y1,I)
% title('CCD intensity')
% axis square;xlabel('x/m');ylabel('y/m');colormap('gray')
%% ER算法
%% 初始物体估计值
siza = size(I);                                     
am_p = sqrt(I);
am_p_lens = sqrt(I_lens);
%透镜的初始相位
ph_p_lens = rand(siza)*2*pi;% rand(N(1),N(1))*pi; %随机相位
% ph_p_lens = ph_p_orig_lens; %透镜的原图相位
% ph_p_lens =imnoise(ph_p_orig_lens,'salt & Pepper');  %原图相位加椒盐噪声
% ph_p_lens =imnoise(ph_p_orig_lens,'gaussian',0.4);  %原图相位加高斯噪声 
% ph_p_lens =imnoise(ph_p_orig_lens,'poisson');  %原图相位加泊松噪声 
%物体的初始相位
ph_p = rand(siza)*2*pi;% rand(N(1),N(1))*pi; %随机相位
% ph_p = ph_p_orig; %物体的原图相位
% ph_p = imnoise(ph_p_orig,'salt & Pepper');  %原图相位加椒盐噪声
% ph_p = imnoise(ph_p_orig,'gaussian',0.4);  %原图相位加高斯噪声 
% ph_p = imnoise(ph_p_orig,'poisson');  %原图相位加泊松噪声 
if_kk_lens = 0;
if_kk = 0;
tic;
[initial_lens_field,initial_lens_diffraction,ws] = initial_guess_Fraunhoffer(x1,y1,if_kk_lens ,ulens_2,I_lens,L1,lambda,M,ph_p_lens,am_p_lens,siza,z2);
[initial_field,initial_diffraction,ws] = initial_guess_Fraunhoffer(x1,y1,if_kk ,ulens_2,I,L1,lambda,M,ph_p,am_p,siza,z2);
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
% figure('color',[1 1 1]);plot(x1,unwrap(ph_p_orig(M/2+10,:)));xlabel('x(m)');ylabel('y(rad)');axis square;
% export_fig(gcf,'-eps','-r300','-painters','./原图相位.eps');
% figure('color',[1 1 1]);plot(x1,unwrap(phase(initial_diffraction(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;
% export_fig(gcf,'-eps','-r300','-painters','./kk相位.eps');

%% 计算相似度
% image1 = ph_p;
% image2 = ph_p_orig;
% glcms1 = graycomatrix(image1);
% glcms2 = graycomatrix(image2);
% similarity = cosin_similarity(glcms1,glcms2);
% cos1 = acos(similarity);
% angle = cos1*180/pi; %角度
% fprintf('余弦值 = %f\n',similarity)
% fprintf('余弦角度 = %f°\n',angle) 
%% 初始物体零填充和支持域
size_pixel = size(I);
N = size_pixel(1);
oversamplingRatio = 3;
newDimx = size(I,1)*oversamplingRatio;%size of oversampled array
newDimy = size(I,2)*oversamplingRatio;%size of oversampled array
paddingx = round((newDimx-size(I,1))/2);%how many zeros to add
paddingy = round((newDimy-size(I,2))/2);%how many zeros to add
initial_lens_field_pad = padarray(initial_lens_field,[paddingx paddingy]);
initial_field_pad = padarray(initial_field,[paddingx paddingy]);

%% OSS重建，支持域为200*200的正方形，rec（2048*2048）为样品过采样后的结果
%% er重建结果在实空间进行过采样，进行衍射图相加
rec_modelimg = kcdi_real; %过采样后的物体模型
lens_modelimg = kcdi_lens_real; %过采样后样品面上的透镜模型
% rec_modelimg_diffraction = propTF(kcdi_real,L1,lambda,z2);
% lens_modelimg_diffraction = propTF(kcdi_lens_real,L1,lambda,z2);
% F2D_rec = abs(rec_modelimg_diffraction); %过采样后的物体的衍射幅度
% F2D_lens = abs(lens_modelimg_diffraction); %过采样后的透镜的衍射幅度
F2D_rec = am_p; %过采样后的物体的衍射幅度
F2D_lens = am_p_lens;%过采样后的透镜的衍射幅度
supp = [200,200]; %物体的支持域大小
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
%%
[mask_sample,RfacR_sample,RfacF_sample,rec]= OSS_samplecode4(F2D_rec,supp,iter,beta,showim,rec_modelimg,hiofirst,initial_diffraction);  
[mask_lens,RfacR_lens,RfacF_lens,rec_lens]= OSS_lenscode4(F2D_lens,supp,iter,beta,showim,lens_modelimg,hiofirst,initial_lens_field);  
%重建时间
reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('Fraunhoffer OSS: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);
rec_sample_dif = fft2(rec);
rec_lens_dif = fft2(rec_lens);
rec_dif = rec_sample_dif - rec_lens_dif;
rec = ifft2(rec_dif);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./重建结果.eps');
rec = rec(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
figure('color',[1 1 1]);imagesc(abs(rec));axis square;axis off;colormap('gray');
% export_fig(gcf,'-eps','-r300','-painters','./放大的重建结果.eps');
%% 重建的相位和原图相位的比较,带噪声时可不开
% figure
% subplot(1,2,1);plot(x1,unwrap(ph_p_orig(M/2+10,:)));xlabel('x(m)');ylabel('y(rad)');axis square;title('ccd free noise E(r) angle M/2+10')
% subplot(1,2,2);plot(x1,unwrap(phase(rec_sample_dif(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;title('rec E(r) angle M/2+10')
