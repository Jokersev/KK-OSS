%% KCDI模拟实验，TF模拟
close all
clear
%% 单位全是米 改变CCD的采样边长、改变相位求法
% L1=2867e-6;      %side length(m)  28672微米=28.672mm失败
L1=28672e-6; 
M=2048;       %number of samples
dx1=L1/M;    %src sample interval
x1=-L1/2:dx1:L1/2-dx1;    %src coords
y1=x1;
lambda=632e-9;     %wavelength(m)
k=2*pi/lambda;      %wavenumber
z1=0.601;             %从透镜出射后菲涅尔传输的距离474+127mm
z2=0.463;             %到ccd的菲涅尔传输距离463mm
zf=0.474;           %光束汇聚的距离（焦距）474mm
wl=400;              %透镜镜头口径半径(像素为单位)
% wl=400;    %常用值透镜口径(像素为单位)
[X1,Y1]=meshgrid(x1,x1);
% u1=circ(sqrt(X1.^2+Y1.^2)/wl);   %圆形波
[x_array,y_array] = meshgrid(1:M,1:M); 
x_array = x_array - floor(max(x_array(:))/2+1); % center of image to be zero 
y_array = y_array - floor(max(y_array(:))/2+1); % center of image to be zero 
u1=(x_array./wl).^2+(y_array./wl).^2 <= 1; 
% figure;imagesc(u1);axis square;title('透镜大小');

%% 透镜出射光场
uout=u1.*exp(-1i*k/(2*zf)*(X1.^2+Y1.^2));
I1=(abs(u1)).^2;
% figure;imagesc(x1,y1,I1);colormap('gray');axis square;title('透镜前')
% figure;imagesc(x1,y1,angle(u1));colormap('gray');axis square;title('通过透镜前的相位');

Iout=(abs(uout)).^2;
figure;imagesc(x1,y1,Iout);colormap('gray');title('通过透镜后');axis square
figure;imagesc(x1,y1,angle(uout));colormap('gray');axis square;title('通过透镜后的相位');

u2=propTF(u1,L1,lambda,z1); %propagation, TF情况下u2的边长L2=L1
I2=abs(u2.^2);             %obs irrad
x2=x1;
y2=y1;
% figure;
%imagesc(x2,y2,I2);axis xy;axis square;colormap('gray');xlabel('x(m)');ylabel('y(m)');title(['没有通过透镜z=',num2str(z1),'m']);
%% 通过透镜后菲涅尔传输
%透镜出射后的菲涅尔传输propagation, TF情况下u2的边长L2=L1
u3=propTF(uout,L1,lambda,z1); 
% u3=propIR(u1,L1,lambda,z1);   %去掉透镜的作用
I3=abs(u3.^2);             %obs irrad
x3=x1;
y3=y1;
% figure;
%imagesc(x2,y2,I3);axis xy;axis square;colormap('gray');xlabel('x(m)');ylabel('y(m)');title(['菲涅尔衍射z=',num2str(z1),'m','到达样品前']);
%% 加入样品
obj1=imread('USAF1951B250_2048','png');   %read image file
obj1=rgb2gray(obj1);
obj=obj1;
Iobj=single(obj);                     %integer to floating 整数转化为单精度浮点数
Iobj=Iobj/max(max(Iobj));               %normalize ideal image
uobj=sqrt(Iobj);                      %ideal image field
xobj=x1;yobj=y1;
figure%原图
imagesc(xobj,yobj,Iobj);title('原图');colormap('gray');xlabel('u(m)');ylabel('v(m)');axis square;
uobj=circshift(uobj,-100,1);
%% 平面波通过样品，复振幅乘以样品的透射率
u4=u3.*uobj;  
% u4 = u3;   %没有添加样品，只进行透镜的菲涅尔衍射
x4=x1;
y4=y1;
I4=(abs(u4)).^2;
figure('color',[1 1 1]);imagesc(x4,y4,I4);colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');%title('通过样品后的光强');
figure;imagesc(x4,y4,(angle(u4)),[-pi,pi]);axis square;title('通过样品后相位');
%% 菲涅尔传输到达ccd
L5=L1;
M=2048;       %number of samples
dx5=L5/M;    %src sample interval
x5=-L5/2:dx5:L5/2-dx5;    %src coords
y5=x5;
u5=propTF(u4,L5,lambda,z2); %原本是TF
% u5=propIR(u4,L5,lambda,z2); 
I5=(abs(u5)).^2;
% figure;imagesc(x5,y5,I5);title('CCD intensity');colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');
figure('color',[1 1 1]);imagesc(x5,y5,I5);colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');
figure('color',[1 1 1]);imagesc(x5,y5,angle(u5));axis square;xlabel('x(m)');ylabel('y(m)');%title('CCD phase');
% figure;imagesc(x5,y5,I5);title('CCD');axis square
% figure;imagesc(x5,y5,angle(u5));xlabel('x(m)');ylabel('y(m)');title('CCD phase');axis square
figure('color',[1 1 1]);plot(x5,unwrap(angle(u5(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;%title('CCD phase M/2+10');

%% 保存数据
% Intensity = I5;
% save('KCDI_CCDsample_intensity_TF.mat','Intensity')
% lens_Intensity = I5;
% save('KCDI_CCDlens_intensity_TF.mat','lens_Intensity')
% % Angle = angle(u5);
% % save('KCDI_CCDsample_angle_TF.mat','Angle')
% lens_Angle = angle(u5);
% save('KCDI_CCDlens_angle_TF.mat','lens_Angle')
% % % save('CCD_field.mat','u5')
%% 物体单独衍射
% uobj_dif=propTF(uobj,L5,lambda,z2);
% figure
% imagesc(x5,y5,abs(uobj_dif));
% xlabel('x(m)');ylabel('y(rad)');
% title('物体单独衍射')
% axis square;colormap('gray')
% Angle_uobj = angle(uobj_dif);
% save('KCDI_CCD_Angle_uobj.mat','Angle_uobj')
%% inverse IR
% u42=propIR_inverse(u5,L5,lambda,z2); 
% figure;imagesc(x5,y5,abs(u42));xlabel('x(m)');ylabel('y(m)');title('IR逆变换');axis square
% figure;imagesc(x5,y5,abs(u5));xlabel('x(m)');ylabel('y(m)');title('IR正变换');axis square
%% kk关系
ws = 100;
ulens_2=propTF(u3,L1,lambda,z2); 
kr=angle(ulens_2);
kk_imaginary_part = kk2(I5,L1,lambda,M);
initial_diffraction=exp((log(I5)./2-1i*kk_imaginary_part)+1i*kr);   %kr写作透镜经过菲涅尔衍射之后的相位
ph_p_kk=imag(log(initial_diffraction));
figure('color',[1 1 1]);imagesc(x5,y5,angle(initial_diffraction));axis square;xlabel('x(m)');ylabel('y(m)');%title('CCD phase');
figure('color',[1 1 1]);plot(x5,unwrap(angle(initial_diffraction(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;%title('CCD phase M/2+10');
figure
subplot(1,2,1);imagesc(x1,y1,ph_p_kk);xlabel('x(m)');ylabel('y(m)');title('kk E(r) angle','fontsize',22);axis square
subplot(1,2,2);plot(x1,unwrap(phase(initial_diffraction(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');title('kk E(r) angle M/2+10','fontsize',22);axis square