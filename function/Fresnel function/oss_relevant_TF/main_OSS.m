%%
clear
close all
%%
load('10%noise_Diffraction_Pattern.mat')
load('model.mat')
%%
F2D = diffpat; %����������ͼ�������ͼ����ȥ��������ɢ�ߣ�������� 5% �������� (Rnoise)
modelimg = model; %�����������Ĵ�С
supp = [52,28]; %�����֧�����С
iter = 10;
beta = 0.9;
showim = 1;
hiofirst = 0;
[mask RfacR RfacF]= OSS_samplecode (F2D,supp,iter,beta,showim,modelimg,hiofirst)  
% [rec, errK, Rfree_magnitude, bestErrors] = GENFIRE_OSS(numIterations,initialObject,support,measuredK,constraintInd_complex,constraintInd_magnitude,R_freeInd_mag)
