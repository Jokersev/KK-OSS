%%
clear
close all
%%
load('10%noise_Diffraction_Pattern.mat')
load('model.mat')
%%
F2D = diffpat; %过采样囊泡图像的衍射图案，去除了中央散斑，并添加了 5% 泊松噪声 (Rnoise)
modelimg = model; %物体过采样后的大小
supp = [52,28]; %物体的支持域大小
iter = 10;
beta = 0.9;
showim = 1;
hiofirst = 0;
[mask RfacR RfacF]= OSS_samplecode (F2D,supp,iter,beta,showim,modelimg,hiofirst)  
% [rec, errK, Rfree_magnitude, bestErrors] = GENFIRE_OSS(numIterations,initialObject,support,measuredK,constraintInd_complex,constraintInd_magnitude,R_freeInd_mag)
