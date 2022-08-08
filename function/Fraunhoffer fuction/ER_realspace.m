function [initial_object]= ER_realspace(diffraction_update,support_pad,vecX,ncX_big,vecY,ncY_big,paddingx,paddingy) 
%衍射图逆傅里叶变换到实空间
initial_object = ifft2(diffraction_update);
%实空间约束
initial_object = initial_object(vecX + ncX_big, vecY + ncY_big);
initial_object(initial_object<0) = 0;   %小于0 的部分为0，正向约束
initial_object = padarray(initial_object,[paddingx paddingy]);
initial_object = initial_object.*support_pad;  %支撑域约束
end




