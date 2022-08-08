function [initial_field,initial_diffraction,ws] = initial_guess_Fraunhoffer(x1,y1,if_kk ,ulens_2,I,L1,lambda,M,ph_p,am_p,siza,z2)
if if_kk == 1
%     ws = 80;
    ws = 100;
    kr=angle(ulens_2);
    kk_imaginary_part = kk2(I,L1,lambda,M);
    initial_diffraction=exp((log(I)./2-1i*kk_imaginary_part)+1i*kr);   %kr写作透镜经过菲涅尔衍射之后的相位
    ph_p=imag(log(initial_diffraction));
    figure
    imagesc(x1,y1,ph_p)
    xlabel('x(m)');ylabel('y(m)');
    title('kk E(r) angle')
    axis square
    figure
    plot(x1,unwrap(phase(initial_diffraction(M/2+10,:))));
    xlabel('x(m)');ylabel('y(rad)');
    title('kk E(r) angle M/2+10')
    axis square
%     ph_p = angle(initial_diffraction);
else
    initial_diffraction = am_p.*exp(j*ph_p); %初始光场形式
    ws = 100;
end
% initial_field = propTF_inverse(initial_diffraction,L1,lambda,z2);
initial_field = ifft2(initial_diffraction);
end