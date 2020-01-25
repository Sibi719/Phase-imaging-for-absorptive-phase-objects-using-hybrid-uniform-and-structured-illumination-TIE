clc;
clear all;
close all;
lambda=633*10^-9;
Hsize=2^10;
k=(2*pi)/lambda;
R1=4.3*10^-3;
R2=3.4*10^-3;
z=0.4;
delta=0.001;
beta=0.08;
cx1=0;
cy1=0;
cx2=-R1/2;
cy2=0;
Dimesnion=8*10^-3;
dpix=Dimesnion/Hsize;
x=((-Hsize/2):((Hsize/2)-1))*dpix;
y=((-Hsize/2):((Hsize/2)-1))*dpix;
[X,Y]=meshgrid(x,y);
umax = 1/(2*dpix); 
df=1/(Hsize*dpix);
u=(-Hsize/2:Hsize/2-1)*df;
v=(-Hsize/2:Hsize/2-1)*df;
[X,Y]=meshgrid(x,y);
[U,V]=meshgrid(u,v);
spat_sq = X.^2+Y.^2;
freq_sq= U.^2+V.^2;

L1= real(sqrt(R^2-((X-cx).^2 + (Y-cy).^2)));
L2= real(sqrt(R^2-((X-cx).^2 + (Y-cy).^2)));

phase = k*delta*L;

figure
imagesc(x,y,phase);
axis image
title(['Phase profile ']);
colorbar
colormap(jet)

T = exp(-( 70.*(X+max(max(X)))));

figure
imagesc(x,y,T);
axis image
title(['T']);
colorbar
colormap(jet)

%%%Uniform Illumination


Uo=ones(Hsize);
Hcrop=sqrt(Uo).*T.*exp(1j.*phase);

I_d=FresProp(dpix,z,lambda,Hsize,Hcrop);


deld= I_d - ((Uo).*T);

figure
imagesc(x,y,deld);
axis image
title(['del_d without noise']);
colorbar
colormap(jet)


SNR=0.00000000000000001;
g= deld - min(deld(:));
g =g / max(g(:));
Var_bi= var((g (:)));
Var_noise = Var_bi./(10^(SNR/10))
deldn=deld + (sqrt(10^-2.5).*randn(size(deld)));

figure
imagesc(x,y,deldn);
axis image
title(['del_d with noise ']);
colorbar
colormap(jet)

figure
imagesc(x,y,deldn-deld);
axis image
title(['noise']);
colorbar
colormap(jet)

ft_in= (-k./(z.*Uo)).*deld;
T1=fftshift((fft2(ft_in)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
D= ifft2(ifftshift(T2.*T1));

% figure
% imagesc(x,y,D);
% axis image
% title(['D']);
% colorbar
% colormap(jet)

%%%phase equation

[del_invTx,del_invTy] = gradient((1./(T)),dpix,dpix);
[Dx,Dy] = gradient((D),dpix,dpix);

T1=fftshift(fft2((del_invTx.*Dx) + (del_invTy.*Dy) +  (4*del2(D,dpix,dpix)./T)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
phasee= ifft2(ifftshift(T2.*T1));

% figure
% imagesc(x,y,phasee);
% axis image
% title(['phase']);
% colorbar
% colormap(jet)

T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
AA= ifft2(ifftshift(T2.*fftshift(fft2(-k.*deld./z))));

[del_invTx,del_invTy] = gradient((1./(T)),dpix,dpix);
[Dx,Dy] = gradient((AA),dpix,dpix);

T1=fftshift((fft2(del_invTx.*Dx + del_invTy.*Dy +  4*del2(AA,dpix,dpix)./T)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
phasee= ifft2(ifftshift(T2.*T1));


figure
imagesc(x,y,phasee);
axis image
title(['Recovered phase wihtout noise']);
colorbar
colormap(jet)
%caxis([min(min(phase)), max(max(phase))]);


%%%%
ft_in= (-k./(z.*Uo)).*deldn;
T1=fftshift((fft2(ft_in)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
D= ifft2(ifftshift(T2.*T1));

% figure
% imagesc(x,y,D);
% axis image
% title(['D']);
% colorbar
% colormap(jet)

%%%phase equation

[del_invTx,del_invTy] = gradient((1./(T)),dpix,dpix);
[Dx,Dy] = gradient((D),dpix,dpix);

T1=fftshift(fft2((del_invTx.*Dx) + (del_invTy.*Dy) +  (4*del2(D,dpix,dpix)./T)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
phasee= ifft2(ifftshift(T2.*T1));

% figure
% imagesc(x,y,phasee);
% axis image
% title(['phase']);
% colorbar
% colormap(jet)

T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
AA= ifft2(ifftshift(T2.*fftshift(fft2(-k.*deldn./z))));

[del_invTx,del_invTy] = gradient((1./(T)),dpix,dpix);
[Dx,Dy] = gradient((AA),dpix,dpix);

T1=fftshift((fft2(del_invTx.*Dx + del_invTy.*Dy +  4*del2(AA,dpix,dpix)./T)));
T2= (-4.*pi.*pi.*freq_sq);
T2=T2./(T2.^2 + eps^2);
phasee= ifft2(ifftshift(T2.*T1));


figure
imagesc(x,y,phasee);
axis image
title(['Recovered phase with noise']);
colorbar
colormap(jet)
%caxis([min(min(phase)), max(max(phase))]);

function [I]=FresProp(dpix,z,lambda,Hsize,Hcrop)
 
%Spatial frequeHsizecies
Xsize = Hsize*dpix; %Hsize is the Hsizeumber of pixel, dpix is the leHsizegth of oHsizee pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Hsizeyquist cut-off for Sampling Hologram
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
 
%Evanescent cut-off 
uev = 1/lambda; %???????
 
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
 
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
 
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
 
%Truncate kernel
H = W.*H;
clear W;
 
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
 
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
R = ifft2(RR);
I=abs(R).^2;
clear RR;
end

