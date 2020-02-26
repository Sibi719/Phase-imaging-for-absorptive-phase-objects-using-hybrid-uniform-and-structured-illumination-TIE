clc;
clear all;
close all;
Wavelength=850*10^-9;
N=2^10;
k=(2*pi)/Wavelength;
R1=8*10^-3;
R2=5*10^-3;
d=0.005;
Dimesnion=20*10^-3;
del1=Dimesnion/N;
alpha1=0.002%0.0000002;
beta1=0.0001;%0.00000001;
alpha2=alpha1/1.1;
beta2=beta1/2;
cx1=0;
cy1=0;
wn=0;
cx2=-R1/2;
cy2=0*10^-3;
M=zeros(N);
x=((-N/2):((N/2)-1))*del1;
y=((-N/2):((N/2)-1))*del1;
du = 1/(N*del1);
umax = 1/(2*del1);
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
fr_sq=U.^2+V.^2;
cut=0;
n=0.5;
noise=randn(N).*sqrt(10^(-15/10));
[X,Y]=meshgrid(x,y);
L1= real(sqrt(R1^2-((X-cx1).^2 + (Y-cy1).^2)));
L2= real(sqrt(R2^2-((X-cx2).^2 + (Y-cy2).^2)));
phase1 = k*alpha1*L1/1.4;
phase2 = k*alpha2*L2;
phase= phase1+phase2;
T1 = k*beta1*L1/8;
T2 = k*beta2*L2;

T= exp(-(T1+T2)).^2;

figure
imagesc(x,y,phase);
axis image
title(['Phase profile ']);
colorbar
colormap(jet)

figure
imagesc(x,y,T);
axis image
title(['Transmission profile']);
colorbar
colormap(jet)


[phase_dx,phase_dy]= gradient(flipud(phase),del1,del1);
figure
sv = 50;   
sc=2
hold on% ‘Step’ Value
h=quiver(X(1:sv:end,1:sv:end), Y(1:sv:end,1:sv:end), phase_dx(1:sv:end,1:sv:end), phase_dy(1:sv:end,1:sv:end))
set(h,'AutoScale','on', 'AutoScaleFactor', sc)
set(h,'LineWidth',1)

[T_dx,T_dy]= gradient(flipud(T),del1,del1);
h=quiver(X(1:sv:end,1:sv:end), Y(1:sv:end,1:sv:end), -T_dx(1:sv:end,1:sv:end), -T_dy(1:sv:end,1:sv:end),'r')
set(h,'AutoScale','on', 'AutoScaleFactor', sc)
set(h,'LineWidth',1)
hold off

Uo=ones(N);
Hcrop =sqrt(Uo.*T).*exp(1j.*phase);
[phase_u,deld,Io]=uniform(del1,d,Wavelength,N,k,fr_sq,noise,Hcrop,Uo,x,y,wn);

%%%Structured illumination

function [phase_u, Intensity_diff,Io]=uniform(del1,d,Wavelength,N,k,fr_sq,noise,Hcrop,Uo,x,y,wn)

z=0;
[Io]=FresProp(del1,z,Wavelength,N,Hcrop);
% figure
% imagesc(x,y,Io);
% axis image
% title(['Intensity pattern at z= ',num2str(z),'m']);
% colorbar
% colormap(jet)
if wn==1
noise=randn(N).*sqrt(10^(-25/10));
Io = Io +noise;
end


z=d;
[IP]=FresProp(del1,z,Wavelength,N,Hcrop);
if wn==1
noise=randn(N).*sqrt(10^(-25/10));
IP = IP +noise;
end
% figure
% imagesc(x,y,IP);
% axis image
% title(['Intensity pattern at z= ',num2str(z),'m']);
% colorbar
% colormap(jet)

z=-d;
[IM]=FresProp(del1,z,Wavelength,N,Hcrop);
if wn==1
noise=randn(N).*sqrt(10^(-25/10));
IM = IM +noise;
end
% figure
% imagesc(x,y,IM);
% axis image
% title(['Intensity pattern at z= ',num2str(z),'m']);
% colorbar
% colormap(jet)

%%Intensity differential

H= (-4.*pi.*pi.*fr_sq);
HH=H./(H.^2 + eps^2);

Intensity_diff=((IP-Io))./(d);

D= (ifft2(ifftshift(HH.*(fftshift(fft2(Intensity_diff))))));
% figure
% imagesc(x,y,D);
% title("D");
% axis image
% colorbar
% colormap(jet)

[D_dx,D_dy]= gradient(D,del1,del1);

E1=(1./Io).*D_dx;
E2=(1./Io).*D_dy;
[E1_dx,E1_dy]= gradient(E1,del1,del1);
[E2_dx,E2_dy]= gradient(E2,del1,del1);

% figure
% imagesc(x,y,E1_dx+E2_dy);
% title("E1_dx+E2_dy");
% axis image
% colorbar
% colormap(jet)

phase_u= -k.*ifft2(ifftshift(HH.*(fftshift(fft2(E1_dx+E2_dy)))));

phase_u= phase_u + abs(min(min(phase_u)));

figure
imagesc(x,y,phase_u);
title("phase_u");
axis image
colorbar
colormap(jet)
end

function [I]=FresProp(dpix,z,lambda,Hsize,Hcrop)
 
%Spatial frequencies
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Nyquist cut-off for Sampling Hologram
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
