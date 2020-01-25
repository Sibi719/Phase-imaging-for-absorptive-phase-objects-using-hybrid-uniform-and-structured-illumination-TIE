clc;
clear all;
close all
Wavelength=0.0633*10^-6;
N=2^10;
k=(2*pi)/Wavelength;
R1=8*10^-3;
R2=5*10^-3;
Dimesnion=20*10^-3;
del1=Dimesnion/N;
alpha1=0.0002%0.0000002;
beta1=0.00001;%0.00000001;
alpha2=alpha1/1.5;
beta2=beta1/2;
cx1=0;
cy1=0;
cx2=-4*10^-3;
cy2=2*10^-3;
M=zeros(N);
x=((-N/2):((N/2)-1))*del1;
y=((-N/2):((N/2)-1))*del1;
[X,Y]=meshgrid(x,y);
L1= real(sqrt(R1^2-((X-cx1).^2 + (Y-cy1).^2)));
L2= real(sqrt(R2^2-((X-cx2).^2 + (Y-cy2).^2)));
phase1 = k*alpha1*L1/2;
phase2 = k*alpha2*L2;
phase= phase1+phase2;
T1 = k*beta1*L1/8;
T2 = k*beta2*L2;

T= exp(-(T1+T2));

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

% figure
% xi = [0 size(phase,2)];
% yi = [size(phase,1)/2 size(phase,1)/2  ];
% c = improfile(phase,xi,yi);
% plot(x*10^3, c((2:N+1),1));
% xlabel('x (mm)') 
% ylabel('Intensity') 
% axis on
% title(['Cross sectional phase profile']);
% 
% figure
% xi = [0 size(T,2)];
% yi = [size(T,1)/2 size(T,1)/2  ];
% c = improfile(T,xi,yi);
% plot(x*10^3, c((2:N+1),1));
% xlabel('x (mm)') 
% ylabel('Intensity') 
% axis on
% title(['Cross sectional transmission profile']);

[phase_dx,phase_dy]= gradient(phase,del1,del1);
figure
sv = 50;   
sc=2
hold on% ‘Step’ Value
h=quiver(X(1:sv:end,1:sv:end), Y(1:sv:end,1:sv:end), phase_dx(1:sv:end,1:sv:end), phase_dy(1:sv:end,1:sv:end))
set(h,'AutoScale','on', 'AutoScaleFactor', sc)
set(h,'LineWidth',1)

[T_dx,T_dy]= gradient(T,del1,del1);
h=quiver(X(1:sv:end,1:sv:end), Y(1:sv:end,1:sv:end), -T_dx(1:sv:end,1:sv:end), -T_dy(1:sv:end,1:sv:end),'r')
set(h,'AutoScale','on', 'AutoScaleFactor', sc)
set(h,'LineWidth',1)
hold off

U_o=ones(N);
z=0.01;
Hcrop =sqrt(U_o).*T.*exp(1j.*phase);
[IP]=FresProp(del1,z,Wavelength,N,Hcrop);
figure
imagesc(x,y,IP);
axis image
title(['Intensity pattern at z= ',num2str(z),'m']);
colorbar
colormap(jet)

z=-z;
[IM]=FresProp(del1,z,Wavelength,N,Hcrop);
figure
imagesc(x,y,IM);
axis image
title(['Intensity pattern at z= ',num2str(z),'m']);
colorbar
colormap(jet)

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