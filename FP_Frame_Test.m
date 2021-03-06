close all;clear;clc;

theta = 0;
wlength = 5.32e-07;
z = 2.5e-05;
xint = 0;
yint = 0;

NA          = 0.39;      % objective NA
spsize      = 1.462e-6; % pixel size of low-res image on sample plane, in m
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in m

xstart = 8; ystart = 8; % absolute coordinate of initial LED
arraysize = 16; % side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);

H      = 60; % distance between LEDs and sample, in mm
LEDp   = 3.4*0.75;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

%[m1, n1, numim] = [3040, 4056, 256];
m1 = 3040;
n1 = 3040;
numim = 256;
pratio = round(spsize/psize); % upsampling ratio
m = pratio*m1; n = pratio*n1;
k0 = 2*pi/wlength;
kx = k0*kx; ky = k0*ky;
NAfilx = NA*(1/wlength)*n*psize; NAfily = NA*(1/wlength)*m*psize; % m1*spsize = m*psize
kmax = pi/psize; % the max wave vector of the OTF
dkx = 2*pi/(psize*n); dky = 2*pi/(psize*m);
kx2 = -kmax:kmax/((n-1)/2):kmax;
ky2 = -kmax:kmax/((m-1)/2):kmax; % odd N


[kxm, kym] = meshgrid(kx2,ky2); 
kzm = sqrt(k0^2-kxm.^2-kym.^2);

H2 = exp(1j.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm))); % define the defocus aberration if it is known or you want to test it
astigx = 0; astigy = 0; % define the astigmatism aberration if it is known or you want to test it
[M1, N1] = meshgrid(1:m1,1:n1);
zn = astigx*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),2,2)+...
     astigy*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),-2,2);
zn = imresize(zn,[m1,n1]);

fmaskproPT1 = 1.*double(((M1-(n1+1)/2)/NAfilx).^2+((N1-(m1+1)/2)/NAfily).^2<=1); % low-pass filter
fmaskproPT2 = H2(round((m+1)/2-(m1-1)/2):round((m+1)/2+(m1-1)/2),round((n+1)/2-(n1-1)/2):round((n+1)/2+(n1-1)/2)); % defocus aberration
for a =1:m1
    for b=1:n1
        if (a<(m1/5))||(a>(m1*4/5))||(b<(n1/5))||(b>(n1*4/5))
            fmaskproPT2(a,b) = 0;
        end
    end
end
fmaskproPT3 = exp(pi*1j.*zn);
fmaskpro = fmaskproPT1.*fmaskproPT2.*fmaskproPT3;
fmaskTest = fmaskpro;

him = zeros(m,n);
himFT = fftshift(fft2(him));
imgArray = zeros(1,numim);

%% main part to optimize estimate of high-res image
for i = 1:2 % 2 initial iterations to get a rough estimate of high-res image
    for i3 = 1:numim         
        kxc=round((n+1)/2-kx(1,i3)/dkx);
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        
%         imshow(him);
%         hold on;
%         rectangle('Position',[kxl,kyl,kxh,kyh],'LineWidth',2);
%         hold off;
%         pause(1)
        O_j=himFT(kyl:kyh,kxl:kxh);
        lowFT=fmaskpro;
        im_lowFT=ifft2(ifftshift(lowFT));
        im_lowFT=1.*exp(1j.*angle(im_lowFT)); 
        lowFT_p=fftshift(fft2(im_lowFT));
        himFT(kyl:kyh,kxl:kxh)=himFT(kyl:kyh,kxl:kxh)+conj(fmaskpro)./(max(max((abs(fmaskpro)).^2))).*(lowFT_p-lowFT);
        
    end
end
imshow(him);